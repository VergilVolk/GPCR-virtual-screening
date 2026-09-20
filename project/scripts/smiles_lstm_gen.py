"""
方向2: SMILES LSTM 迁移学习生成器 (REINVENT 式, 纯 CPU)
=======================================================
从已知 PAM 阳性分子学习结构分布, 生成类似物:

  预训练  在 PubChem 类药语料上训字符级 LSTM -> 学「药物分子语法」(语言先验)
  fine-tune  用少量 PAM 阳性分子微调 -> 把分布掰向 PAM 化学型
  采样     自回归采样新 SMILES -> RDKit 校验/规范化 -> 类药+药效团过滤

这是「极少阳性 -> 生成类似物」的业界标准范式 (Olivecrona/Segler REINVENT),
比 diffusion 更适合: 纯 CPU、专为 few-shot 设计、不会因 N=2 直接记忆。

子命令:
  pretrain   下载语料上预训练, 存 model.pt + vocab.json
  finetune   加载预训练, 用阳性分子微调, 存 model_ft.pt
  generate   采样生成候选, 过滤后存 CSV

用法: python scripts/smiles_lstm_gen.py pretrain|finetune|generate
"""
from __future__ import annotations
import sys
import json
import math
import random
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Crippen, Lipinski
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

GEN = Path("data/generated")
MODELS = Path("results/generated")
CORPUS = GEN / "corpus.smi"
ACTIVES = GEN / "actives.smi"

# 已知 CHRM4 (M4) PAM 阳性 (优先读 data/generated/actives.smi, 否则用此默认)
_DEFAULT_ACTIVES = [
    "Cc1nnc2sc(C(=O)NC3CN(c4ccncc4Cl)C3)c(N)c2c1C",
    "Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)NC3Cc4ccccc4C3)c(N)c12",
    "Cc1c(Cl)nnc2sc(C(=O)NC3Cc4ccccc4C3)c(N)c12",
    "COc1cc(N2CC(NC(=O)c3sc4nnc(Cl)c(C)c4c3N)C2)c(Cl)cn1",
    "Cc1nnc2sc(C(=O)NC3CN(c4cc(F)nc(F)c4)C3)c(N)c2c1C",
    "Cc1c(Cl)c2nncn2c2sc(C(=O)NC3Cc4ccccc4C3)c(N)c12",
]

SPECIAL = ["<PAD>", "<BOS>", "<EOS>", "<UNK>"]
PAD, BOS, EOS, UNK = 0, 1, 2, 3
MAX_LEN = 120

# 类药边界 (与 fragment_generation 一致)
MW_RANGE = (250, 550)
LOGP_MAX = 5.5
# PAM 药效团 (简化自 M4 别构口袋 VU0467154 接触模式, 见 PIPELINE.md §3)
MIN_AROM, MIN_HBA, MIN_HBD = 2, 2, 1


# ---------------- tokenizer ----------------
def build_vocab(lines):
    chars = set()
    for s in lines:
        chars.update(s)
    return {c: i for i, c in enumerate(SPECIAL + sorted(chars))}


def encode(smi, c2i):
    return [BOS] + [c2i.get(c, UNK) for c in smi] + [EOS]


def decode(ids, i2c):
    return "".join(i2c[i] for i in ids if i not in (PAD, BOS, EOS, UNK))


# ---------------- dataset ----------------
class SmilesDataset(Dataset):
    def __init__(self, smis, c2i, max_len=MAX_LEN):
        self.seqs = [encode(s, c2i)[:max_len] for s in smis]

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, i):
        return self.seqs[i]


def collate_fn(batch):
    """动态 padding: 每个 batch pad 到本 batch 最大长度, 避免全 120 浪费."""
    T = max(len(s) for s in batch) - 1
    xs, ys = [], []
    for s in batch:
        xs.append(s[:-1] + [PAD] * (T - len(s) + 1))
        ys.append(s[1:] + [PAD] * (T - len(s) + 1))
    return torch.tensor(xs, dtype=torch.long), torch.tensor(ys, dtype=torch.long)


# ---------------- model ----------------
class SmilesLSTM(nn.Module):
    def __init__(self, vocab_size, embed=128, hidden=256, layers=2, dropout=0.2):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, embed, padding_idx=PAD)
        self.lstm = nn.LSTM(embed, hidden, layers, batch_first=True, dropout=dropout)
        self.fc = nn.Linear(hidden, vocab_size)

    def forward(self, x):
        return self.fc(self.lstm(self.embed(x))[0])


def train_model(model, loader, steps, lr, device, log_every=50):
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    crit = nn.CrossEntropyLoss(ignore_index=PAD)
    model.train()
    step = 0
    running = 0.0
    while step < steps:
        for x, y in loader:
            if step >= steps:
                break
            x, y = x.to(device), y.to(device)
            logits = model(x)
            loss = crit(logits.reshape(-1, logits.size(-1)), y.reshape(-1))
            opt.zero_grad()
            loss.backward()
            opt.step()
            running += loss.item()
            step += 1
            if step % log_every == 0:
                print(f"    step {step}/{steps}  loss {running/log_every:.3f}", flush=True)
                running = 0.0
    return step


@torch.no_grad()
def sample(model, c2i, i2c, n, temp=1.0, device="cpu", seed=None, batch=64):
    """批量自回归采样: 维护 LSTM 隐藏状态 (h/c) 逐 token 生成.
    关键: 必须把上一时刻的 hidden state 传入下一时刻, 否则模型每步失忆,
    只能生成碎片 (曾导致采样全是短分子)."""
    if seed is not None:
        torch.manual_seed(seed)
    model.eval()
    outs = []
    seen = set()
    B = batch
    cur = torch.full((B, 1), BOS, dtype=torch.long, device=device)
    ids = [[BOS] for _ in range(B)]
    alive = [True] * B
    h = torch.zeros(model.lstm.num_layers, B, model.lstm.hidden_size, device=device)
    c = torch.zeros_like(h)
    alive_mask = torch.ones(B, 1, dtype=torch.bool, device=device)

    while len(outs) < n:
        emb = model.embed(cur)
        out, (h, c) = model.lstm(emb, (h, c))
        logits = model.fc(out[:, -1, :]) / temp
        p = torch.softmax(logits, dim=-1)
        nxt = torch.multinomial(p, 1).squeeze(-1)
        next_ids = torch.full((B, 1), EOS, dtype=torch.long, device=device)
        any_alive = False
        for b in range(B):
            if not alive[b]:
                continue
            t = nxt[b].item()
            if t in (EOS, PAD) or len(ids[b]) >= MAX_LEN:
                alive[b] = False
                alive_mask[b, 0] = False
                smi = decode(ids[b], i2c)
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    cs = Chem.MolToSmiles(mol)
                    if cs not in seen:
                        seen.add(cs)
                        outs.append(cs)
            else:
                ids[b].append(t)
                next_ids[b, 0] = t
                any_alive = True
        if not any_alive:
            # 整批结束 -> 重置
            cur = torch.full((B, 1), BOS, dtype=torch.long, device=device)
            ids = [[BOS] for _ in range(B)]
            alive = [True] * B
            alive_mask = torch.ones(B, 1, dtype=torch.bool, device=device)
            h = torch.zeros_like(h)
            c = torch.zeros_like(c)
        else:
            # 死亡的序列 hidden 清零, 避免干扰后续
            mask3d = alive_mask.t().unsqueeze(-1).to(h.dtype)  # (1, B, 1)
            h = h * mask3d
            c = c * mask3d
            cur = next_ids
    return outs[:n]


# ---------------- filters ----------------
def druglike_ok(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    return MW_RANGE[0] <= mw <= MW_RANGE[1] and logp <= LOGP_MAX and hbd <= 5 and hba <= 10 and tpsa <= 140


def pharmacophore_ok(mol):
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    return rings >= MIN_AROM and hba >= MIN_HBA and hbd >= MIN_HBD


# ---------------- helpers ----------------
def load_smis(path):
    if not Path(path).exists():
        return []
    return [l.strip() for l in Path(path).read_text(encoding="utf-8").splitlines() if l.strip()]


def save_artifacts(path, model, c2i):
    torch.save(model.state_dict(), path)
    with open(path.with_suffix(".json"), "w", encoding="utf-8") as fh:
        json.dump(c2i, fh)


def load_artifacts(path, model_cls):
    c2i = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    model = model_cls(len(c2i))
    model.load_state_dict(torch.load(path, map_location="cpu"))
    return model, c2i


def _actives():
    smis = load_smis(ACTIVES)
    if not smis:
        smis = _DEFAULT_ACTIVES
        ACTIVES.parent.mkdir(parents=True, exist_ok=True)
        ACTIVES.write_text("\n".join(smis) + "\n", encoding="utf-8")
    return [s for s in smis if Chem.MolFromSmiles(s) is not None]


def cmd_pretrain():
    corpus = load_smis(CORPUS)
    if not corpus:
        raise SystemExit(f"无语料 {CORPUS}, 先跑 download_pubchem.py")
    actives = _actives()
    c2i = build_vocab(corpus + actives)
    i2c = {v: k for k, v in c2i.items()}
    print(f"语料 {len(corpus)} | 阳性 {len(actives)} | vocab {len(c2i)}")

    ds = SmilesDataset(corpus, c2i)
    loader = DataLoader(ds, batch_size=256, shuffle=True, num_workers=0, collate_fn=collate_fn)
    model = SmilesLSTM(len(c2i))
    print(f"模型参数 {sum(p.numel() for p in model.parameters()):,}")
    steps = 10 * len(loader)
    print(f"预训练 {steps} 步 (~10 epoch) ...")
    train_model(model, loader, steps, lr=1e-3, device="cpu")
    MODELS.mkdir(parents=True, exist_ok=True)
    save_artifacts(MODELS / "smiles_lstm.pt", model, c2i)
    print(f"已存 {MODELS/'smiles_lstm.pt'}")


def cmd_finetune():
    pt = MODELS / "smiles_lstm.pt"
    if not pt.exists():
        raise SystemExit("先跑 pretrain")
    model, c2i = load_artifacts(pt, SmilesLSTM)
    actives = _actives()
    ds = SmilesDataset(actives, c2i)
    loader = DataLoader(ds, batch_size=max(1, len(actives)), shuffle=True, num_workers=0, collate_fn=collate_fn)
    steps = 25
    print(f"fine-tune 到 {len(actives)} 个阳性, {steps} 步 (低lr防遗忘) ...")
    train_model(model, loader, steps, lr=1e-4, device="cpu", log_every=20)
    save_artifacts(MODELS / "smiles_lstm_ft.pt", model, c2i)
    print(f"已存 {MODELS/'smiles_lstm_ft.pt'}")


def cmd_generate():
    ft = MODELS / "smiles_lstm_ft.pt"
    pt = MODELS / "smiles_lstm.pt"
    src = ft if ft.exists() else pt
    if not src.exists():
        raise SystemExit("先跑 pretrain (最好也 finetune)")
    model, c2i = load_artifacts(src, SmilesLSTM)
    i2c = {v: k for k, v in c2i.items()}
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
    batch = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    print(f"采样 (来自 {src.name}, n={n}, batch={batch}) ...")
    raw = sample(model, c2i, i2c, n=n, temp=0.7, batch=batch)
    valid = [s for s in raw if Chem.MolFromSmiles(s) is not None]
    print(f"  采样 3000 -> 有效 {len(valid)} ({100*len(valid)/3000:.0f}%)")

    passed = []
    for s in valid:
        mol = Chem.MolFromSmiles(s)
        if druglike_ok(mol) and pharmacophore_ok(mol):
            passed.append({
                "smiles": s, "mw": round(Descriptors.MolWt(mol), 1),
                "logp": round(Crippen.MolLogP(mol), 2),
                "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "hbd": Lipinski.NumHDonors(mol),
            })
    passed.sort(key=lambda r: (abs(r["logp"] - 3.0), abs(r["mw"] - 350)))

    import csv
    out = MODELS / "lstm_generated_pam_analogs.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "mw", "logp", "aromatic_rings", "hba", "hbd"])
        w.writeheader()
        w.writerows(passed)
    print(f"  过滤后 {len(passed)} 个候选 -> {out}")

    # 阳性种子是否在候选里 (非必须, 仅参考)
    actives = set(_actives())
    hit = sum(1 for r in passed if r["smiles"] in actives)
    print(f"  阳性种子回填 {hit}/{len(actives)}")

    print("\nTop 12:")
    for i, r in enumerate(passed[:12], 1):
        print(f"  {i:>2} {r['smiles'][:70]}  MW {r['mw']} logP {r['logp']}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "generate"
    {"pretrain": cmd_pretrain, "finetune": cmd_finetune, "generate": cmd_generate}[cmd]()
