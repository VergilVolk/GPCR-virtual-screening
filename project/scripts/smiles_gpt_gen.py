"""
SMILES GPT 生成器 (字符级 decoder-only Transformer, 纯 CPU)
=============================================================
与 smiles_lstm_gen.py / smiles_diffusion_gen.py 对齐的迁移学习范式,
用于生成式架构对比: LSTM (RNN) vs GPT (Transformer) vs Diffusion (masked).

  预训练  PubChem 类药语料上训 next-token prediction (学化学语法)
  fine-tune  用少量 PAM 阳性微调 (把分布掰向 PAM 化学型)
  generate  自回归采样 -> RDKit 校验 -> 类药+药效团过滤

模型: 2 层 causal Transformer (d_model=128, 4 heads), 字符级 token.

用法:
  python scripts/smiles_gpt_gen.py pretrain
  python scripts/smiles_gpt_gen.py finetune
  python scripts/smiles_gpt_gen.py generate [N] [batch]
输出: results/generated/gpt_generated_pam_analogs.csv
"""
from __future__ import annotations
import sys
import json
import csv
from pathlib import Path

import torch
import torch.nn as nn

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
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

SPECIAL = ["<PAD>", "<BOS>", "<EOS>", "<UNK>"]
PAD, BOS, EOS, UNK = 0, 1, 2, 3
MAX_LEN = 120

# 类药边界 (与 LSTM/扩散一致)
MW_RANGE = (250, 550)
LOGP_MAX = 5.5
MIN_AROM, MIN_HBA, MIN_HBD = 2, 2, 1


def build_vocab(lines):
    chars = set()
    for s in lines:
        chars.update(s)
    return {c: i for i, c in enumerate(SPECIAL + sorted(chars))}


def encode(smi, c2i):
    return [BOS] + [c2i.get(c, UNK) for c in smi] + [EOS]


def decode(ids, i2c):
    return "".join(i2c.get(int(i), "") for i in ids if int(i) not in (PAD, BOS, EOS, UNK))


def load_smis(path):
    if not Path(path).exists():
        return []
    return [l.strip() for l in Path(path).read_text(encoding="utf-8").splitlines() if l.strip()]


def _actives():
    smis = load_smis(ACTIVES)
    if not smis:
        raise SystemExit(f"无阳性 {ACTIVES}")
    return [s for s in smis if Chem.MolFromSmiles(s) is not None]


# ---------------- 模型 ----------------
class SmilesGPT(nn.Module):
    """字符级 decoder-only Transformer (causal mask)."""

    def __init__(self, vocab_size, d_model=128, nhead=4, layers=2, max_len=MAX_LEN):
        super().__init__()
        self.tok = nn.Embedding(vocab_size, d_model, padding_idx=PAD)
        self.pos = nn.Parameter(torch.zeros(1, max_len, d_model))
        nn.init.normal_(self.pos, std=0.02)
        enc_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=512,
            dropout=0.1, batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=layers)
        self.out = nn.Linear(d_model, vocab_size)

    def forward(self, x):
        B, L = x.shape
        mask = torch.triu(torch.full((L, L), float("-inf"), device=x.device), diagonal=1)
        h = self.tok(x) + self.pos[:, :L, :]
        h = self.encoder(h, mask=mask)
        return self.out(h)


# ---------------- 数据/训练 ----------------
class SmilesDS:
    def __init__(self, smis, c2i):
        self.seqs = [encode(s, c2i)[:MAX_LEN] for s in smis]

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, i):
        s = self.seqs[i]
        x = s[:-1] + [PAD] * (MAX_LEN - len(s) + 1)
        y = s[1:] + [PAD] * (MAX_LEN - len(s) + 1)
        return torch.tensor(x, dtype=torch.long), torch.tensor(y, dtype=torch.long)


def train_epoch(model, ds, opt, batch=256, steps=None, log_every=200):
    model.train()
    crit = nn.CrossEntropyLoss(ignore_index=PAD)
    n = len(ds)
    order = torch.randperm(n)
    done = 0
    total_loss = 0.0
    n_batch = 0
    i = 0
    while i < n:
        idx = order[i:i + batch]
        xs = torch.stack([ds[j][0] for j in idx])
        ys = torch.stack([ds[j][1] for j in idx])
        logits = model(xs)
        loss = crit(logits.reshape(-1, logits.size(-1)), ys.reshape(-1))
        opt.zero_grad()
        loss.backward()
        opt.step()
        total_loss += loss.item()
        n_batch += 1
        done += len(idx)
        if n_batch % log_every == 0:
            print(f"    {done}/{n}  loss {total_loss/n_batch:.3f}", flush=True)
        if steps is not None and n_batch >= steps:
            break
        i += batch
    return total_loss / max(1, n_batch)


@torch.no_grad()
def sample(model, c2i, i2c, n, batch=128, temp=0.9, seed=None):
    """自回归采样: 每步 forward 全部序列 (causal), 取最后 token."""
    if seed is not None:
        torch.manual_seed(seed)
    model.eval()
    outs, seen = [], set()
    B = batch
    while len(outs) < n:
        ids = [[BOS] for _ in range(B)]
        alive = [True] * B
        for step in range(MAX_LEN):
            max_len = max(len(s) for s in ids)
            x = torch.tensor([s + [PAD] * (max_len - len(s)) for s in ids], dtype=torch.long)
            logits = model(x)[:, -1, :] / temp
            p = torch.softmax(logits, dim=-1)
            nxt = torch.multinomial(p, 1).squeeze(-1)
            any_alive = False
            for b in range(B):
                if not alive[b]:
                    continue
                t = nxt[b].item()
                if t in (EOS, PAD) or len(ids[b]) >= MAX_LEN:
                    alive[b] = False
                    smi = decode(ids[b], i2c)
                    mol = Chem.MolFromSmiles(smi)
                    if mol is not None:
                        cs = Chem.MolToSmiles(mol)
                        if cs not in seen:
                            seen.add(cs)
                            outs.append(cs)
                else:
                    ids[b].append(t)  # 存 token id (int), decode 时再转字符
                    any_alive = True
            if not any_alive:
                break
    return outs[:n]


# ---------------- filters (与 LSTM/扩散一致) ----------------
def druglike_ok(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    return (MW_RANGE[0] <= mw <= MW_RANGE[1] and logp <= LOGP_MAX
            and hbd <= 5 and hba <= 10 and tpsa <= 140)


def pharmacophore_ok(mol):
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    return rings >= MIN_AROM and hba >= MIN_HBA and hbd >= MIN_HBD


# ---------------- 命令 ----------------
def cmd_pretrain():
    corpus = load_smis(CORPUS)
    actives = _actives()
    c2i = build_vocab(corpus + actives)
    i2c = {v: k for k, v in c2i.items()}
    print(f"语料 {len(corpus)} | 阳性 {len(actives)} | vocab {len(c2i)}")
    ds = SmilesDS(corpus, c2i)
    model = SmilesGPT(len(c2i))
    print(f"模型参数 {sum(p.numel() for p in model.parameters()):,}")
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    n_epochs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    for epoch in range(n_epochs):
        loss = train_epoch(model, ds, opt)
        print(f"  epoch {epoch+1}/{n_epochs}  loss {loss:.3f}", flush=True)
    torch.save(model.state_dict(), MODELS / "smiles_gpt.pt")
    with open(MODELS / "smiles_gpt.json", "w", encoding="utf-8") as fh:
        json.dump(c2i, fh)
    print(f"已存 {MODELS/'smiles_gpt.pt'}")


def cmd_finetune():
    pt = MODELS / "smiles_gpt.pt"
    if not pt.exists():
        raise SystemExit("先跑 pretrain")
    c2i = json.loads((MODELS / "smiles_gpt.json").read_text(encoding="utf-8"))
    model = SmilesGPT(len(c2i))
    model.load_state_dict(torch.load(pt, map_location="cpu"))
    actives = _actives()
    ds = SmilesDS(actives, c2i)
    opt = torch.optim.AdamW(model.parameters(), lr=5e-5, weight_decay=1e-4)
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    loss = train_epoch(model, ds, opt, batch=len(ds), steps=steps, log_every=20)
    print(f"  ft loss {loss:.3f}")
    torch.save(model.state_dict(), MODELS / "smiles_gpt_ft.pt")
    print(f"已存 {MODELS/'smiles_gpt_ft.pt'}")


def cmd_generate():
    ft = MODELS / "smiles_gpt_ft.pt"
    pt = MODELS / "smiles_gpt.pt"
    src = ft if ft.exists() else pt
    c2i = json.loads((MODELS / "smiles_gpt.json").read_text(encoding="utf-8"))
    model = SmilesGPT(len(c2i))
    model.load_state_dict(torch.load(src, map_location="cpu"))
    i2c = {v: k for k, v in c2i.items()}
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    batch = int(sys.argv[3]) if len(sys.argv) > 3 else 128
    print(f"采样 (来自 {src.name}, n={n}) ...")
    raw = sample(model, c2i, i2c, n=n, batch=batch, temp=0.9)
    valid = [s for s in raw if Chem.MolFromSmiles(s) is not None]
    print(f"  采样 {n} -> 有效 {len(valid)} ({100*len(valid)/max(1,n):.0f}%)")

    passed = []
    for s in valid:
        mol = Chem.MolFromSmiles(s)
        if druglike_ok(mol) and pharmacophore_ok(mol):
            passed.append({"smiles": s, "mw": round(Descriptors.MolWt(mol), 1),
                           "logp": round(Crippen.MolLogP(mol), 2),
                           "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
                           "hba": Lipinski.NumHAcceptors(mol),
                           "hbd": Lipinski.NumHDonors(mol)})
    passed.sort(key=lambda r: (abs(r["logp"] - 3.0), abs(r["mw"] - 350)))
    out = MODELS / "gpt_generated_pam_analogs.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "mw", "logp", "aromatic_rings", "hba", "hbd"])
        w.writeheader()
        w.writerows(passed)
    print(f"  过滤后 {len(passed)} 个候选 -> {out}")
    print("\nTop 12:")
    for i, r in enumerate(passed[:12], 1):
        print(f"  {i:>2} {r['smiles'][:70]}  MW {r['mw']} logP {r['logp']}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "generate"
    {"pretrain": cmd_pretrain, "finetune": cmd_finetune, "generate": cmd_generate}[cmd]()
