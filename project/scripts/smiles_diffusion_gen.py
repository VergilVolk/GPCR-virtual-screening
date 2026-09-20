"""
SMILES 离散扩散生成器 (absorbing/masked diffusion, D3PM 简化版, 纯 CPU)
========================================================================
与 smiles_lstm_gen.py 完全对齐的迁移学习范式, 用于生成式方法学对比:

  预训练  PubChem 类药语料上训练 masked diffusion (学药物分子语法)
  fine-tune  用少量 PAM 阳性微调 (把分布掰向 PAM 化学型)
  generate  迭代去噪采样 -> RDKit 校验 -> 类药+药效团过滤

模型: 小 Transformer encoder (d_model=128, 2 层, 4 heads) + 时间 embedding
前向: 按扩散时间表随机 mask token; 反向: 预测被 mask 位置的原始 token
采样: 从全 [MASK] 开始, T 步迭代揭晓 (每步揭晓剩余位置的一部分)

用法:
  python scripts/smiles_diffusion_gen.py pretrain
  python scripts/smiles_diffusion_gen.py finetune
  python scripts/smiles_diffusion_gen.py generate [N] [batch]
输出: results/generated/diffusion_generated_pam_analogs.csv + 对比 JSON
"""
from __future__ import annotations
import sys
import json
import math
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

# 特殊 token: PAD / MASK / UNK
PAD, MASK, UNK = 0, 1, 2
MAX_LEN = 80  # 固定长度 (比 LSTM 120 紧凑, 覆盖中位长度 31 的 ~2.5 倍)

# 类药边界 (与 LSTM/片段生成一致)
MW_RANGE = (250, 550)
LOGP_MAX = 5.5
MIN_AROM, MIN_HBA, MIN_HBD = 2, 2, 1

N_TIMESTEPS = 50  # 训练用离散时间步 (t=1..N_TIMESTEPS, 越大 mask 越多)


def build_vocab(lines):
    chars = set()
    for s in lines:
        chars.update(s)
    return {c: i + 3 for i, c in enumerate(sorted(chars))}  # 0=PAD 1=MASK 2=UNK


def encode(smi, c2i):
    ids = [c2i.get(c, UNK) for c in smi[:MAX_LEN]]
    ids = ids + [PAD] * (MAX_LEN - len(ids))
    return torch.tensor(ids, dtype=torch.long)


def decode(ids, i2c):
    return "".join(i2c.get(int(i), "") for i in ids if int(i) not in (PAD, MASK, UNK))


# ---------------- 时间 embedding ----------------
class TimeEmbed(nn.Module):
    def __init__(self, d_model, n_timesteps=N_TIMESTEPS):
        super().__init__()
        self.emb = nn.Embedding(n_timesteps + 1, d_model)

    def forward(self, t):
        return self.emb(t)  # (B, d)


# ---------------- 模型 ----------------
class SmilesDiffusion(nn.Module):
    def __init__(self, vocab_size, d_model=128, nhead=4, layers=2, max_len=MAX_LEN):
        super().__init__()
        self.token_emb = nn.Embedding(vocab_size, d_model, padding_idx=PAD)
        self.time_emb = TimeEmbed(d_model)
        self.pos = nn.Parameter(torch.zeros(1, max_len, d_model))
        nn.init.normal_(self.pos, std=0.02)
        enc_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=512,
            dropout=0.1, batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=layers)
        self.out = nn.Linear(d_model, vocab_size)

    def forward(self, x, t):
        # x: (B, L) token ids; t: (B,) 时间步
        B, L = x.shape
        h = self.token_emb(x) + self.pos[:, :L, :] + self.time_emb(t).unsqueeze(1)
        h = self.encoder(h)
        return self.out(h)  # (B, L, V)


# ---------------- 数据 ----------------
def load_smis(path):
    if not Path(path).exists():
        return []
    return [l.strip() for l in Path(path).read_text(encoding="utf-8").splitlines() if l.strip()]


def _actives():
    smis = load_smis(ACTIVES)
    if not smis:
        raise SystemExit(f"无阳性 {ACTIVES}, 请先确认 actives.smi")
    return [s for s in smis if Chem.MolFromSmiles(s) is not None]


def mask_and_target(ids, mask_ratio, rng):
    """按 mask_ratio 随机 mask 位置, 返回 (输入, 目标, mask标记)."""
    B, L = ids.shape
    mask = (torch.rand(B, L, generator=rng) < mask_ratio.unsqueeze(1))
    # 保证每个样本至少 mask 1 个 (避免无监督退化)
    for b in range(B):
        if not mask[b].any():
            mask[b, torch.randint(0, L, (1,), generator=rng)] = True
    x = ids.clone()
    x[mask] = MASK
    return x, ids, mask


def train_step(model, opt, batch_ids, rng):
    """一个 batch: 随机 t -> mask_ratio, 预测被 mask 位置的 token."""
    B, L = batch_ids.shape
    t = torch.randint(1, N_TIMESTEPS + 1, (B,), generator=rng)
    mask_ratio = t / N_TIMESTEPS  # 均匀调度
    x, target, mask = mask_and_target(batch_ids, mask_ratio, rng)
    logits = model(x, t)
    loss_f = nn.CrossEntropyLoss(reduction="none")
    loss = loss_f(logits.reshape(-1, logits.size(-1)), target.reshape(-1))
    loss = (loss.reshape(B, L) * mask).sum() / mask.sum().clamp(min=1)
    opt.zero_grad()
    loss.backward()
    opt.step()
    return loss.item()


@torch.no_grad()
def sample(model, c2i, i2c, n, batch=128, T=10, temp=1.0, seed=None, device="cpu", max_rounds=100):
    """从全 [MASK] 迭代去噪: 每步揭晓约 1/T 的位置 (向量化, 按置信度)."""
    if seed is not None:
        torch.manual_seed(seed)
    model.eval()
    outs, seen = [], set()
    rounds = 0
    attempted = 0
    while len(outs) < n and rounds < max_rounds:
        rounds += 1
        attempted += batch
        cur = torch.full((batch, MAX_LEN), MASK, dtype=torch.long, device=device)
        for step in range(T):
            t = torch.full((batch,), N_TIMESTEPS - step, dtype=torch.long, device=device)
            logits = model(cur, t) / temp
            probs = torch.softmax(logits, dim=-1)
            # 每步揭晓: 所有 mask 位置中, 每分子取置信度 top ceil(rem/T')
            mask_pos = (cur == MASK)
            n_mask = mask_pos.sum(dim=1).clamp(min=1)  # (B,)
            k = torch.ceil(n_mask / (T - step)).clamp(min=1)  # (B,)
            for b in range(batch):
                if n_mask[b] == 0:
                    continue
                idx = mask_pos[b].nonzero(as_tuple=True)[0]
                p_sel = probs[b, idx].max(dim=-1).values
                topk = idx[p_sel.argsort(descending=True)[:int(k[b].item())]]
                cur[b, topk] = probs[b, topk].argmax(dim=-1)
        # 解码校验
        for b in range(batch):
            smi = decode(cur[b].tolist(), i2c)
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                cs = Chem.MolToSmiles(mol)
                if cs not in seen:
                    seen.add(cs)
                    outs.append(cs)
        if rounds % 5 == 0:
            print(f"  round {rounds}: 已累积 {len(outs)}/{n}", flush=True)
    return outs[:n], {"attempted": attempted, "valid_unique": len(outs),
                      "max_rounds": max_rounds, "terminated_by_cap": len(outs) < n}


# ---------------- filters (与 LSTM 一致) ----------------
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
    print(f"语料 {len(corpus)} | 阳性 {len(actives)} | vocab {len(c2i)+3}")

    ids = torch.stack([encode(s, c2i) for s in corpus])
    model = SmilesDiffusion(len(c2i) + 3)
    print(f"模型参数 {sum(p.numel() for p in model.parameters()):,}")
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    rng = torch.Generator().manual_seed(42)
    B = 256
    n_epochs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    steps = 0
    for epoch in range(n_epochs):
        perm = torch.randperm(len(ids), generator=rng)
        tot = 0.0
        for i in range(0, len(ids), B):
            batch = ids[perm[i:i + B]]
            tot += train_step(model, opt, batch, rng)
            steps += 1
        print(f"  epoch {epoch+1}/{n_epochs}  loss {tot/(len(ids)/B):.3f}", flush=True)
    torch.save(model.state_dict(), MODELS / "smiles_diff.pt")
    with open(MODELS / "smiles_diff.json", "w", encoding="utf-8") as fh:
        json.dump(c2i, fh)
    print(f"已存 {MODELS/'smiles_diff.pt'}")


def cmd_finetune():
    pt = MODELS / "smiles_diff.pt"
    if not pt.exists():
        raise SystemExit("先跑 pretrain")
    c2i = json.loads((MODELS / "smiles_diff.json").read_text(encoding="utf-8"))
    model = SmilesDiffusion(len(c2i) + 3)
    model.load_state_dict(torch.load(pt, map_location="cpu"))
    actives = _actives()
    ids = torch.stack([encode(s, c2i) for s in actives])
    opt = torch.optim.AdamW(model.parameters(), lr=5e-5, weight_decay=1e-4)
    rng = torch.Generator().manual_seed(7)
    B = max(1, len(ids))
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    for step in range(steps):
        batch = ids[torch.randperm(len(ids), generator=rng)[:B]]
        loss = train_step(model, opt, batch, rng)
        if (step + 1) % 100 == 0:
            print(f"  ft step {step+1}/{steps}  loss {loss:.3f}", flush=True)
    torch.save(model.state_dict(), MODELS / "smiles_diff_ft.pt")
    print(f"已存 {MODELS/'smiles_diff_ft.pt'}")


def cmd_generate():
    ft = MODELS / "smiles_diff_ft.pt"
    pt = MODELS / "smiles_diff.pt"
    src = ft if ft.exists() else pt
    c2i = json.loads((MODELS / "smiles_diff.json").read_text(encoding="utf-8"))
    model = SmilesDiffusion(len(c2i) + 3)
    model.load_state_dict(torch.load(src, map_location="cpu"))
    i2c = {v: k for k, v in c2i.items()}
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    batch = int(sys.argv[3]) if len(sys.argv) > 3 else 128
    print(f"采样 (来自 {src.name}, n={n}) ...")
    raw, sampling_audit = sample(model, c2i, i2c, n=n, batch=batch, temp=0.9, seed=42)
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
    out = MODELS / "diffusion_generated_pam_analogs.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "mw", "logp", "aromatic_rings", "hba", "hbd"])
        w.writeheader()
        w.writerows(passed)
    (MODELS / "diffusion_sampling_audit.json").write_text(
        json.dumps({**sampling_audit, "requested": n, "passed_filters": len(passed)}, indent=2),
        encoding="utf-8")
    print(f"  过滤后 {len(passed)} 个候选 -> {out}")
    print("\nTop 12:")
    for i, r in enumerate(passed[:12], 1):
        print(f"  {i:>2} {r['smiles'][:70]}  MW {r['mw']} logP {r['logp']}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "generate"
    {"pretrain": cmd_pretrain, "finetune": cmd_finetune, "generate": cmd_generate}[cmd]()
