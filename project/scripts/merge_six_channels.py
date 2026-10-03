"""六通道终合并：Vina×{7TRQ,7TRP,7TRS} + GNINA-vinardo×{同} → rank-percentile 共识。

诚实标注：开源引擎六通道事后补充分析（非冻结 PACER-XR 六路）。
通道语义：每通道取该受体上的 best docking affinity → 通道内 rank percentile
（能量越低越好 → percentile 越高越好）→ 六通道等权平均 = XR 同款数学。
"""
import glob, os, re
import numpy as np, pandas as pd

R = r'D:\CLC\project\results\pacer_rerun_v01'
PAT = re.compile(r'REMARK VINA RESULT:\s+(-?\d+\.?\d*)')
RECEPTORS = ['7TRQ', '7TRP', '7TRS']

base = pd.read_csv(R + r'\module2_top200_output.csv')

# --- Vina 三通道（CSV 已恢复/将恢复）---
vina_scores = {}
for name in RECEPTORS:
    csv = f'{R}\\module3_{name}_results.csv'
    if os.path.exists(csv):
        d = pd.read_csv(csv)
        col = [c for c in d.columns if c.startswith('vina_')]
        if col:
            vina_scores[name] = d[col[0]].to_numpy()
            # top200 CSV 行序 == base 行序（同一读取源）
print('Vina 通道:', {k: int(np.isnan(v).sum()) for k, v in vina_scores.items() if isinstance(v, np.ndarray)} if vina_scores else '无',
      '| 有效:', {k: int(np.sum(~np.isnan(v))) for k, v in vina_scores.items()})

# --- GNINA 三通道（解析 pose 的 minimizedAffinity = vinardo 能量）---
PAT_GNINA = re.compile(r'REMARK minimizedAffinity\s+(-?\d+\.?\d*)')
gnina_scores, gnina_cnn = {}, {}
for name in RECEPTORS:
    scores = np.full(len(base), np.nan)
    cnn = np.full(len(base), np.nan)
    for f in glob.glob(f'{R}\\gnina_tmp\\{name}_out*.pdbqt'):
        idx = int(re.search(r'out(\d+)', os.path.basename(f)).group(1))
        txt = open(f, errors='replace').read()
        m = PAT_GNINA.search(txt)
        mc = re.search(r'REMARK CNNaffinity\s+(-?\d+\.?\d*)', txt)
        if m and idx < len(base):
            scores[idx] = float(m.group(1))
        if mc and idx < len(base):
            cnn[idx] = float(mc.group(1))
    gnina_scores[name] = scores
    gnina_cnn[name] = cnn
    print(f'GNINA-{name}: vinardo {np.sum(~np.isnan(scores))}/200  CNNaff {np.sum(~np.isnan(cnn))}/200')

# --- 组装与共识 ---
out = base.copy()
channels = {}
for name in RECEPTORS:
    if name in vina_scores:
        out[f'vina_{name}'] = vina_scores[name]
        channels[f'vina_{name}'] = vina_scores[name]
    if np.sum(~np.isnan(gnina_scores[name])) > 0:
        out[f'gnina_{name}'] = gnina_scores[name]
        channels[f'gnina_{name}'] = gnina_scores[name]

# --- GNINA rescore 通道（score_only 的 CNNaffinity / vinardo）---
PAT_CNN = re.compile(r'CNNaffinity:?\s*(-?\d+\.?\d*)')
PAT_MIN = re.compile(r'minimizedAffinity:?\s*(-?\d+\.?\d*)')
for name in RECEPTORS:
    cnn = np.full(len(base), np.nan)
    via = np.full(len(base), np.nan)
    for f in glob.glob(f'{R}\\gnina_tmp\\{name}_rescore*.txt'):
        idx = int(re.search(r'rescore(\d+)', os.path.basename(f)).group(1))
        txt = open(f, errors='replace').read()
        mc, mv = PAT_CNN.search(txt), PAT_MIN.search(txt)
        if mc and idx < len(base):
            cnn[idx] = float(mc.group(1))
        if mv and idx < len(base):
            via[idx] = float(mv.group(1))
    if np.sum(~np.isnan(cnn)) > 0:
        channels[f'gninaCNN_{name}'] = cnn
        out[f'gninaCNN_{name}'] = cnn
    if np.sum(~np.isnan(via)) > 0:
        channels[f'gninaVia_{name}'] = via
        out[f'gninaVia_{name}'] = via
    print(f'GNINA-rescore-{name}: CNN {np.sum(~np.isnan(cnn))}  vinardo {np.sum(~np.isnan(via))}')

pct = {}
for ch, arr in channels.items():
    s = pd.Series(arr)
    pct[ch] = s.rank(pct=True).to_numpy()          # 能量低 -> rank 高 -> pct 高
stack = np.stack([p for p in pct.values()])
cons = np.nanmean(stack, axis=0) if len(pct) else np.full(len(base), np.nan)
out['six_channel_consensus'] = cons
out['consensus_rank'] = pd.Series(cons).rank(ascending=False, method='first').astype(int)
out = out.sort_values('consensus_rank').reset_index(drop=True)
out.to_csv(R + r'\module3_six_channel_merged.csv', index=False)

# --- MD 三候选定位 ---
from rdkit import Chem
from rdkit.RDLogger import DisableLog
DisableLog('rdApp.*')
canon = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s))) or str(s)
md3 = {'PACER0010': canon('Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)Nc3cnn(C4Cc5ccccc5C4)c3)c(N)c12'),
       'PACER0073': canon('Cc1c(Cl)c2nnc(C)n2c2sc(N3CC(NC(=O)c4cc(F)nc(F)c4)C3)c(N)c12'),
       'PACER0014': canon('Cc1c(Cl)c2nnc(C)n2c2sc(N3CC(NC(=O)c4ccncc4Cl)C3)c(N)c12'),
       'PACER0027': canon('Cc1c(Cl)c2nncn2c2sc(N3CC(N4CC(NC5Cc6ccccc6C5)C4)C3)c(N)c12')}
cs = {v: k for k, v in md3.items()}
out['canon'] = [canon(s) for s in out.smiles]
print()
print('===== MD 候选六通道共识位置 =====')
for _, r in out.iterrows():
    if r.canon in cs:
        print('%s -> 共识排名 #%d/%d（共识分 %.3f）' % (cs[r.canon], r.consensus_rank, len(out), r.six_channel_consensus))
print()
print('===== 共识 top-15 =====')
for _, r in out.head(15).iterrows():
    tag = ' <== ' + cs[r.canon] if r.canon in cs else ''
    print('#%-3d 共识%.3f  模块2#%-4d %s%s' % (r.consensus_rank, r.six_channel_consensus, r.final_rank, str(r.smiles)[:40], tag))
out.head(15).to_csv(R + r'\final_dozen_candidates.csv', index=False)
