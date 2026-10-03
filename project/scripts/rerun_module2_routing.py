"""模块2 真算：2,605 生成分子 × 双权重路由 → top-100 输出库。

2023 侧（M4-safe 主路由）：官方 ckpt 表征 × GPCR-LOTO(M4留出) 3种子投影
  vs 冻结 M4R_cluster0 口袋表征（gpcr_drugclip_screening_v01 权威 npz）。
2026 侧（第二意见，若有口袋表征）：Science2026 表征 × famaug 3种子 vs 晶体三口袋。
"""
import numpy as np, pandas as pd, torch, sys
from pathlib import Path
sys.path.insert(0, r'D:\CLC\project\scripts')
from finetune_drugclip_muscarinic_triplet import Proj

R = Path(r'D:\CLC\project\results\pacer_rerun_v01')
SEEDS = (20260925, 20260926, 20260927)

z23 = np.load(R / 'gen2605_2023.npz', allow_pickle=False)
ids23 = [str(v) for v in z23['molecule_ids']]
mol23 = torch.as_tensor(z23['molecule_representations'].astype(np.float32))
lib = np.load(r'D:\CLC\project\results\gpcr_drugclip_screening_v01\ensemble_embeddings.npz', allow_pickle=False)
gp = [str(v) for v in lib['pocket_ids']]
pk23 = torch.as_tensor(lib['pocket_representations'][gp.index('M4R_cluster0')].astype(np.float32))

per_seed = []
with torch.inference_mode():
    for s in SEEDS:
        d = torch.load(rf'D:\CLC\project\results\gpcr_drugclip_retrieval_v01\finetune_result.seed{s}.projection.pt',
                       map_location='cpu', weights_only=False)
        pm, pp = Proj(d['mol_project']).eval(), Proj(d['pocket_project']).eval()
        sc = (pm(mol23) @ pp(pk23)).numpy()
        per_seed.append(sc)
s23 = np.mean(per_seed, axis=0)

# ---- 2026 第二意见（找到 2026 空间晶体口袋表征才跑）----
s26 = None
cand_files = [
    r'D:\CLC\project\results\drugclip_science2026\m4_crystal_pockets_v01\pockets.npz',
    r'D:\CLC\project\results\pacer_prospective_run_20261002_v01\crystal_pockets_2026.npz',
]
import os
pk26_path = next((f for f in cand_files if os.path.exists(f)), None)
if pk26_path is None:
    # 归档里找
    import subprocess
    hit = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', '-r', '--name-only',
                          'origin/incident/pacer-xr-protocol-deviation-20261003'],
                         capture_output=True, text=True).stdout.splitlines()
    hits = [h for h in hit if '2026' in h and 'pocket' in h.lower() and h.endswith('.npz')]
    print('归档 2026 口袋 npz 候选:', hits[:4] if hits else '无')
else:
    print('本地 2026 口袋:', pk26_path)

z26 = np.load(R / 'gen2605_2026.npz', allow_pickle=False)
ids26 = [str(v) for v in z26['molecule_ids']]
assert ids26 == ids23, '两侧分子顺序不一致!'
mol26 = torch.as_tensor(z26['molecule_representations'].astype(np.float32))

# ---- M4-safe 路由：仅 2023 分数定排名 ----
df = pd.DataFrame({'molecule_id': ids23, 'score_2023_gpcr_loto': s23})
df['rankpct_2023'] = df.score_2023_gpcr_loto.rank(pct=True)
df['pacer_binding_score'] = df.rankpct_2023
df['final_rank'] = df.rankpct_2023.rank(ascending=False, method='first').astype(int)
df = df.sort_values('final_rank').reset_index(drop=True)
gen = pd.read_csv(R / 'gen2605.csv')
from rdkit import Chem
from rdkit.RDLogger import DisableLog
DisableLog('rdApp.*')
canon = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s))) or str(s)
gmap = {}
for mid, smi in zip(gen.molecule_id, gen.smiles):
    gmap[canon(smi)] = (mid, smi)
df['canon'] = [canon(v) for v in df.molecule_id]          # npz 的 ID 就是 SMILES
df['smiles'] = [gmap.get(c, (None, None))[1] for c in df.canon]
df['gen_id'] = [gmap.get(c, (None, None))[0] for c in df.canon]
assert df.smiles.notna().all(), '仍有未映射分子!'
df = df.drop(columns=['canon'])
top100 = df.head(100).copy()
top100.to_csv(R / 'module2_top100_library.csv', index=False)
df.to_csv(R / 'module2_routed2605_full.csv', index=False)

# 与冻结 200 分子对账：top100 里有多少是当时的 200？
import subprocess, io
h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', 'origin/incident/pacer-xr-protocol-deviation-20261003', '--',
                    'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01/predock_portfolio.csv'],
                   capture_output=True, text=True).stdout.split()
raw = subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True).stdout
port = pd.read_csv(io.BytesIO(raw))
port_smiles = set(port.canonical_smiles.astype(str)) if 'canonical_smiles' in port else set(port.smiles.astype(str))
arch = port  # 200 组合表同时用于三候选对账
canon_ids = [canon(v) for v in df.smiles]
in_top100 = sum(canon(v) in {canon(x) for x in port_smiles} for v in df.head(100).smiles)
top10_in_port = sum(canon(v) in {canon(x) for x in port_smiles} for v in df.head(10).smiles)
print(f'2,605 → top-100 完成；top-100 中属于当年 200 入栈的: {in_top100}')
print(f'top-10 中属于当年 200 的: {top10_in_port}')
# 三候选的直接对账（用 corrected_v1 权威排名表，有 pair_id）
h2 = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', 'origin/incident/pacer-xr-protocol-deviation-20261003', '--',
                     'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv'],
                    capture_output=True, text=True).stdout.split()
raw2 = subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h2[2]], capture_output=True).stdout
ranking = pd.read_csv(io.BytesIO(raw2))
arch_canon = {canon(s): pid for pid, s in zip(ranking.pair_id, ranking.canonical_smiles)
              if pid in ('PACER0010', 'PACER0073', 'PACER0027', 'PACER0014')}
for cs, pid in arch_canon.items():
    row = df[[canon(v) == cs for v in df.smiles]]
    if len(row):
        r = row.iloc[0]
        print('%s -> 新分数 %.4f  新排名 %d / 2605' % (pid, r.score_2023_gpcr_loto, r.final_rank))
    else:
        print(pid, '未找到!')
print('新 top-10:')
for _, r in df.head(10).iterrows():
    print('  %s %.4f %s' % (r.gen_id, r.score_2023_gpcr_loto,
                            'IN-200' if canon(r.smiles) in {canon(x) for x in port_smiles} else 'new'))
