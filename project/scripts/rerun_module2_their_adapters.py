import subprocess, io, json
import numpy as np, pandas as pd, torch, sys
from pathlib import Path
sys.path.insert(0, r'D:\CLC\project\scripts')
from finetune_drugclip_muscarinic_triplet import Proj
from rdkit import Chem
from rdkit.RDLogger import DisableLog
DisableLog('rdApp.*')

BR = 'origin/incident/pacer-xr-protocol-deviation-20261003'
R = Path(r'D:\CLC\project\results\pacer_rerun_v01')

def blob(path):
    h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', BR, '--', path],
                       capture_output=True, text=True).stdout.split()
    return subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True).stdout

A = 'project/artifacts/drugclip2023_m4_loto_pacer200_v01'
adapters = []
for s in (20260925, 20260926, 20260927):
    raw = blob(f'{A}/drugclip2023_m4_loto_seed{s}.projection.pt')
    adapters.append(torch.load(io.BytesIO(raw), map_location='cpu', weights_only=False))
print('adapters loaded:', len(adapters))

pk_raw = blob(f'{A}/pacer200_m4_2023_frozen_representations.npz')
z = np.load(io.BytesIO(pk_raw), allow_pickle=False)
print('pocket npz keys:', z.files)
pk = torch.as_tensor(z['m4_pocket_representation'].astype(np.float32))
print('pocket rep shape:', tuple(pk.shape))

mol = torch.as_tensor(np.load(R / 'gen2605_2023.npz', allow_pickle=False)['molecule_representations'].astype(np.float32))
per = []
with torch.inference_mode():
    for d in adapters:
        pm, pp = Proj(d['mol_project']).eval(), Proj(d['pocket_project']).eval()
        per.append((pm(mol) @ pp(pk)).numpy())
s23 = np.mean(per, axis=0)

canon = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s))) or str(s)
gen = pd.read_csv(R / 'gen2605.csv')
gmap = {canon(s): (mid, s) for mid, s in zip(gen.molecule_id, gen.smiles)}
ids = [str(v) for v in np.load(R / 'gen2605_2023.npz', allow_pickle=False)['molecule_ids']]
rows = {'gen_id': [], 'smiles': [], 'score': []}
for i in ids:
    c = canon(i)
    mid, smi = gmap.get(c, (None, None))
    rows['gen_id'].append(mid); rows['smiles'].append(smi); rows['score'].append(s23[len(rows['score'])])
df = pd.DataFrame(rows)
assert df.smiles.notna().all()
df['rankpct'] = df.score.rank(pct=True)
df['final_rank'] = df.rankpct.rank(ascending=False, method='first').astype(int)
df = df.sort_values('final_rank').reset_index(drop=True)
df.to_csv(R / 'module2_routed2605_theirAdapters.csv', index=False)
df.head(100).to_csv(R / 'module2_top100_theirAdapters.csv', index=False)

# 200 内排序对账
h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', BR, '--',
                    'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv'],
                   capture_output=True, text=True).stdout.split()
arch = pd.read_csv(io.BytesIO(subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True).stdout))
acanon = [canon(s) for s in arch.canonical_smiles]
amap = {c: i for i, c in enumerate(acanon)}
mcanon = [canon(s) for s in df.smiles]
mask = [c in amap for c in mcanon]
sub = df[mask].copy()
sub['arch_rank'] = [amap[c] for c, m in zip(mcanon, mask) if m]
from scipy.stats import spearmanr
print('200交集:', len(sub), ' Spearman(全库名次, 当年final_rank) = %.3f' %
      spearmanr(sub.final_rank, sub.arch_rank).statistic)
print('当年三候选在全库的名次:')
for pid in ('PACER0010', 'PACER0073', 'PACER0014', 'PACER0027'):
    i = arch[arch.pair_id == pid].index[0]
    c = acanon[i]
    row = df[[x == c for x in mcanon]]
    if len(row):
        print('  %s -> 全库名次 %d / 2605  新分 %.4f (当年 %.4f)' %
              (pid, row.final_rank.iloc[0], row.score.iloc[0], arch[arch.pair_id == pid].score_2023_gpcr_loto.iloc[0]))
print('新 top-5:', [(r.gen_id, round(r.score, 4)) for _, r in df.head(5).iterrows()])
