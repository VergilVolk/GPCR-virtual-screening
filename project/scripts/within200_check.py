import pandas as pd, subprocess, io
from rdkit import Chem
from rdkit.RDLogger import DisableLog
DisableLog('rdApp.*')
canon = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s))) or str(s)

mine = pd.read_csv(r'D:\CLC\project\results\pacer_rerun_v01\module2_routed2605_full.csv')
h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', 'origin/incident/pacer-xr-protocol-deviation-20261003', '--',
                    'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv'],
                   capture_output=True, text=True).stdout.split()
arch = pd.read_csv(io.BytesIO(subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True).stdout))
acanon = [canon(s) for s in arch.canonical_smiles]
mcanon = [canon(s) for s in mine.smiles]
mask = [c in set(acanon) for c in mcanon]
sub = mine[mask].copy()
sub['arch_canon'] = [c for c, m in zip(mcanon, mask) if m]
arch_idx = {c: i for i, c in enumerate(acanon)}
sub['my_rank_within200'] = range(1, len(sub) + 1)
sub['arch_final_rank'] = [arch_idx[c] for c in sub.arch_canon]
from scipy.stats import spearmanr
rho = spearmanr(sub.sort_values('my_rank_within200').arch_final_rank,
                range(1, len(sub) + 1)).statistic
print('200交集大小:', len(sub))
print('Spearman(我的200内排名, 当年final_rank) = %.3f' % rho)
top = sub.sort_values('my_rank_within200').head(8)
print('我的200内前八 -> 当年排名:')
for _, r in top.iterrows():
    pid = arch[arch.canonical_smiles.notna()]
    print('  我 #%d  当年 #%-3d  新分 %.4f' % (r.my_rank_within200, r.arch_final_rank, r.score_2023_gpcr_loto))
sub.to_csv(r'D:\CLC\project\results\pacer_rerun_v01\within200_ordering_check.csv', index=False)
