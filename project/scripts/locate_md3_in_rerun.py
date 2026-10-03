import pandas as pd, subprocess, io
from rdkit import Chem

def canon(s):
    try:
        m = Chem.MolFromSmiles(str(s))
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None

h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', 'origin/incident/pacer-xr-protocol-deviation-20261003', '--',
                    'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv'],
                   capture_output=True, text=True).stdout.split()
raw = subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True).stdout
arch = pd.read_csv(io.BytesIO(raw))
targets = {}
for pid, s in zip(arch.pair_id, arch.canonical_smiles):
    if pid in ('PACER0010', 'PACER0073', 'PACER0027', 'PACER0014'):
        c = canon(s)
        if c:
            targets[c] = pid

mine = pd.read_csv(r'D:\CLC\project\results\pacer_rerun_v01\module2_routed2605_full.csv')
mine['canon'] = [canon(s) for s in mine.smiles]
for cs, pid in targets.items():
    row = mine[mine.canon == cs]
    if len(row):
        r = row.iloc[0]
        print('%s -> 新分数 %.4f  新排名 %d / 2605' % (pid, r.score_2023_gpcr_loto, r.final_rank))
    else:
        print(pid, '不在 2605 里!')
