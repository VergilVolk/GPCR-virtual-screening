import subprocess, io, pandas as pd
from rdkit import Chem
from rdkit.RDLogger import DisableLog
DisableLog('rdApp.*')

def blob(path):
    h = subprocess.run(['git', '-C', r'D:\CLC', 'ls-tree', 'origin/incident/pacer-xr-protocol-deviation-20261003', '--', path],
                       capture_output=True, text=True).stdout.split()
    return subprocess.run(['git', '-C', r'D:\CLC', 'cat-file', 'blob', h[2]], capture_output=True).stdout

BR = 'origin/incident/pacer-xr-protocol-deviation-20261003'
arch = pd.read_csv(io.BytesIO(blob(
    'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv')))

gen = pd.read_csv(io.BytesIO(blob(
    'project/archive/incident_20261003/external_evidence/C/projects/GPCR-virtual-screening-final/project/results/generated/generated_pam_analogs.csv')))
print('gen rows:', len(gen), 'nan smiles:', gen.smiles.isna().sum())
gc = {}
for i, s in enumerate(gen.smiles.dropna()):
    m = Chem.MolFromSmiles(str(s))
    if m:
        gc.setdefault(Chem.MolToSmiles(m), i)

hit = 0
for pid in ('PACER0010', 'PACER0073', 'PACER0027', 'PACER0014'):
    smi = arch[arch.pair_id == pid].canonical_smiles.iloc[0]
    m = Chem.MolFromSmiles(str(smi))
    c = Chem.MolToSmiles(m)
    if c in gc:
        i = gc[c]
        print('%s = gen 行 %d  %s' % (pid, i, str(gen.smiles.iloc[i])[:60]))
        hit += 1
    else:
        # 就近找：与归档 SMILES 的 MCS 级相似
        from rdkit.Chem import rdFMCS
        print('%s 精确不匹配! 归档: %s' % (pid, str(smi)[:70]))
print('精确命中 %d / 4' % hit)
