"""模块3 真跑（静态晶体臂）：top-200 × 7TRQ 别构口袋 Vina 对接。

协议（对齐本地冻结基准 dock_benchmark_chrm4.py）：
  受体 = 7TRQ_R_meeko.pdbqt（meeko 制备）
  盒子 = results/structure/allosteric_box.config（VU0467154 别构位，22 Å 立方）
  exhaustiveness = 8, num_modes = 9, seed = 42
输出：每分子 best Vina affinity + 排名合并表。
"""
import subprocess, sys, tempfile, json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

ROOT = Path(r'D:\CLC\project')
VINA = ROOT / 'tools' / 'vina.exe'
REC = ROOT / 'results' / 'structure' / 'ensemble' / '7TRQ_R_meeko.pdbqt'
TOP200 = ROOT / 'results' / 'pacer_rerun_v01' / 'module2_top200_output.csv'
OUT = ROOT / 'results' / 'pacer_rerun_v01' / 'module3_static7TRQ_results.csv'
BOX = dict(center_x=107.948, center_y=85.737, center_z=70.412,
           size_x=22.0, size_y=22.0, size_z=22.0)

def prep_ligand(smiles):
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(str(smiles)))
    if mol is None or AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        if mol is None or AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42) != 0:
            return None
    try:
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return None
        s, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        return s if ok else None
    except Exception:
        return None

def dock_one(args):
    idx, smiles, tmpdir = args
    pdbqt = prep_ligand(smiles)
    if pdbqt is None:
        return idx, None
    lig = Path(tmpdir) / f'lig{idx}.pdbqt'
    lig.write_text(pdbqt)
    outp = Path(tmpdir) / f'out{idx}.pdbqt'
    cmd = [str(VINA), '--receptor', str(REC), '--ligand', str(lig), '--out', str(outp),
           '--center_x', str(BOX['center_x']), '--center_y', str(BOX['center_y']), '--center_z', str(BOX['center_z']),
           '--size_x', str(BOX['size_x']), '--size_y', str(BOX['size_y']), '--size_z', str(BOX['size_z']),
           '--exhaustiveness', '8', '--num_modes', '9', '--seed', '42', '--cpu', '1']
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    best = None
    for line in r.stdout.splitlines():
        if line.strip().startswith('1 ') and '       ' in line:
            best = float(line.split()[1]); break
    for line in r.stdout.splitlines():
        if 'REMARK VINA RESULT' in line and best is None:
            best = float(line.split()[3]); break
    return idx, best

def main():
    df = pd.read_csv(TOP200)
    print('ligands:', len(df), flush=True)
    tmpdir = tempfile.mkdtemp(prefix='vina7trq_')
    jobs = [(i, s, tmpdir) for i, s in enumerate(df.smiles)]
    results = {}
    with ProcessPoolExecutor(max_workers=7) as ex:
        for n, (idx, best) in enumerate(ex.map(dock_one, jobs), 1):
            results[idx] = best
            if n % 25 == 0:
                print(f'{n}/{len(jobs)} done', flush=True)
    df['vina_7TRQ_best'] = [results.get(i) for i in range(len(df))]
    ok = df.vina_7TRQ_best.notna().sum()
    print(f'docked {ok}/{len(df)}')
    d2 = df.dropna(subset=['vina_7TRQ_best']).copy()
    d2['vina_rank'] = d2.vina_7TRQ_best.rank(method='first').astype(int)
    d2 = d2.sort_values('vina_rank')
    d2.to_csv(OUT, index=False)
    print('saved', OUT)
    print('Vina top-8:')
    for _, r in d2.head(8).iterrows():
        print('  %d  %.2f  模块2排名#%d  %s' % (r.vina_rank, r.vina_7TRQ_best, r.final_rank, str(r.smiles)[:46]))

if __name__ == '__main__':
    main()
