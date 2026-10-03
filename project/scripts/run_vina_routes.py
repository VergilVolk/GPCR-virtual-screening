"""通用 Vina 路由 runner：共享配体目录 + 指定受体。"""
import subprocess, sys, tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

ROOT = Path(r'D:\CLC\project')
VINA = ROOT / 'tools' / 'vina.exe'
LIGDIR = ROOT / 'results' / 'pacer_rerun_v01' / 'ligands'
TOP200 = ROOT / 'results' / 'pacer_rerun_v01' / 'module2_top200_output.csv'

def prep_all():
    LIGDIR.mkdir(parents=True, exist_ok=True)
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    df = pd.read_csv(TOP200)
    n = 0
    for i, smi in enumerate(df.smiles):
        p = LIGDIR / f'lig{i:04d}.pdbqt'
        if p.exists():
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(str(smi)))
        if mol is None or AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            if mol is None or AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42) != 0:
                continue
        try:
            setups = MoleculePreparation().prepare(mol)
            if setups:
                s, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
                if ok:
                    p.write_text(s); n += 1
        except Exception:
            pass
    print(f'ligands prepared (+{n})', flush=True)

def dock_one(args):
    idx, rec, tmpdir = args
    lig = LIGDIR / f'lig{idx:04d}.pdbqt'
    if not lig.exists():
        return idx, None
    outp = Path(tmpdir) / f'out{idx:04d}.pdbqt'
    cmd = [str(VINA), '--receptor', str(rec), '--ligand', str(lig), '--out', str(outp),
           '--center_x', '107.948', '--center_y', '85.737', '--center_z', '70.412',
           '--size_x', '22.0', '--size_y', '22.0', '--size_z', '22.0',
           '--exhaustiveness', '8', '--num_modes', '9', '--seed', '42', '--cpu', '1']
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    # vina 结果表打 stdout，REMARK 只写 out 文件：两处都解析
    for line in r.stdout.splitlines():
        if 'REMARK VINA RESULT' in line:
            return idx, float(line.split()[3])
    for line in r.stdout.splitlines():
        ls = line.strip()
        if ls and ls[0].isdigit() and len(ls.split()) >= 2:
            try:
                return idx, float(ls.split()[1])
            except ValueError:
                pass
    if outp.exists():
        for line in outp.read_text(errors='replace').splitlines():
            if 'REMARK VINA RESULT' in line:
                return idx, float(line.split()[3])
    return idx, None

def run(rec_name):
    rec = ROOT / 'results' / 'structure' / 'ensemble' / f'{rec_name}_R_meeko.pdbqt'
    out_csv = ROOT / 'results' / 'pacer_rerun_v01' / f'module3_{rec_name}_results.csv'
    if out_csv.exists():
        print(rec_name, 'already done'); return
    prep_all()
    df = pd.read_csv(TOP200)
    tmpdir = tempfile.mkdtemp(prefix=f'vina{rec_name}_')
    jobs = [(i, str(rec), tmpdir) for i in range(len(df))]
    results = {}
    with ProcessPoolExecutor(max_workers=7) as ex:
        for n, (idx, best) in enumerate(ex.map(dock_one, jobs), 1):
            results[idx] = best
            if n % 50 == 0:
                print(f'{rec_name} {n}/{len(jobs)}', flush=True)
    df[f'vina_{rec_name}'] = [results.get(i) for i in range(len(df))]
    df.to_csv(out_csv, index=False)
    ok = df[f'vina_{rec_name}'].notna().sum()
    print(f'{rec_name} DONE {ok}/{len(df)} -> {out_csv.name}', flush=True)

if __name__ == '__main__':
    for name in sys.argv[1:]:
        run(name)
