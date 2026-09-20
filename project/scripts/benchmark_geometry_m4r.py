# -*- coding: utf-8 -*-
"""DEPRECATED / INVALID exploratory benchmark.

The historical implementation associates a hub contact with the last receptor
atom visited rather than the atom defining the minimum distance. Its saved AUC
must not be cited. Use dock_potency_coupling_features.py plus
evaluate_structure_loso.py for the audited experiment.

Large-n benchmark: our geometric contact evidence on the Miao 2026 M4R open data.

Reproduces the M4R screening benchmark (2314 PAMs / ~115k decoys) and scores it with
OUR geometric features instead of (or in addition to) force-field scoring:

  pipeline:
    subsample actives+decoys
      -> 3D embed (RDKit ETKDG) -> meeko PDBQT
      -> pose generation into receptor conformation(s) [Vina as POSE GENERATOR only;
         scores are never used for ranking]
      -> geometric features per pose:
            n_contacts   : protein-ligand heavy-atom contacts (<4.5 A)
            hub_contacts : contacts with the allosteric hub residue (TYR439 in 7TRQ numbering)
            buried_ratio : fraction of ligand heavy atoms within 6 A of protein
            pocket_density: contacts per ligand heavy atom
      -> rank by each feature; report AUC / EF(0.5%) / EF(1%) vs is_active
      -> compare with published Vina PDB / BEmin / BEavg baselines

Usage:
  python benchmark_geometry_m4r.py --n-act 200 --n-dec 200 --receptor 7TRQ
"""
import argparse, csv, os, sys, json, subprocess, time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT = r'D:\CLC\project'
DATA = os.path.join(ROOT, 'data', 'miao2026_m4r')
OUT = os.path.join(ROOT, 'results', 'structure')
VINA = os.path.join(ROOT, 'tools', 'vina.exe')
SEED = 42

# ---- pocket definition (7TRQ numbering, validated: contact residues 13/13) ----
# allosteric pocket residues from validate_allosteric_pocket.py (13 residues)
POCKET_RESIDUES = []  # filled at runtime from receptor PDBQT or passed explicitly

def load_smiles(path, n=None):
    out = []
    with open(path, encoding='utf-8') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if not parts:
                continue
            out.append((parts[1] if len(parts) > 1 else str(i), parts[0]))
            if n and len(out) >= n:
                break
    return out

def embed_smiles(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    return m

def mol_to_pdbqt(mol, out_pdbqt):
    """meeko PDBQTWriterLegacy (same as prep_receptor_7trq.py)."""
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)
    with open(out_pdbqt, 'w') as fh:
        for setup in mol_setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup)
            if not is_ok:
                return False
            fh.write(pdbqt_string)
    return True

def parse_pdbqt_model(pdbqt_path):
    """Return list of (atom_name, resnum, resname, x, y, z, element) for MODEL 1."""
    atoms = []
    in_model = False
    with open(pdbqt_path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('MODEL'):
                in_model = True
                continue
            if line.startswith('ENDMDL'):
                break
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append({
                    'name': line[12:16].strip(), 'res': int(line[22:26]),
                    'resname': line[17:20].strip(),
                    'x': float(line[30:38]), 'y': float(line[38:46]), 'z': float(line[46:54]),
                    'elem': line[76:78].strip().upper(),
                })
    return atoms

def load_receptor_pocket(path, radius_around_center=None, center=None):
    """Load receptor atoms; return (atoms, pocket_atom_indices)."""
    atoms = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append({
                    'name': line[12:16].strip(), 'res': int(line[22:26]),
                    'resname': line[17:20].strip(),
                    'x': float(line[30:38]), 'y': float(line[38:46]), 'z': float(line[46:54]),
                    'elem': line[76:78].strip().upper(),
                })
    if center is not None and radius_around_center is not None:
        idx = []
        for i, a in enumerate(atoms):
            d = ((a['x']-center[0])**2 + (a['y']-center[1])**2 + (a['z']-center[2])**2) ** 0.5
            if d <= radius_around_center:
                idx.append(i)
        return atoms, idx
    return atoms, list(range(len(atoms)))

def geometric_features(lig_atoms, rec_atoms, pocket_idx, hub_residues=()):
    """Compute geometric contact features for one pose."""
    n_contacts = 0
    hub_contacts = 0
    buried = 0
    contact_dists = []
    for la in lig_atoms:
        min_d = 1e9
        for ri in pocket_idx:
            ra = rec_atoms[ri]
            d = ((la['x']-ra['x'])**2 + (la['y']-ra['y'])**2 + (la['z']-ra['z'])**2) ** 0.5
            if d < min_d:
                min_d = d
        if min_d < 4.5:
            n_contacts += 1
            contact_dists.append(min_d)
            if ra['res'] in hub_residues:
                hub_contacts += 1
        if min_d < 6.0:
            buried += 1
    n = max(len(lig_atoms), 1)
    return {
        'n_contacts': n_contacts,
        'hub_contacts': hub_contacts,
        'buried_ratio': buried / n,
        'pocket_density': n_contacts / n,
    }

def run_vina(lig_pdbqt, rec_pdbqt, center, box, out_pdbqt, seed=SEED):
    cfg = os.path.join(OUT, '_tmp_bench.cfg')
    with open(cfg, 'w') as f:
        f.write(f'receptor = {rec_pdbqt}\nligand = {lig_pdbqt}\n')
        f.write(f'center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\ncenter_z = {center[2]:.3f}\n')
        f.write(f'size_x = {box[0]}\nsize_y = {box[1]}\nsize_z = {box[2]}\n')
        f.write(f'exhaustiveness = 8\nnum_modes = 1\nseed = {seed}\n')
        f.write(f'out = {out_pdbqt}\n')
    r = subprocess.run([VINA, '--config', cfg], capture_output=True, text=True)
    return r.returncode

def auc_rank(scores, labels):
    """AUC with higher score = better (for ranking we use feature magnitude)."""
    order = np.argsort(scores)
    ranks = np.empty(len(scores))
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and scores[order[j+1]] == scores[order[i]]:
            j += 1
        ranks[order[i:j+1]] = (i + j) / 2.0 + 1
        i = j + 1
    r1 = ranks[np.array(labels) == 1].sum()
    n1 = sum(labels); n0 = len(labels) - n1
    return (r1 - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else 0.5

def ef_top(scores, labels, pct):
    n = len(scores)
    k = max(1, int(round(pct / 100.0 * n)))
    order = np.argsort(scores)[::-1][:k]
    hits = sum(1 for i in order if labels[i])
    total_hits = sum(labels)
    return hits / k / (total_hits / n) if total_hits else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-act', type=int, default=200)
    ap.add_argument('--n-dec', type=int, default=200)
    ap.add_argument('--rec', default='7TRQ')  # receptor source for poses
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    actives = load_smiles(os.path.join(DATA, 'M4R_actives.smi'), args.n_act)
    decoys = load_smiles(os.path.join(DATA, 'M4R_decoys.smi'), args.n_dec)
    print(f'actives: {len(actives)}, decoys: {len(decoys)}', flush=True)

    # receptor + pocket (7TRQ validated allosteric box center/size)
    rec_pdbqt = os.path.join(ROOT, 'results', 'structure', '7trq_R_receptor.pdbqt')
    if not os.path.exists(rec_pdbqt):
        # fall back to .pdb (should have been converted in prep)
        rec_pdbqt = os.path.join(ROOT, 'results', 'structure', '7trq_R_receptor.pdb')
    center = (107.95, 85.74, 70.41)  # allosteric box from prep_receptor_7trq.py
    box = (22, 22, 22)
    hub_residues = {439}  # TYR439 hub in 7TRQ numbering

    rec_atoms, pocket_idx = load_receptor_pocket(rec_pdbqt, radius_around_center=8.0, center=center)
    print(f'receptor atoms: {len(rec_atoms)}, pocket atoms within 8A: {len(pocket_idx)}', flush=True)

    tmp = os.path.join(OUT, '_tmp_bench')
    os.makedirs(tmp, exist_ok=True)
    rows = []
    t0 = time.time()
    for i, (cid, smi) in enumerate(actives + decoys):
        active = 1 if i < len(actives) else 0
        mol = embed_smiles(smi)
        if mol is None:
            continue
        lig_pdbqt = os.path.join(tmp, f'lig_{i}.pdbqt')
        out_pdbqt = os.path.join(tmp, f'pose_{i}.pdbqt')
        if not mol_to_pdbqt(mol, lig_pdbqt):
            continue
        rc = run_vina(lig_pdbqt, rec_pdbqt, center, box, out_pdbqt)
        if rc != 0 or not os.path.exists(out_pdbqt):
            continue
        lig_atoms = parse_pdbqt_model(out_pdbqt)
        feats = geometric_features(lig_atoms, rec_atoms, pocket_idx, hub_residues)
        feats.update({'id': cid, 'active': active, 'smiles': smi})
        rows.append(feats)
        if (i + 1) % 50 == 0:
            print(f'  {i+1}/{len(actives)+len(decoys)} done, elapsed {time.time()-t0:.0f}s', flush=True)

    # metrics per feature
    res = {'n': len(rows), 'n_active': sum(1 for r in rows if r['active']),
           'seed': args.seed, 'features': {}}
    for feat in ['n_contacts', 'hub_contacts', 'buried_ratio', 'pocket_density']:
        scores = [r[feat] for r in rows]
        labels = [r['active'] for r in rows]
        res['features'][feat] = {
            'AUC': round(auc_rank(scores, labels), 3),
            'EF(0.5%)': round(ef_top(scores, labels, 0.5), 3),
            'EF(1%)': round(ef_top(scores, labels, 1.0), 3),
        }
    with open(os.path.join(OUT, 'benchmark_geometry_m4r.json'), 'w', encoding='utf-8') as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(json.dumps(res, ensure_ascii=False, indent=2))
    print('Pose-generation engine used only for poses; ranking used geometric features only.')

if __name__ == '__main__':
    main()
