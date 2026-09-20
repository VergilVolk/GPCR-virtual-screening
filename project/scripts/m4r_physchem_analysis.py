# -*- coding: utf-8 -*-
"""Zero-cost reproduction: Nguyen 2025 lipophilicity prior on the M4R open benchmark.

Compares physicochemical properties of M4R PAM actives (n=2314) vs property-matched
decoys (n~115k) from the Miao 2026 open data, and cross-checks against their Table 1
(median MW 385.684, logP 3.443, HBD 1, HBA 6, RB 4, rings 4). Also saves per-set
distributions for figures.
"""
import os, json, time
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

DATA = r'D:\CLC\project\data\miao2026_m4r'
OUT = r'D:\CLC\project\results\structure\m4r_physchem.json'

DESC = {
    'MW': lambda m: Descriptors.MolWt(m),
    'logP': lambda m: Crippen.MolLogP(m),
    'HBD': lambda m: rdMolDescriptors.CalcNumHBD(m),
    'HBA': lambda m: rdMolDescriptors.CalcNumHBA(m),
    'RB': lambda m: rdMolDescriptors.CalcNumRotatableBonds(m),
    'Rings': lambda m: rdMolDescriptors.CalcNumRings(m),
    'TPSA': lambda m: rdMolDescriptors.CalcTPSA(m),
    'Charge': lambda m: Chem.GetFormalCharge(m),
}

def load(smi_path, limit=None):
    mols = []
    with open(smi_path, encoding='utf-8') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if not parts:
                continue
            smi, cid = parts[0], parts[1] if len(parts) > 1 else str(i)
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            mols.append((cid, m))
            if limit and len(mols) >= limit:
                break
    return mols

def describe(mols, n_max=100000):
    rows = []
    for cid, m in mols[:n_max]:
        rows.append({k: fn(m) for k, fn in DESC.items()})
    return rows

def summarize(rows):
    import statistics as st
    out = {}
    for k in DESC:
        vals = [r[k] for r in rows]
        out[k] = {'median': round(st.median(vals), 3), 'mean': round(st.mean(vals), 3),
                  'p25': round(sorted(vals)[len(vals) // 4], 3), 'p75': round(sorted(vals)[3 * len(vals) // 4], 3)}
    return out

def main():
    t0 = time.time()
    act = load(os.path.join(DATA, 'M4R_actives.smi'))
    dec = load(os.path.join(DATA, 'M4R_decoys.smi'), limit=20000)  # sample for speed
    print(f'actives parsed: {len(act)}, decoys parsed (sample 20k): {len(dec)}', flush=True)
    a_rows = describe(act)
    d_rows = describe(dec)
    sa, sd = summarize(a_rows), summarize(d_rows)
    res = {
        'n_actives': len(a_rows), 'n_decoys_sampled': len(d_rows),
        'actives': sa, 'decoys': sd,
        'diff_median': {k: round(sa[k]['median'] - sd[k]['median'], 3) for k in DESC},
        'note': 'Nguyen2025 A1R: PAMs more lipophilic -> LogP 2.5-4.5 filter. Test on M4R.',
        'paper_table1_medians': {'MW': 385.684, 'logP': 3.443, 'HBD': 1, 'HBA': 6, 'RB': 4, 'Rings': 4},
        'runtime_s': round(time.time() - t0, 1),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(json.dumps(res, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
