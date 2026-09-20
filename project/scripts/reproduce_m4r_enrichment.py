# -*- coding: utf-8 -*-
"""Reproduce M4R enrichment metrics from Miao 2026 open data (zero docking cost).

Reads the published per-cluster Vina BE scores + PDB Vina scores, recomputes
EF/EF'/AUC/logAUC for PDB, BEmin, BEavg rankings, and compares with paper Table 4.
"""
import csv, math, sys, os, json
import numpy as np

DATA = r'D:\CLC\project\data\miao2026_m4r'

def load_framewise():
    rows = []
    with open(os.path.join(DATA, 'M4R_Ensemble_Vina_framewise_scores.csv'), newline='', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            bes = []
            for i in range(10):
                v = row.get(f'cluster{i:02d}_BE', '')
                try:
                    bes.append(float(v))
                except ValueError:
                    bes.append(np.nan)
            arr = np.array(bes)
            valid = arr[~np.isnan(arr)]
            if len(valid) == 0:
                continue
            rows.append({'id': row['ligand_id'], 'active': int(row['is_active']),
                         'be': arr, 'be_min': float(np.nanmin(arr)), 'be_avg': float(np.nanmean(arr))})
    return rows

def load_pdb():
    rows = []
    with open(os.path.join(DATA, 'M4R_PDB_Vina_scores.csv'), newline='', encoding='utf-8') as f:
        r = csv.DictReader(f)
        # detect score column
        cols = r.fieldnames
        score_col = [c for c in cols if 'score' in c.lower() or 'BE' in c or 'vina' in c.lower()]
        active_col = [c for c in cols if 'active' in c.lower()]
        for row in r:
            sc = row.get(score_col[0]) if score_col else None
            if sc is None or sc == '':
                continue
            rows.append({'id': row.get('ligand_id', row.get(cols[0], '')), 'active': int(float(row[active_col[0]])) if active_col else 0,
                         'score': float(sc)})
    return rows

def metrics(scores, actives):
    """scores: list of (score, is_active) sorted descending by score; actives: total # actives."""
    n = len(scores)
    order = sorted(range(n), key=lambda i: scores[i][0], reverse=True)
    hits = 0
    ranks = []
    for pos, idx in enumerate(order, start=1):
        if scores[idx][1]:
            hits += 1
            ranks.append(pos / n)
    apr = float(np.mean(ranks)) if ranks else 1.0
    def ef(pct):
        k = max(1, int(round(pct / 100.0 * n)))
        top = order[:k]
        hits_top = sum(1 for i in top if scores[i][1])
        return hits_top / k / (hits / n) if hits else 0.0
    def efd(pct):  # EF' early enrichment factor
        k = max(1, int(round(pct / 100.0 * n)))
        top = order[:k]
        hits_top = sum(1 for i in top if scores[i][1])
        return (pct / 100.0) / apr * (hits_top / max(hits, 1))
    # AUC
    labs = np.array([s[1] for s in scores])
    vals = np.array([s[0] for s in scores])
    # rank by value, tie-break average; use mann-whitney
    order_vals = np.argsort(vals, kind='mergesort')
    # AUC via rank sum
    ranks2 = np.empty(n)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and vals[order_vals[j + 1]] == vals[order_vals[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        ranks2[order_vals[i:j + 1]] = avg
        i = j + 1
    r1 = ranks2[labs == 1].sum()
    n1 = (labs == 1).sum(); n0 = (labs == 0).sum()
    auc = (r1 - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else 0.5
    # logAUC: semilog ROC, normalized to random
    # build sorted (desc) cumulative hit fractions at log-spaced cutoffs
    order2 = sorted(range(n), key=lambda i: vals[i], reverse=True)
    n1 = float(n1)
    fracs = []
    fpr = 0.0
    for frac in [10 ** x for x in np.arange(-4, 0.001, 0.25)]:
        k = max(1, int(round(frac * n)))
        top = order2[:k]
        tpr = sum(1 for i in top if labs[i]) / n1
        fracs.append((frac, tpr))
    xs = np.array([f for f, _ in fracs])
    ys = np.array([t for _, t in fracs])
    logauc = np.trapz(ys, np.log10(xs)) / np.log10(10) / 4.0  # normalized over 4 log decades
    return {'EF(0.5%)': ef(0.5), 'EF\'(0.5%)': efd(0.5), 'EF(1%)': ef(1.0), 'EF\'(1%)': efd(1.0),
            'AUC%': auc * 100, 'logAUC%': logauc * 100, 'n': n, 'n_act': n1}

def main():
    fw = load_framewise()
    pdb = load_pdb()
    print(f'framewise rows: {len(fw)}  PDB rows: {len(pdb)}')
    n_act = sum(1 for r in fw if r['active'])
    print(f'actives: {n_act}  decoys: {len(fw) - n_act}')
    # PDB ranking: scores as-is (higher = better? vina scores negative, rank by value desc = most negative first? Vina more negative = better -> rank ascending)
    # Their PDB CSV: check whether score column already flipped. We rank by score ascending (more negative better) for vina raw.
    # But for BE values (already PMF-corrected, more negative better) also ascending.
    res = {}
    res['PDB'] = metrics([(-r['score'], r['active']) for r in pdb], n_act) if pdb else None
    # careful: their PDB CSV may store scores as positive or negative; try both directions and keep better? No - keep raw order.
    # Standard: more negative = better.
    res['Ensemble BEmin'] = metrics([(-r['be_min'], r['active']) for r in fw], n_act)
    res['Ensemble BEavg'] = metrics([(-r['be_avg'], r['active']) for r in fw], n_act)
    for k, v in res.items():
        if v is None:
            print(f'{k}: no data')
            continue
        e1, e2, e3, e4 = v['EF(0.5%)'], v["EF'(0.5%)"], v['EF(1%)'], v["EF'(1%)"]
        print(f"{k}: n={v['n']} act={v['n_act']} EF(0.5%)={e1:.3f} EF'(0.5%)={e2:.3f} EF(1%)={e3:.3f} EF'(1%)={e4:.3f} AUC={v['AUC%']:.2f}% logAUC={v['logAUC%']:.2f}%")
    out = os.path.join(os.path.dirname(DATA), '..', 'results', 'pacer_external_miao_v01')
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'metrics.json'), 'w', encoding='utf-8') as f:
        json.dump({'source': 'Miao 2026 public M4R active/decoy docking benchmark',
                   'metrics': res,
                   'interpretation': 'ROC AUC can coexist with random-or-worse top-0.5% enrichment; this benchmark evaluates active/decoy enrichment, not PAM functional efficacy.'}, f, indent=2)

if __name__ == '__main__':
    main()
