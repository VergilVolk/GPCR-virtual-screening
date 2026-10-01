import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

old_p = r'D:\CLC\project\results\litpcba_drugclip_external_v01\combined_13target_loso_3seed_summary.predictions.csv'
new_dir = r'D:\CLC\project\results\drugclip_science2026\sweep_v01'
out_p = r'D:\CLC\project\results\drugclip_science2026\xcheckpoint_ensemble_v01\xens_3seed.json'

old = pd.read_csv(old_p)
keys = ['target', 'label', 'canonical_smiles', 'murcko_scaffold']
tabs = [pd.read_csv(new_dir + r'\epochs80_seed%s.predictions.csv' % s) for s in (20260925, 20260926, 20260927)]
new = tabs[0][keys + ['official']].copy()
new['tuned'] = np.mean([t['tuned'].values for t in tabs], axis=0)
m = old[keys + ['tuned']].rename(columns={'tuned': 'old_tuned'}).merge(
    new[keys + ['tuned', 'official']].rename(columns={'tuned': 'new_tuned'}), on=keys, how='inner')
m['cluster'] = m.target + '|' + m.murcko_scaffold.astype(str)
targets = sorted(m.target.unique())
rng = np.random.default_rng(20260926)

def per_target_scores(col):
    return {t: roc_auc_score(s.label, s[col]) for t, s in m.groupby('target') if s.label.nunique() == 2}

def fast_auc(y, s):
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if not n1 or not n0:
        return None
    r = np.argsort(np.argsort(s))
    return (r[y == 1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)

def macro_of(col, tab):
    rows = []
    for t in targets:
        sub = tab[tab.target == t]
        y, s = sub.label.values, sub[col].values
        order = np.argsort(-s, kind='mergesort')
        k = max(1, len(y) // 100)
        rows.append((roc_auc_score(y, s), average_precision_score(y, s), y[order[:k]].mean() / y.mean()))
    return tuple(np.mean(rows, axis=0))

# 主配置 w=0.5（对称，无挑权重嫌疑）；敏感性 0.4/0.6
r_old, r_new = m.old_tuned.rank(), m.new_tuned.rank()
m['xens05'] = 0.5 * r_old + 0.5 * r_new
m['xens04'] = 0.4 * r_old + 0.6 * r_new
m['xens06'] = 0.6 * r_old + 0.4 * r_new

report = {'protocol': 'cross-checkpoint ensemble: rank-mean of old-2023 3-seed LOSO tuned '
                      'and Science-2026 ep80 3-seed ensemble tuned scores; identical 37,524 rows',
          'n_pairs': int(len(m)), 'n_targets': len(targets)}
for c in ('old_tuned', 'new_tuned', 'official', 'xens04', 'xens05', 'xens06'):
    r, p, e = macro_of(c, m)
    report[c] = {'roc_auc': r, 'pr_auc': p, 'ef1pct': e}
    print('%-10s ROC=%.4f PR=%.4f EF1%%=%.2f' % (c, r, p, e))

# 多样性证据：两路分数的逐靶 Spearman
from scipy.stats import spearmanr
rhos = [spearmanr(s.old_tuned, s.new_tuned).statistic for _, s in m.groupby('target')]
print('per-target spearman(old,new): median=%.3f min=%.3f max=%.3f' % (np.median(rhos), np.min(rhos), np.max(rhos)))
report['diversity'] = {'per_target_spearman_median': float(np.median(rhos)),
                       'per_target_spearman_min': float(np.min(rhos)),
                       'per_target_spearman_max': float(np.max(rhos))}

# 靶分层 cluster bootstrap：xens05 - official 与 xens05 - new_tuned
y_arr = m.label.values.astype(int)
tgt_arr = m.target.values
cols = {'xens05': m.xens05.values, 'official': m.official.values, 'new_tuned': m.new_tuned.values}
m['cluster'] = m.target + '|' + m.murcko_scaffold.astype(str)
clusters = sorted(m.cluster.unique())
cl_rows = {c: np.flatnonzero((m.cluster == c).values) for c in clusters}
tgt_clusters = {t: np.unique(m.cluster.values[tgt_arr == t]) for t in targets}

def boot_ci(a, b, n=2000):
    d = []
    for _ in range(n):
        aa, bb = [], []
        for t, cls in tgt_clusters.items():
            pick = rng.choice(cls, size=len(cls), replace=True)
            sel = np.concatenate([cl_rows[c] for c in pick])
            ya = y_arr[sel]
            if ya.sum() == 0 or ya.min() == ya.max():
                continue
            ua, ub = fast_auc(ya, a[sel]), fast_auc(ya, b[sel])
            if ua is not None:
                aa.append(ua); bb.append(ub)
        if aa:
            d.append(np.mean(aa) - np.mean(bb))
    d = np.array(d)
    return [float(np.quantile(d, .025)), float(np.quantile(d, .975))]

report['bootstrap_95ci'] = {
    'xens05_minus_official_roc': boot_ci(cols['xens05'], cols['official']),
    'xens05_minus_newens_roc': boot_ci(cols['xens05'], cols['new_tuned']),
}
print('CI xens-official ROC:', report['bootstrap_95ci']['xens05_minus_official_roc'])
print('CI xens-newens  ROC:', report['bootstrap_95ci']['xens05_minus_newens_roc'])

import os
os.makedirs(os.path.dirname(out_p), exist_ok=True)
json.dump(report, open(out_p, 'w', encoding='utf-8'), indent=2)
m.to_csv(out_p.replace('.json', '.predictions.csv'), index=False)
print('saved:', out_p)
