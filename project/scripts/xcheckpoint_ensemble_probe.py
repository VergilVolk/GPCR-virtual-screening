import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

old_p = r'D:\CLC\project\results\litpcba_drugclip_external_v01\combined_13target_loso_3seed_summary.predictions.csv'
new_dir = r'D:\CLC\project\results\drugclip_science2026\sweep_v01'
old = pd.read_csv(old_p)
print('old cols:', list(old.columns)[:8])
keys = ['target', 'label', 'canonical_smiles', 'murcko_scaffold']
tabs = [pd.read_csv(new_dir + r'\epochs80_seed%s.predictions.csv' % s) for s in (20260925, 20260926, 20260927)]
new = tabs[0][keys + ['official']].copy()
new['tuned'] = np.mean([t['tuned'].values for t in tabs], axis=0)
m = old[keys + ['tuned']].rename(columns={'tuned': 'old_tuned'}).merge(
    new[keys + ['tuned', 'official']].rename(columns={'tuned': 'new_tuned'}), on=keys, how='inner')
print('joined rows:', len(m), ' old:', len(old), ' new:', len(new))
targets = sorted(m.target.unique())

def macro(col):
    rows = []
    for t in targets:
        sub = m[m.target == t]
        y, s = sub.label.values, sub[col].values
        order = np.argsort(-s, kind='mergesort')
        k = max(1, len(y) // 100)
        rows.append((roc_auc_score(y, s), average_precision_score(y, s), y[order[:k]].mean() / y.mean()))
    return tuple(np.mean(rows, axis=0))

for c in ('old_tuned', 'new_tuned', 'official'):
    print('%-10s ROC=%.4f PR=%.4f EF1%%=%.2f' % (c, *macro(c)))
for w in (0.4, 0.5, 0.6):
    r_old = m.old_tuned.rank(); r_new = m.new_tuned.rank()
    m['xens'] = w * r_old + (1 - w) * r_new
    print('xens w=%.1f  ROC=%.4f PR=%.4f EF1%%=%.2f' % (w, *macro('xens')))
