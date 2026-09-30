import json, numpy as np, pandas as pd

out = r'D:\CLC\project\results\drugclip_science2026\extended20_finetune_v01'
ens = json.load(open(out + r'\ensemble20_ep80.json', encoding='utf-8'))
ecfp = json.load(open(out + r'\ecfp20_loso.json', encoding='utf-8'))
print('20-target macro:')
print('  ours(ens)  ROC=%.4f PR=%.4f BEDROC=%.4f EF1%%=%.2f' % (ens['ensemble_macro']['roc_auc'], ens['ensemble_macro']['pr_auc'], ens['ensemble_macro']['bedroc_alpha20'], ens['ensemble_macro']['ef1pct']))
print('  official   ROC=%.4f' % ens['official_macro']['roc_auc'])
m = ecfp['metrics']['macro']
print('  ECFP4      ROC=%.4f PR=%.4f BEDROC=%.4f EF1%%=%.2f' % (m['roc_auc'], m['pr_auc'], m['bedroc_alpha20'], m['ef1pct']))

# fusion probe: rank-average ours (ensemble mean) + ECFP on identical rows
files = [out + r'\loso20_ep80_seed20260925.predictions.csv',
         out + r'\loso20_ep80_seed20260926.predictions.csv',
         out + r'\loso20_ep80_seed20260927.predictions.csv']
tabs = [pd.read_csv(f) for f in files]
keys = ['target', 'label', 'canonical_smiles', 'murcko_scaffold']
t = tabs[0][keys + ['official']].copy()
t['tuned'] = np.mean([x['tuned'].values for x in tabs], axis=0)
e = pd.read_csv(out + r'\ecfp20_loso.predictions.csv')
m2 = t.merge(e[keys + ['ecfp4_logistic']], on=keys, how='inner')
print('merged rows:', len(m2), '/', len(t))
targets = sorted(m2.target.unique())

def metrics_row(sub, col):
    from sklearn.metrics import roc_auc_score, average_precision_score
    y, s = sub.label.values, sub[col].values
    order = np.argsort(-s, kind='mergesort')
    k = max(1, len(y) // 100)
    return (roc_auc_score(y, s), average_precision_score(y, s),
            y[order[:k]].mean() / y.mean())

def rank_avg(sub, cols):
    r = None
    for c in cols:
        rr = sub[c].rank().values
        r = rr if r is None else r + rr
    return r / len(cols)

res = {}
for w in (0.3, 0.5, 0.7):
    rows = []
    for tg in targets:
        sub = m2[m2.target == tg]
        # weighted rank fusion: w on ours, (1-w) on ECFP
        r_ours = sub['tuned'].rank().values
        r_ecfp = sub['ecfp4_logistic'].rank().values
        score = w * r_ours + (1 - w) * r_ecfp
        y = sub.label.values
        order = np.argsort(-score, kind='mergesort')
        k = max(1, len(y) // 100)
        from sklearn.metrics import roc_auc_score, average_precision_score
        rows.append((roc_auc_score(y, score), average_precision_score(y, score), y[order[:k]].mean() / y.mean()))
    res[w] = tuple(np.mean(rows, axis=0))
    print('fusion w=%.1f  ROC=%.4f PR=%.4f EF1%%=%.2f' % (w, *res[w]))
# pure baselines on same merged rows
for col in ('tuned', 'ecfp4_logistic', 'official'):
    rows = [metrics_row(m2[m2.target == tg], col) for tg in targets]
    print('%-13s ROC=%.4f PR=%.4f EF1%%=%.2f' % (col, *np.mean(rows, axis=0)))
