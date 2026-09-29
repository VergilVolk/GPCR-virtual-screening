import json, glob
import numpy as np
out = r'D:\CLC\project\results\drugclip_science2026\migration_13target_v01'
names = ['roc_auc', 'pr_auc', 'bedroc_alpha20', 'ef1pct', 'ef5pct']

# ECFP structure
e = json.load(open(out + r'\ecfp_loso.json', encoding='utf-8'))
print('ECFP metrics keys:', list(e['metrics'].keys()))
for k, v in e['metrics'].items():
    if isinstance(v, dict) and 'macro' in v:
        print('ECFP macro:', {n: round(v['macro'][n], 4) for n in names})
        break

# per-target mean across 3 no-preserve seeds
per = {}
for f in sorted(glob.glob(out + r'\loso_preserve0_seed*.json')):
    d = json.load(open(f, encoding='utf-8'))
    for t in d['metrics']['bce_retrieval']['per_target']:
        per.setdefault(t, {'official': [], 'tuned': []})
        for k in ['official', 'bce_retrieval']:
            per[t]['official' if k == 'official' else 'tuned'].append(
                [d['metrics'][k]['per_target'][t][n] for n in names])
print(f"\n{'target':10s} {'off ROC':>8s} {'tun ROC':>8s} {'off EF1%':>9s} {'tun EF1%':>9s} {'off BED':>8s} {'tun BED':>8s}")
for t in per:
    o = np.mean(per[t]['official'], axis=0); u = np.mean(per[t]['tuned'], axis=0)
    print(f"{t:10s} {o[0]:8.3f} {u[0]:8.3f} {o[3]:9.2f} {u[3]:9.2f} {o[2]:8.3f} {u[2]:8.3f}")
