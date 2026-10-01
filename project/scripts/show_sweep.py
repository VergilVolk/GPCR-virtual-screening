import json, glob
def macro(path):
    d = json.load(open(path, encoding='utf-8'))
    m = d['metrics']['bce_retrieval']['macro']
    o = d['metrics']['official']['macro']
    return m['roc_auc'], m.get('pr_auc', m.get('average_precision', 0)), m.get('bedroc_alpha20', m.get('bedroc', 0)), o['roc_auc']

r, p, b, o = macro(r'D:\CLC\project\results\drugclip_science2026\migration_13target_v01\loso_preserve0_seed20260925.json')
print('baseline(seed25, rank8, ep40): ROC=%.4f PR=%.4f BEDROC=%.4f (official %.4f)' % (r, p, b, o))
for f in sorted(glob.glob(r'D:\CLC\project\results\drugclip_science2026\sweep_v01\*_seed20260925.json')):
    r, p, b, o = macro(f)
    name = f.split('\\')[-1].split('_')[0]
    print('%-10s ROC=%.4f PR=%.4f BEDROC=%.4f (official %.4f)' % (name, r, p, b, o))
