import numpy as np, pandas as pd, torch, sys, json
sys.path.insert(0, r'D:\CLC\project\scripts')
from finetune_drugclip_muscarinic_triplet import Proj
from finetune_drugclip_gpcr_screening import binary_metrics

# 同一测试集：GaMD M4R 筛选对（screening_pairs.csv 中 target==M4R 的行）
z = np.load(r'D:\CLC\project\results\gpcr_drugclip_screening_v01\ensemble_embeddings.npz', allow_pickle=False)
ids = [str(x) for x in z['molecule_ids']]; mi = {v: i for i, v in enumerate(ids)}
gp = [str(x) for x in z['pocket_ids']]
m4 = [i for i, p in enumerate(gp) if p.startswith('M4R')]
mol = torch.as_tensor(z['molecule_representations'].astype(np.float32))
pocket = torch.as_tensor(z['pocket_representations'][m4].astype(np.float32))
pairs = pd.read_csv(r'D:\CLC\project\data\benchmarks\gpcr_drugclip_screening_v01\screening_pairs.csv')
m4p = pairs[pairs.target == 'M4R'].copy()
m4p = m4p[m4p.canonical_smiles.astype(str).isin(mi)]
m4p['idx'] = [mi[v] for v in m4p.canonical_smiles.astype(str)]
y = m4p.label.to_numpy(int)
rows = m4p.idx.to_numpy()
print('M4R test set: %d rows, %d actives' % (len(y), y.sum()))

res = {}
# --- 2023 线：triplet external_result x3（训练于 ChEMBL 毒蕈碱族，未见此库）---
d0 = torch.load(r'D:\CLC\project\results\muscarinic_triplet_expansion_v01\external_result.seed20260924.projection.pt',
                map_location='cpu', weights_only=False)
pm0, pp0 = Proj(d0['mol_project']).eval(), Proj(d0['pocket_project']).eval()
with torch.inference_mode():
    sc0 = (pm0(mol[rows]) @ pp0(pocket).T).numpy()
c0 = gp.index('M4R_cluster0')
pocket_c0 = torch.as_tensor(z['pocket_representations'][[c0]].astype(np.float32))
with torch.inference_mode():
    sc_c0 = (pm0(mol[rows]) @ pp0(pocket_c0).T).numpy()[:, 0]
res['2023-triplet (10pock max)'] = binary_metrics(y, sc0.max(axis=1))
res['2023-triplet (10pock mean)'] = binary_metrics(y, sc0.mean(axis=1))
res['2023-triplet (cluster0)'] = binary_metrics(y, sc_c0)
res['2023-triplet (10pock LSE8)'] = binary_metrics(y, np.log(np.exp(sc0 * 8).mean(axis=1)) / 8)
# --- 2026 线：LOSO 留出 M4R 的 tuned 预测（同样未见 M4R 对）---
ens = None
for s in (20260925, 20260926, 20260927):
    p = pd.read_csv(rf'D:\CLC\project\results\drugclip_science2026\family_aug_v01\famaug_ep80_seed{s}.predictions.csv')
    q = p[p.target == 'M4R']
    q = q.assign(idx=[mi[v] for v in q.canonical_smiles.astype(str)])
    q = q.set_index('idx').loc[rows]
    res[f'2026-famaug-LOSO seed{s}'] = binary_metrics(y, q.tuned.to_numpy(float))
    res[f'2026-official-raw seed{s}'] = binary_metrics(y, q.official.to_numpy(float))
    col = q.tuned.to_numpy(float)
    ens = col if ens is None else ens + col
res['2026-famaug-LOSO 3seed-avg'] = binary_metrics(y, ens / 3)

order = ['roc_auc', 'ef1pct', 'bedroc_alpha20', 'pr_auc']
print()
print('%-28s %8s %8s %10s %8s' % ('model', 'ROC', 'EF1%', 'BEDROC20', 'PR'))
for k, v in res.items():
    print('%-28s %8.4f %8.2f %10.4f %8.4f' % (k, v['roc_auc'], v['ef1pct'], v['bedroc_alpha20'], v['pr_auc']))
json.dump(res, open(r'D:\CLC\project\results\drugclip_science2026\family_aug_v01\m4r_head2head_2023vs2026.json', 'w'), indent=2)
