import numpy as np, torch, csv, sys
sys.path.insert(0, r'D:\CLC\project\scripts')
from finetune_drugclip_muscarinic_triplet import Proj

z = np.load(r'D:\CLC\project\results\gpcr_drugclip_screening_v01\ensemble_embeddings.npz', allow_pickle=False)
ids = [str(x) for x in z['molecule_ids']]
mol = torch.as_tensor(z['molecule_representations'].astype(np.float32))
gp = [str(x) for x in z['pocket_ids']]
m4 = [i for i, p in enumerate(gp) if p.startswith('M4R')]
pocket = torch.as_tensor(z['pocket_representations'][m4].astype(np.float32))

per_seed, mean_seed = [], []
for s in (20260924, 20260925, 20260926):
    d = torch.load(rf'D:\CLC\project\results\muscarinic_triplet_expansion_v01\external_result.seed{s}.projection.pt',
                   map_location='cpu', weights_only=False)
    pm, pp = Proj(d['mol_project']).eval(), Proj(d['pocket_project']).eval()
    with torch.inference_mode():
        sc = (pm(mol) @ pp(pocket).T).numpy()          # n_mol x 10
    per_seed.append(sc.max(axis=1))
    mean_seed.append(sc.mean(axis=1))                   # head2head 胜出配置（mean 池化）

ens = np.mean(mean_seed, axis=0); ens_max = np.mean(per_seed, axis=0)
order = np.argsort(-ens)
out = r'D:\CLC\project\results\drugclip_science2026\family_aug_v01\m4_library_oldtriplet_ranking.csv'
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'molecule', 'score_mean_ens', 'score_max_ens',
                'seed24_mean', 'seed25_mean', 'seed26_mean'])
    for r, i in enumerate(order, 1):
        w.writerow([r, ids[i], round(float(ens[i]), 6), round(float(ens_max[i]), 6),
                    round(float(mean_seed[0][i]), 6), round(float(mean_seed[1][i]), 6), round(float(mean_seed[2][i]), 6)])
sd = np.std([s[order[:100]] for s in per_seed], axis=0).mean()
print('saved:', out)
print('top-1% (285) 与 famaug 排名的重叠:')
new = open(r'D:\CLC\project\results\drugclip_science2026\family_aug_v01\m4_library_famaug_ranking.csv', encoding='utf-8').read().splitlines()[1:286]
new_set = {line.split(',')[1] for line in new}
old_top = {ids[i] for i in order[:285]}
print('  overlap = %d / 285' % len(new_set & old_top))
print('top10:')
for r, i in enumerate(order[:10], 1):
    print(' %2d %s %.4f' % (r, ids[i][:60], ens[i]))
