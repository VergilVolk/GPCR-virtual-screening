import numpy as np, torch, pandas as pd, json, sys
sys.path.insert(0, r'D:\CLC\project\scripts')
from finetune_drugclip_muscarinic_triplet import Proj

c = r'D:\CLC\project\data\drugclip_muscarinic_family_aug_v01\cand200'
lib = np.load(r'D:\CLC\project\results\gpcr_drugclip_screening_v01\science2026_ensemble_embeddings.npz', allow_pickle=False)
gp_ids = [str(x) for x in lib['pocket_ids']]
m4_idx = [i for i, p in enumerate(gp_ids) if p.startswith('M4R_')]
pocket_rep = torch.as_tensor(lib['pocket_representations'][m4_idx].astype(np.float32))

cand = np.load(c + r'\cands_science2026.npz', allow_pickle=False)
ids = [str(x) for x in cand['molecule_ids']]
mol_rep = torch.as_tensor(cand['molecule_representations'].astype(np.float32))
old = pd.read_csv(r'D:\CLC\project\results\drugclip_m4_candidate_screen_v01\ranked_candidates.csv')
smi2id = old.set_index('canonical_smiles')

rows = []
scores_by_seed = {}
for s in (20260925, 20260926, 20260927):
    d = torch.load(rf'D:\CLC\project\results\drugclip_science2026\family_aug_v01\science2026_13target_ep80_famaug.seed{s}.projection.pt', map_location='cpu', weights_only=False)
    pm, pp = Proj(d['mol_project']).eval(), Proj(d['pocket_project']).eval()
    with torch.inference_mode():
        sc = (pm(mol_rep) @ pp(pocket_rep).T).max(axis=1).values.numpy()
    scores_by_seed[s] = sc
ens = np.mean(list(scores_by_seed.values()), axis=0)
order = np.argsort(-ens)
out = pd.DataFrame({'candidate_smiles': ids,
                    'famaug_m4_score': ens,
                    'seed25': scores_by_seed[20260925], 'seed26': scores_by_seed[20260926], 'seed27': scores_by_seed[20260927]})
out = out.iloc[order].reset_index(drop=True)
out['new_rank'] = out.index + 1
out['old_fusion_rank'] = out.candidate_smiles.map(smi2id['fixed_rank_fusion'] if 'fixed_rank_fusion' in smi2id else smi2id.iloc[:, 0])
out['candidate_id'] = out.candidate_smiles.map(old.set_index('canonical_smiles')['candidate_id'])
out['drugclip_m4_probability_old'] = out.candidate_smiles.map(old.set_index('canonical_smiles')['drugclip_m4_probability'])
out.to_csv(c + r'\cands200_famaug_scores.csv', index=False)
rho = out[['new_rank', 'old_fusion_rank']].dropna()
from scipy.stats import spearmanr
print('n matched:', len(rho))
print('spearman(new famaug rank, old fusion rank) = %.3f' % spearmanr(rho.new_rank, rho.old_fusion_rank).statistic)
top50_old = rho[rho.old_fusion_rank <= 50]
print('old top50 -> new rank median %.0f, in new top50: %d, top100: %d' % (top50_old.new_rank.median(), (top50_old.new_rank <= 50).sum(), (top50_old.new_rank <= 100).sum()))
print('new top10 candidate_ids:', out.head(10).candidate_id.tolist())
