import pandas as pd
from scipy.stats import spearmanr
out = pd.read_csv(r'D:\CLC\project\data\drugclip_muscarinic_family_aug_v01\cand200\cands200_famaug_scores.csv')
old = pd.read_csv(r'D:\CLC\project\results\drugclip_m4_candidate_screen_v01\ranked_candidates.csv')
out['old_drugclip_rank'] = out.candidate_smiles.map(old.set_index('canonical_smiles')['drugclip_rank'])
rho = out.dropna(subset=['old_drugclip_rank'])
print('n =', len(rho))
print('spearman(new famaug rank, old drugclip rank) = %+.3f' % spearmanr(rho.new_rank, rho.old_drugclip_rank).statistic)
for k in (20, 50, 100):
    old_top = set(rho.nsmallest(k, 'old_drugclip_rank').candidate_id)
    new_top = set(rho.nsmallest(k, 'new_rank').candidate_id)
    print('top%d overlap: %d/%d (jaccard %.2f)' % (k, len(old_top & new_top), k, len(old_top & new_top) / len(old_top | new_top)))
old_top50 = rho.nsmallest(50, 'old_drugclip_rank')
print('old drugclip-top50 -> new rank median %.0f (of 200)' % old_top50.new_rank.median())
both = old_top50.nsmallest(25, 'old_drugclip_rank').nsmallest(25, 'new_rank')
print('dual-top candidates (old<=50 & new<=25):', both.candidate_id.tolist())
