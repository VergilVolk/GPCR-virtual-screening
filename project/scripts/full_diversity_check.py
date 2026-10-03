import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.RDLogger import DisableLog
DisableLog('rdApp.*')

df = pd.read_csv(r'D:\CLC\project\results\pacer_rerun_v01\module2_routed2605_theirAdapters.csv')
df['scaffold'] = [MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(str(s))) or str(s)
                  for s in df.smiles]
seen, reps = set(), []
for _, r in df.iterrows():
    if r.scaffold not in seen:
        seen.add(r.scaffold); reps.append(r)
    if len(reps) >= 12:
        break
print('全库骨架多样性前 12:')
for i, r in enumerate(reps, 1):
    tag = ''
    if str(r.smiles).startswith('Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)Nc3cnn(C4'):
        tag = ' <- PACER0010 同分子'
    print('  #%02d 全库排名%-5d %.4f  %s%s' % (i, r.final_rank, r.score, str(r.smiles)[:52], tag))
# MD 三候选各自的多样性序号
canon = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s)))
md3 = {'PACER0010': 'Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)Nc3cnn(C4Cc5ccccc5C4)c3)c(N)c12',
       'PACER0073': 'Cc1c(Cl)c2nnc(C)n2c2sc(N3CC(NC(=O)c4cc(F)nc(F)c4)C3)c(N)c12',
       'PACER0027': 'Cc1c(Cl)c2nncn2c2sc(N3CC(N4CC(NC5Cc6ccccc6C5)C4)C3)c(N)c12'}
cs = {canon(v): k for k, v in md3.items()}
df['canon'] = [canon(s) for s in df.smiles]
seen, div_rank = set(), {}
for _, r in df.iterrows():
    if r.scaffold not in seen:
        seen.add(r.scaffold)
        if r.canon in cs:
            div_rank[cs[r.canon]] = len(seen)
print()
for pid, d in sorted(div_rank.items(), key=lambda x: x[1]):
    print('%s -> 骨架多样性序号 #%d' % (pid, d))
