# -*- coding: utf-8 -*-
"""Select a diverse, non-dominated PACER-M4 hypothesis portfolio.

The selector deliberately excludes Vina affinity and unvalidated potency
predictions from the Pareto objectives.  Output molecules are computational
hypotheses, not predicted-effective or experimentally validated PAMs.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs,RDConfig,RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
sys.path.append(str(Path(RDConfig.RDContribDir)/"SA_Score"));import sascorer
RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1];BASE=P/"results"/"pacer_candidates_v01";OUT=BASE/"final";OUT.mkdir(parents=True,exist_ok=True)

def fp(s):return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s),2,nBits=2048)
def pareto(a):
 # all columns are maximized
 n=len(a);keep=np.ones(n,bool)
 for i in range(n):
  if np.any(np.all(a>=a[i],axis=1)&np.any(a>a[i],axis=1)):keep[i]=False
 return keep
def clusters(fps,cutoff=.58):
 d=[]
 for i in range(1,len(fps)):d.extend(1-x for x in DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i]))
 return Butina.ClusterData(d,len(fps),cutoff,isDistData=True)
def choose_diverse(d,n):
 if d.empty:return d
 fps=[fp(x) for x in d.canonical_smiles];cs=clusters(fps);rank={idx:k for k,c in enumerate(cs) for idx in c};d=d.copy();d["diversity_cluster"]=[rank[i] for i in range(len(d))]
 # first pass: one member per cluster in deterministic Pareto/QED order
 picked=[]
 for _,x in d.sort_values(["pareto_front","pocket_residue_coverage","QED","SA_score"],ascending=[True,False,False,True]).iterrows():
  if x.diversity_cluster not in {z.diversity_cluster for z in picked}:picked.append(x)
  if len(picked)>=n:break
 if len(picked)<n:
  used={x.candidate_id for x in picked}
  for _,x in d.sort_values(["pareto_front","pocket_residue_coverage","QED"],ascending=[True,False,False]).iterrows():
   if x.candidate_id not in used:picked.append(x);used.add(x.candidate_id)
   if len(picked)>=n:break
 return pd.DataFrame(picked)
def fronts(d,cols):
 left=d.copy();left["pareto_front"]=-1;f=0
 while (left.pareto_front<0).any():
  ix=left.index[left.pareto_front<0];mask=pareto(left.loc[ix,cols].to_numpy(float));left.loc[ix[mask],"pareto_front"]=f;f+=1
 return left
def main():
 d=pd.read_csv(BASE/"predock_portfolio.csv");s=pd.read_csv(BASE/"7trs_docking"/"features.csv");s=s[s.error.fillna("")==""]
 # deterministic single-seed structural gate; future multi-seed MD/ensemble stage is separate.
 d=d.merge(s.sort_values("seed").drop_duplicates("candidate_id"),on="candidate_id",validate="one_to_one")
 d["SA_score"]=[float(sascorer.calculateScore(Chem.MolFromSmiles(x))) for x in d.canonical_smiles]
 d["novelty"]=1-d.max_tanimoto;d["inactive_safety_ref"]=1-d.strict_inactive_risk_ref
 # Independent objectives: state-pocket compatibility, medicinal chemistry,
 # bounded novelty, and a weak inactive-risk flag. Potency references excluded.
 cols=["pocket_residue_coverage","QED","novelty","inactive_safety_ref"]
 d=fronts(d,cols)
 allocations=[("local_exploitation",d[(d.domain_tier=="local")&(d.max_tanimoto>=.65)],12),
              ("local_diversification",d[(d.domain_tier=="local")&(d.max_tanimoto<.65)],6),
              ("exploratory_hypothesis",d[d.domain_tier=="exploratory"],6)]
 parts=[]
 for tier,x,n in allocations:
  q=choose_diverse(x,n);q["portfolio_role"]=tier;parts.append(q)
 final=pd.concat(parts,ignore_index=True).drop_duplicates("candidate_id");final.insert(0,"final_rank",np.arange(1,len(final)+1))
 final["claim_level"]="computational PAM hypothesis; requires functional ternary-complex assay"
 final["murcko_scaffold"]=[Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(x))) for x in final.canonical_smiles]
 # Add nearest known compound potency for transparent provenance.
 known=pd.read_csv(P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv").set_index("canonical_molecule_id")
 final["nearest_known_pEC50"]=[known.loc[x,"pEC50"] if x in known.index else np.nan for x in final.nearest_known_id]
 keep=["final_rank","candidate_id","canonical_smiles","murcko_scaffold","portfolio_role","generator","claim_level","pareto_front","diversity_cluster","nearest_known_id","nearest_known_pEC50","max_tanimoto","n_neighbors_045","local_knn_pEC50","local_neighbor_sd","pocket_residue_coverage","unique_pocket_residue_contacts","contacted_residues","vina_affinity","QED","SA_score","MW","logP","TPSA","PAINS","druglike","strict_inactive_risk_ref"]
 final[keep].to_csv(OUT/"final_candidate_hypotheses.csv",index=False);d.to_csv(OUT/"all_predock_with_structure.csv",index=False)
 audit={"n_docked":len(d),"n_final":len(final),"roles":final.portfolio_role.value_counts().to_dict(),"unique_ecfp_clusters":int(final.diversity_cluster.nunique()),"unique_murcko_scaffolds":int(final.murcko_scaffold.nunique()),"similarity_range":[float(final.max_tanimoto.min()),float(final.max_tanimoto.max())],"pocket_coverage_range":[float(final.pocket_residue_coverage.min()),float(final.pocket_residue_coverage.max())],"policy":"Pareto objectives exclude Vina affinity and potency predictions. Local kNN and inactive classifier are references only. Candidates are hypotheses requiring a functional PAM assay."}
 (OUT/"selection_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8");print(json.dumps(audit,indent=2))
if __name__=="__main__":main()
