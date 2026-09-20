# -*- coding: utf-8 -*-
"""Build a deduplicated scaffold-held-out M4R early-enrichment benchmark."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem,RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1];SRC=P/"data"/"miao2026_m4r";OUT=P/"data"/"benchmarks"/"m4r_enrichment_v1";OUT.mkdir(parents=True,exist_ok=True)
def canon(s):
 m=Chem.MolFromSmiles(s);return Chem.MolToSmiles(m,isomericSmiles=True) if m else None
def scaffold(s):
 m=Chem.MolFromSmiles(s);sc=MurckoScaffold.GetScaffoldForMol(m);x=Chem.MolToSmiles(sc,isomericSmiles=True)
 return x if x else "ACYCLIC_"+hashlib.sha1(s.encode()).hexdigest()[:16]
def assign(groups,k=5):
 state=[{"n":0,"p":0} for _ in range(k)];ans={}
 positive=groups[groups.pos>0].sort_values(["pos","n","scaffold"],ascending=[False,False,True])
 decoy=groups[groups.pos==0].sort_values(["n","scaffold"],ascending=[False,True])
 # Longest-processing-time bin packing: balance active counts first.
 for r in positive.itertuples():
  f=min(range(k),key=lambda j:(state[j]["p"],state[j]["n"],j));ans[r.scaffold]=f;state[f]["n"]+=r.n;state[f]["p"]+=r.pos
 # Pure-decoy scaffolds then balance the total screening-set size.
 for r in decoy.itertuples():
  f=min(range(k),key=lambda j:(state[j]["n"],j));ans[r.scaffold]=f;state[f]["n"]+=r.n
 return ans,state
def main():
 pdb=pd.read_csv(SRC/"M4R_PDB_Vina_scores.csv");fw=pd.read_csv(SRC/"M4R_Ensemble_Vina_framewise_scores.csv")
 pdb["canonical_smiles"]=[canon(x) for x in pdb.SMILES];invalid=pdb.canonical_smiles.isna().sum();pdb=pdb.dropna(subset=["canonical_smiles"])
 conflicts=pdb.groupby("canonical_smiles").is_active.nunique();bad=set(conflicts[conflicts>1].index);pdb=pdb[~pdb.canonical_smiles.isin(bad)].copy()
 # Duplicate structures are one molecule. Keep median PDB score and deterministic ID.
 base=pdb.groupby("canonical_smiles",as_index=False).agg(ligand_id=("ligand_id","min"),is_active=("is_active","first"),pdb_vina=("vina_score","median"),n_source_rows=("ligand_id","size"))
 bes=[f"cluster{i:02d}_BE" for i in range(10)];fw["canonical_smiles"]=[canon(x) for x in fw.SMILES];ens=fw.groupby("canonical_smiles")[bes].median();base=base.merge(ens,left_on="canonical_smiles",right_index=True,how="left")
 base["ensemble_BEmin"]=base[bes].min(axis=1,skipna=True);base["ensemble_BEavg"]=base[bes].mean(axis=1,skipna=True);base["ensemble_BEsd"]=base[bes].std(axis=1,skipna=True);base["n_states_available"]=base[bes].notna().sum(axis=1)
 base["murcko_scaffold"]=[scaffold(x) for x in base.canonical_smiles];g=base.groupby("murcko_scaffold").is_active.agg(n="size",pos="sum").reset_index().rename(columns={"murcko_scaffold":"scaffold"});mp,state=assign(g);base["scaffold_fold"]=[mp[x] for x in base.murcko_scaffold]
 base.insert(0,"molecule_id",[f"M4R{i:06d}" for i in range(1,len(base)+1)]);base.to_csv(OUT/"molecules.csv",index=False)
 fold=base.groupby("scaffold_fold").is_active.agg(n="size",actives="sum");fold["decoys"]=fold.n-fold.actives
 audit={"source_pdb_rows":len(pd.read_csv(SRC/"M4R_PDB_Vina_scores.csv")),"invalid_smiles":int(invalid),"conflicting_canonical_smiles_excluded":len(bad),"deduplicated_molecules":len(base),"actives":int(base.is_active.sum()),"decoys":int((base.is_active==0).sum()),"unique_scaffolds":int(base.murcko_scaffold.nunique()),"folds":fold.reset_index().to_dict("records"),"rules":"Canonical duplicates collapsed; conflicting labels excluded; entire Murcko scaffolds assigned to one deterministic fold for both actives and decoys."}
 (OUT/"audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8");print(json.dumps(audit,indent=2))
if __name__=="__main__":main()
