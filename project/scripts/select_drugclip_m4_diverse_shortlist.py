#!/usr/bin/env python3
"""Select a scaffold-diverse, multi-evidence shortlist from frozen candidate ranks."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def scaffold(s):
 m=Chem.MolFromSmiles(str(s));return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else ""

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--ranked',type=Path,required=True);ap.add_argument('--final',type=Path,required=True)
 ap.add_argument('--round2',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--n',type=int,default=12);args=ap.parse_args()
 d=pd.read_csv(args.ranked);final=pd.read_csv(args.final);round2=pd.read_csv(args.round2)
 d['murcko_scaffold_recomputed']=d.canonical_smiles.map(scaffold);d=d.merge(final[['candidate_id','final_rank','pareto_front','diversity_cluster','vina_affinity','pocket_residue_coverage']],on='candidate_id',how='left',suffixes=('','_old'))
 rcols=['candidate_id','dynamic_pareto','BE_min','BE_avg','vina_sd','pocket_coverage_mean','pocket_coverage_min']
 d=d.merge(round2[rcols],on='candidate_id',how='left');d['prior_final24']=d.final_rank.notna();d['prior_dynamic16']=d.dynamic_pareto.notna()
 selected=[];seen=set()
 # First pass: one per scaffold, ordered by frozen fusion. Second pass only if needed.
 for row in d.sort_values(['fusion_rank','drugclip_seed_sd','candidate_id']).itertuples():
  if row.murcko_scaffold_recomputed in seen:continue
  selected.append(row.Index);seen.add(row.murcko_scaffold_recomputed)
  if len(selected)>=args.n:break
 if len(selected)<args.n:
  selected += [i for i in d.index if i not in selected][:args.n-len(selected)]
 out=d.loc[selected].copy().sort_values('fusion_rank');out.insert(0,'shortlist_rank',range(1,len(out)+1))
 out['claim']='multi-evidence computational hypothesis; requires ACh-conditioned functional PAM assay'
 args.output.parent.mkdir(parents=True,exist_ok=True);out.to_csv(args.output,index=False)
 audit={'evidence_level':'prospective_computational_diverse_shortlist','n_selected':len(out),'n_unique_scaffolds':out.murcko_scaffold_recomputed.nunique(),
  'n_prior_final24':int(out.prior_final24.sum()),'n_prior_dynamic16':int(out.prior_dynamic16.sum()),
  'selected':out[['shortlist_rank','candidate_id','fusion_rank','drugclip_rank','pam2d_rank','prior_final24','prior_dynamic16','drugclip_seed_sd']].to_dict('records'),
  'claim_boundary':'Computational prioritization only; candidates are not confirmed PAMs.'}
 args.output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
