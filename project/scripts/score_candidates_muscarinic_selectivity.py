#!/usr/bin/env python3
"""Score candidate M1-M5 subtype preference with the frozen ECFP3 adapter."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from select_molecular_adapter_on_non_m4 import Adapter,feature


def main():
    p=argparse.ArgumentParser();p.add_argument('--candidates',type=Path,required=True);p.add_argument('--pocket-representations',type=Path,required=True);p.add_argument('--selection-report',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    table=pd.read_csv(a.candidates);report=json.loads(a.selection_report.read_text(encoding='utf-8'));arc=np.load(a.pocket_representations,allow_pickle=False);pocket=F.normalize(torch.as_tensor(arc['pocket_embeddings'].astype(np.float32)),dim=-1);names=[str(v).split('_')[0] for v in arc['pocket_ids']];pi={v:i for i,v in enumerate(names)};x=feature(table.canonical_smiles.astype(str).tolist(),report['selected']);values=[]
    for ckpath in report['checkpoints']:
        ck=torch.load(ckpath,map_location='cpu');model=Adapter(x.shape[1],int(ck['seed']));model.load_state_dict(ck['state_dict']);model.eval()
        with torch.inference_mode():values.append((model(x)@pocket.T).numpy())
    scores=np.mean(values,axis=0);m4=scores[:,pi['M4']];others=np.delete(scores,pi['M4'],axis=1);margin=m4-others.max(axis=1);rank=1+(scores>m4[:,None]).sum(axis=1);seed_margin=np.asarray([v[:,pi['M4']]-np.delete(v,pi['M4'],axis=1).max(axis=1) for v in values]);out=table.copy()
    for j,name in enumerate(names):out[f'ecfp3_{name}_score']=scores[:,j]
    out['ecfp3_m4_vs_best_other_margin']=margin;out['ecfp3_m4_rank_of_5']=rank;out['ecfp3_margin_seed_sd']=seed_margin.std(axis=0);out['ecfp3_selectivity_rank']=pd.Series(-margin).rank(method='min').astype(int);out=out.sort_values(['ecfp3_selectivity_rank','fusion_rank']).reset_index(drop=True)
    top=out.head(20)[['candidate_id','ecfp3_selectivity_rank','ecfp3_m4_rank_of_5','ecfp3_m4_vs_best_other_margin','ecfp3_margin_seed_sd','fusion_rank']].to_dict('records');audit={'method':'Frozen three-seed ECFP3 molecule adapter against frozen DrugCLIP M1-M5 pockets','n_candidates':len(out),'m4_top1_count':int((out.ecfp3_m4_rank_of_5==1).sum()),'m4_top2_count':int((out.ecfp3_m4_rank_of_5<=2).sum()),'top20':top,'claim_boundary':'Subtype-selectivity hypotheses only; scores do not establish binding, allostery, PAM function, or efficacy.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);out.to_csv(a.output,index=False);a.output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
