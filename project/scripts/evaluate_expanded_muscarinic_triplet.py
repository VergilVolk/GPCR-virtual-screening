#!/usr/bin/env python3
"""Train expanded muscarinic triplets and evaluate on molecule-disjoint strict pairs."""
from __future__ import annotations
import argparse, copy, json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import torch
from finetune_drugclip_muscarinic_triplet import Model, fit, clustered_bootstrap

def score_pairs(model,m,p,d,mi,pi):
    s=model.matrix(m,p).detach().numpy()
    return np.asarray([s[mi[r.canonical_smiles],pi[r.positive_subtype]]-s[mi[r.canonical_smiles],pi[r.negative_subtype]] for r in d.itertuples()])

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--representations',type=Path,required=True); ap.add_argument('--projection',type=Path,required=True); ap.add_argument('--train',type=Path,required=True); ap.add_argument('--test',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--seeds',default='20260924,20260925,20260926'); ap.add_argument('--epochs',type=int,default=500); ap.add_argument('--lr',type=float,default=1e-3); ap.add_argument('--margin',type=float,default=.1); ap.add_argument('--preserve-weight',type=float,default=.2); args=ap.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(args.representations,allow_pickle=False); initial=torch.load(args.projection,map_location='cpu'); mids=list(map(str,arc['molecule_ids'])); pids=list(map(str,arc['pocket_ids'])); mi={v:i for i,v in enumerate(mids)}; pi={v.split('_')[0]:i for i,v in enumerate(pids)}; train=pd.read_csv(args.train); test=pd.read_csv(args.test)
    train_missing=sorted(set(train.canonical_smiles)-set(mi)); test_missing=sorted(set(test.canonical_smiles)-set(mi)); train=train[train.canonical_smiles.isin(mi)].reset_index(drop=True); test=test[test.canonical_smiles.isin(mi)].reset_index(drop=True); overlap=set(train.canonical_smiles)&set(test.canonical_smiles)
    m=torch.tensor(arc['molecule_representations'],dtype=torch.float32); p=torch.tensor(arc['pocket_representations'],dtype=torch.float32); frozen=Model(copy.deepcopy(initial)).eval(); base=score_pairs(frozen,m,p,test,mi,pi)
    cfg=SimpleNamespace(lr=args.lr,weight_decay=1e-3,epochs=args.epochs,margin=args.margin,preserve_weight=args.preserve_weight,train_scope='pocket'); rows=[]; targeted=[]; controls=[]; trained=[]
    for seed in map(int,args.seeds.split(',')):
        target_model=fit(m,p,train,mi,pi,initial,seed,cfg,False); random_model=fit(m,p,train,mi,pi,initial,seed,cfg,True); t=score_pairs(target_model,m,p,test,mi,pi); r=score_pairs(random_model,m,p,test,mi,pi); targeted.append(t); controls.append(r)
        trained.append((seed,target_model))
        for i,row in test.iterrows(): rows.append({'canonical_smiles':row.canonical_smiles,'positive_subtype':row.positive_subtype,'negative_subtype':row.negative_subtype,'seed':seed,'frozen_delta':float(base[i]),'random_control_delta':float(r[i]),'targeted_delta':float(t[i])})
    targeted_mean=np.mean(targeted,axis=0); control_mean=np.mean(controls,axis=0)
    report={'method':'Expanded same-molecule muscarinic triplet with molecule-disjoint strict external test','train_pairs':len(train),'train_molecules':train.canonical_smiles.nunique(),'test_pairs':len(test),'test_molecules':test.canonical_smiles.nunique(),'train_test_molecule_overlap':len(overlap),'train_missing_molecules':len(train_missing),'test_missing_molecules':len(test_missing),'frozen_accuracy':float((base>0).mean()),'random_control_accuracy':float((control_mean>0).mean()),'targeted_accuracy':float((targeted_mean>0).mean()),'targeted_minus_frozen_molecule_cluster_bootstrap_95ci':clustered_bootstrap(test,base,targeted_mean),'targeted_minus_random_molecule_cluster_bootstrap_95ci':clustered_bootstrap(test,control_mean,targeted_mean),'per_seed_targeted_accuracy':[float((x>0).mean()) for x in targeted],'per_seed_random_accuracy':[float((x>0).mean()) for x in controls],'claim_boundary':'Molecule-disjoint M1/M2/M3/M5 activity ranking; not M4 PAM efficacy.'}
    args.output.parent.mkdir(parents=True,exist_ok=True); checkpoints=[]
    for seed,model in trained:
        ck={'mol_project':model.mol_project.state_dict(),'pocket_project':model.pocket_project.state_dict(),'logit_scale':initial['logit_scale'],'source_checkpoint':initial.get('source_checkpoint'),'metadata':{'stage':'expanded muscarinic triplet','train_pairs':len(train),'strict_test_excluded':True,'seed':seed}}
        ckpath=args.output.with_name(f'{args.output.stem}.seed{seed}.projection.pt'); torch.save(ck,ckpath); checkpoints.append(str(ckpath))
    # Retain the original path for backward compatibility; it is the first
    # fixed seed, while deployment uses the unweighted three-seed ensemble.
    torch.save(torch.load(checkpoints[0],map_location='cpu'),args.output.with_suffix('.projection.pt'))
    report['checkpoints']=checkpoints; report['deployment_rule']='Unweighted average of all fixed-seed scores.'; args.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); pd.DataFrame(rows).to_csv(args.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
