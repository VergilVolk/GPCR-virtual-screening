#!/usr/bin/env python3
"""Fine-tune DrugCLIP's internal shared space on real M1--M5 activity data."""
from __future__ import annotations
import argparse, copy, json, random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F

class Proj(nn.Module):
    def __init__(self,s):
        super().__init__(); self.linear1=nn.Linear(s['linear1.weight'].shape[1],s['linear1.weight'].shape[0]); self.linear2=nn.Linear(s['linear2.weight'].shape[1],s['linear2.weight'].shape[0]); self.load_state_dict(s)
    def forward(self,x): return F.normalize(self.linear2(F.relu(self.linear1(x))),dim=-1)
class Model(nn.Module):
    def __init__(self,s):
        super().__init__(); self.mol_project=Proj(s['mol_project']); self.pocket_project=Proj(s['pocket_project']); self.log_scale=nn.Parameter(s['logit_scale'].float().reshape(()).clone()); self.bias=nn.Parameter(torch.tensor(0.))
    def forward(self,m,p):
        zm,zp=self.mol_project(m),self.pocket_project(p); sim=(zm*zp).sum(1); return zm,zp,sim,self.log_scale.exp().clamp(max=100)*sim+self.bias

def summarize(y,s,sub):
    per={}
    for name in sorted(set(sub)):
        q=sub==name
        per[name]={'n':int(q.sum()),'roc_auc':float(roc_auc_score(y[q],s[q])),'average_precision':float(average_precision_score(y[q],s[q]))}
    return {'micro_roc_auc':float(roc_auc_score(y,s)),'micro_average_precision':float(average_precision_score(y,s)),'macro_subtype_auc':float(np.mean([v['roc_auc'] for v in per.values()])),'per_subtype':per}

def train(mrep,prep,pairs,initial,mode,seed,args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); rng=np.random.default_rng(seed)
    mrep=torch.as_tensor(mrep,dtype=torch.float32); prep=torch.as_tensor(prep,dtype=torch.float32)
    model=Model(copy.deepcopy(initial)); frozen=Model(copy.deepcopy(initial)).eval()
    for q in frozen.parameters(): q.requires_grad_(False)
    with torch.no_grad(): z0m=frozen.mol_project(mrep); z0p=frozen.pocket_project(prep)
    pools={}
    for subtype,g in pairs.groupby('subtype'):
        pos=g.index[g.target.eq(1)].to_numpy(); neg=g.index[g.target.eq(0)].to_numpy()
        if mode=='hard':
            with torch.no_grad():
                mi=torch.as_tensor(pairs.loc[neg,'mol_index'].to_numpy(),dtype=torch.long); pi=torch.as_tensor(pairs.loc[neg,'pocket_index'].to_numpy(),dtype=torch.long)
                score=(z0m[mi]*z0p[pi]).sum(1).numpy(); neg=neg[np.argsort(-score)[:min(args.hard_pool,len(neg))]]
        pools[subtype]=(pos,neg)
    opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        for _ in range(args.steps_per_epoch):
            chosen=[]
            for pos,neg in pools.values():
                n=args.per_subtype
                chosen.extend(rng.choice(pos,n,replace=len(pos)<n)); chosen.extend(rng.choice(neg,n,replace=len(neg)<n))
            batch=pairs.loc[np.asarray(chosen,int)]; mi=torch.as_tensor(batch.mol_index.to_numpy(),dtype=torch.long); pi=torch.as_tensor(batch.pocket_index.to_numpy(),dtype=torch.long); y=torch.as_tensor(batch.target.to_numpy(),dtype=torch.float32)
            opt.zero_grad(); zm,zp,sim,logit=model(mrep[mi],prep[pi]); loss=F.binary_cross_entropy_with_logits(logit,y)
            if mode!='bce':
                rank=[]
                offset=0
                for _subtype in pools:
                    block=sim[offset:offset+2*args.per_subtype]; rank.append(F.relu(args.margin-block[:args.per_subtype,None]+block[args.per_subtype:][None,:]).mean()); offset+=2*args.per_subtype
                loss=loss+args.rank_weight*torch.stack(rank).mean()
            preserve=(1-(zm*z0m[mi]).sum(1)).mean()+(1-(zp*z0p[pi]).sum(1)).mean(); loss=loss+args.preserve_weight*preserve
            loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),1.); opt.step()
    return model

def predict(model,mrep,prep,pairs):
    mi=torch.as_tensor(pairs.mol_index.to_numpy(),dtype=torch.long); pi=torch.as_tensor(pairs.pocket_index.to_numpy(),dtype=torch.long)
    with torch.no_grad(): return model(torch.as_tensor(mrep)[mi],torch.as_tensor(prep)[pi])[2].numpy()

def main():
    p=argparse.ArgumentParser(); p.add_argument('--representations',type=Path,required=True); p.add_argument('--projection',type=Path,required=True); p.add_argument('--pairs',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--seed',type=int,default=20260924); p.add_argument('--epochs',type=int,default=12); p.add_argument('--steps-per-epoch',type=int,default=15); p.add_argument('--per-subtype',type=int,default=24); p.add_argument('--lr',type=float,default=1e-4); p.add_argument('--weight-decay',type=float,default=1e-3); p.add_argument('--margin',type=float,default=.05); p.add_argument('--rank-weight',type=float,default=.5); p.add_argument('--preserve-weight',type=float,default=.5); p.add_argument('--hard-pool',type=int,default=96); a=p.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(a.representations,allow_pickle=False); initial=torch.load(a.projection,map_location='cpu'); pairs=pd.read_csv(a.pairs)
    mol_lookup={str(v):i for i,v in enumerate(arc['molecule_ids'])}; pocket_lookup={str(v).split('_')[0]:i for i,v in enumerate(arc['pocket_ids'])}
    keep=pairs.canonical_smiles.astype(str).isin(mol_lookup); omitted=int((~keep).sum()); pairs=pairs[keep].reset_index(drop=True); pairs['mol_index']=[mol_lookup[str(v)] for v in pairs.canonical_smiles]; pairs['pocket_index']=[pocket_lookup[str(v)] for v in pairs.subtype]
    mrep=arc['molecule_representations'].astype(np.float32); prep=arc['pocket_representations'].astype(np.float32); folds=pairs.scaffold_fold.to_numpy(int); y=pairs.target.to_numpy(int); sub=pairs.subtype.to_numpy(str); modes=['frozen','bce','random','hard']; pred={m:np.full(len(pairs),np.nan) for m in modes}; frozen=Model(copy.deepcopy(initial)).eval()
    for fold in sorted(np.unique(folds)):
        tr=pairs[folds!=fold].copy(); te=pairs[folds==fold].copy(); pred['frozen'][folds==fold]=predict(frozen,mrep,prep,te)
        for mode in modes[1:]:
            model=train(mrep,prep,tr,initial,mode,a.seed+101*int(fold),a); pred[mode][folds==fold]=predict(model,mrep,prep,te)
    report={'method':'DrugCLIP internal M1-M5 activity-domain fine-tuning','split':'5-fold molecule Murcko scaffold OOF','n_pairs':len(pairs),'omitted_pairs':omitted,'metrics':{m:summarize(y,pred[m],sub) for m in modes},'claim_boundary':'ChEMBL muscarinic activity, not allosteric mechanism or M4 PAM efficacy.'}
    final=train(mrep,prep,pairs,initial,'hard',a.seed,a); checkpoint={'mol_project':final.mol_project.state_dict(),'pocket_project':final.pocket_project.state_dict(),'logit_scale':final.log_scale.detach().cpu().reshape(1),'classification_bias':final.bias.detach().cpu(),'source_checkpoint':initial.get('source_checkpoint'),'metadata':{'stage':'M1-M5 activity domain','n_pairs':len(pairs)}}
    a.output.parent.mkdir(parents=True,exist_ok=True); ck=a.output.with_suffix('.projection.pt'); torch.save(checkpoint,ck); report['checkpoint']=str(ck); a.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); pd.DataFrame({'canonical_molecule_id':pairs.canonical_molecule_id,'subtype':sub,'target':y,'fold':folds,**pred}).to_csv(a.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
