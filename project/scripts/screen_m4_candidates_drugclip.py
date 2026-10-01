#!/usr/bin/env python3
"""Fit frozen M4 pocket-ranking adapters and prioritize generated candidates."""
from __future__ import annotations
import argparse,json,random
from pathlib import Path
import numpy as np,pandas as pd,torch
from rdkit import Chem
from scipy.stats import rankdata
from torch.nn import functional as F
from train_drugclip_m4_pocket_realign import PocketRealign

def can(s):
 m=Chem.MolFromSmiles(str(s));
 if m is None: raise ValueError(s)
 return Chem.MolToSmiles(m,canonical=True)

def load(path,smiles):
 a=np.load(path,allow_pickle=False); ids=list(map(str,a['molecule_ids'])); ix={v:i for i,v in enumerate(ids)}; cs=smiles.map(can)
 return a,np.asarray([ix[v] for v in cs])

def fit(train_mol,pockets,y,seed,args):
 random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
 x=torch.tensor(train_mol,dtype=torch.float32);p=torch.tensor(pockets,dtype=torch.float32);yt=torch.tensor(y,dtype=torch.float32)
 model=PocketRealign(x.shape[1],args.rank,1);opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
 pw=torch.tensor(float((y==0).sum()/(y==1).sum()));ox,op=F.normalize(x,dim=1),F.normalize(p,dim=1)
 for _ in range(args.epochs):
  opt.zero_grad();zm,zp,_,retr,logit,_=model(x,p);loss=F.binary_cross_entropy_with_logits(logit,yt,pos_weight=pw)
  positive=retr[yt.bool()];negative=retr[~yt.bool()];hard=torch.topk(negative,min(args.hard_negative_top_k,len(negative))).values
  if len(positive)>args.max_positive: positive=positive[torch.randperm(len(positive))[:args.max_positive]]
  loss=loss+args.rank_weight*F.relu(args.margin-positive[:,None]+hard[None,:]).mean()
  loss=loss+args.preserve_weight*((1-(zm*ox).sum(1)).mean()+(1-(zp*op).sum(1)).mean());loss.backward();opt.step()
 return model

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--train-embeddings',type=Path,required=True);ap.add_argument('--benchmark',type=Path,required=True)
 ap.add_argument('--candidate-embeddings',type=Path,required=True);ap.add_argument('--candidates',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--seeds',default='20260924,20260925,20260926');ap.add_argument('--epochs',type=int,default=300);ap.add_argument('--rank',type=int,default=16)
 ap.add_argument('--lr',type=float,default=1e-3);ap.add_argument('--weight-decay',type=float,default=1e-3);ap.add_argument('--margin',type=float,default=.05)
 ap.add_argument('--rank-weight',type=float,default=.5);ap.add_argument('--preserve-weight',type=float,default=.2);ap.add_argument('--hard-negative-top-k',type=int,default=16);ap.add_argument('--max-positive',type=int,default=128)
 args=ap.parse_args();torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
 bench=pd.read_csv(args.benchmark);cand=pd.read_csv(args.candidates);ta,ti=load(args.train_embeddings,bench.canonical_smiles);ca,ci=load(args.candidate_embeddings,cand.canonical_smiles)
 tm=ta['molecule_embeddings'].astype(np.float32)[ti];cm=ca['molecule_embeddings'].astype(np.float32)[ci];p=ta['pocket_embeddings'].astype(np.float32);y=bench.target.to_numpy(int)
 probs=[];states=[]
 for seed in map(int,args.seeds.split(',')):
  model=fit(tm,p,y,seed,args);model.eval()
  with torch.no_grad(): _,_,sim,_,logit,_=model(torch.tensor(cm),torch.tensor(p))
  probs.append(torch.sigmoid(logit).numpy());states.append(sim.numpy())
 prob=np.mean(probs,axis=0);state=np.mean(states,axis=0);unc=np.std(probs,axis=0)
 out=cand.copy();out['drugclip_m4_probability']=prob;out['drugclip_seed_sd']=unc
 for j,name in enumerate(map(str,ta['pocket_ids'])):out[f'drugclip_{name}']=state[:,j]
 out['drugclip_rank']=rankdata(-prob,method='min').astype(int)
 pam2d=1-out.strict_inactive_risk_ref.to_numpy(float);out['pam2d_rank']=rankdata(-pam2d,method='min').astype(int)
 out['fixed_rank_fusion']=.5*rankdata(prob)/len(out)+.5*rankdata(pam2d)/len(out);out['fusion_rank']=rankdata(-out.fixed_rank_fusion,method='min').astype(int)
 out=out.sort_values(['fusion_rank','drugclip_seed_sd','candidate_id']);args.output.parent.mkdir(parents=True,exist_ok=True);out.to_csv(args.output,index=False)
 audit={'evidence_level':'prospective_computational_candidate_prioritization','n_candidates':len(out),'seeds':args.seeds,'fusion':'fixed 50:50 percentile rank: DrugCLIP M4 pocket-ranking + (1-inactive-risk)',
 'top10':out[['candidate_id','fusion_rank','drugclip_rank','pam2d_rank','drugclip_m4_probability','drugclip_seed_sd']].head(10).to_dict('records'),
 'claim_boundary':'Ranking hypotheses only; no candidate is experimentally confirmed as a PAM.'}
 args.output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
