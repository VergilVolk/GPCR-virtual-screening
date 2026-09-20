# -*- coding: utf-8 -*-
"""PACER-FS v0.1: strict episodic few-shot series adaptation baselines.

The outer medicinal-chemistry series is absent from all model fitting. Only
0/1/3/5 labelled anchors from that series are revealed at adaptation time.
"""
from __future__ import annotations
import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1];DATA=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv";OUT=P/"results"/"pacer_fs_v01";OUT.mkdir(parents=True,exist_ok=True)
SHOTS=(0,1,3,5);REPEATS=20;MIN_GROUP=12
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]
def features(smiles):
 x=[];ds=[];bv=[]
 for s in smiles:
  m=Chem.MolFromSmiles(s);f=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048);a=np.zeros(2048,np.float32);DataStructs.ConvertToNumpyArray(f,a);x.append(a);ds.append([q(m) for q in DESC]);bv.append(f)
 return np.asarray(x,np.float32),np.asarray(ds,np.float32),bv
def sims(a,b,bv):return np.asarray([[DataStructs.TanimotoSimilarity(bv[i],bv[j]) for j in b] for i in a],float)
def base_model(X,D,y,tr,te):
 sel=np.argsort(X[tr].var(0))[-1024:];xx=np.c_[X[:,sel],D];m=LGBMRegressor(n_estimators=350,learning_rate=.03,num_leaves=15,max_depth=5,min_child_samples=12,reg_alpha=.1,reg_lambda=1,random_state=42,verbosity=-1,n_jobs=6).fit(xx[tr],y[tr]);return m.predict(xx[te])
def pair_matrix(a,b,X,D,bv,sel,dmean,dstd):
 da=X[b][:,sel]-X[a][:,sel];common=X[b][:,sel]*X[a][:,sel];dd=(D[b]-D[a])/dstd;mid=((D[b]+D[a])/2-dmean)/dstd;ss=np.asarray([DataStructs.TanimotoSimilarity(bv[i],bv[j]) for i,j in zip(a,b)])[:,None];return np.c_[da,common,dd,mid,ss].astype(np.float32)
def train_delta(train,X,D,y,g,bv,seed=42,cliff_weight=0.0):
 sel=np.argsort(X[train].var(0))[-512:];mean=D[train].mean(0);std=D[train].std(0);std[std<1e-8]=1;rng=np.random.default_rng(seed);aa=[];bb=[]
 for group in sorted(set(g[train])):
  idx=np.where((g==group)&np.isin(np.arange(len(g)),train))[0];pairs=np.asarray(list(itertools.combinations(idx,2)),int)
  if len(pairs)==0:continue
  if len(pairs)>2500:pairs=pairs[rng.choice(len(pairs),2500,replace=False)]
  aa.extend(pairs[:,0]);bb.extend(pairs[:,1]);aa.extend(pairs[:,1]);bb.extend(pairs[:,0])
 aa=np.asarray(aa,int);bb=np.asarray(bb,int);Z=pair_matrix(aa,bb,X,D,bv,sel,mean,std);target=y[bb]-y[aa]
 # Activity cliffs are scarce but scientifically decisive.  Up-weight only
 # high-similarity pairs with a >=1 log-unit potency jump; cliff_weight=0
 # recovers the original unweighted DeltaSAR model exactly.
 weights=1.0+float(cliff_weight)*((Z[:,-1]>=.5)&(np.abs(target)>=1.0))
 m=LGBMRegressor(n_estimators=400,learning_rate=.025,num_leaves=15,max_depth=6,min_child_samples=20,reg_alpha=.2,reg_lambda=2,random_state=seed,verbosity=-1,n_jobs=6).fit(Z,target,sample_weight=weights);return m,sel,mean,std,len(target)
def choose_diverse(target,k,bv):
 S=sims(target,target,bv);chosen=[int(target[np.argmax(S.mean(1))])]
 while len(chosen)<k:
  rest=[int(x) for x in target if x not in chosen];score=[min(1-DataStructs.TanimotoSimilarity(bv[x],bv[z]) for z in chosen) for x in rest];chosen.append(rest[int(np.argmax(score))])
 return np.asarray(chosen,int)
def choose_facility(target,k,bv,base=None):
 S=sims(target,target,bv)
 if base is not None:
  b=base[target];S=.75*S+.25*np.exp(-np.abs(b[:,None]-b[None,:])/.35)
 covered=np.zeros(len(target));chosen=[]
 for _ in range(k):
  gains=[-np.inf if j in chosen else np.maximum(covered,S[:,j]).sum()-covered.sum() for j in range(len(target))];j=int(np.argmax(gains));chosen.append(j);covered=np.maximum(covered,S[:,j])
 return target[np.asarray(chosen,int)]
def choose_predicted_span(target,k,base):
 order=target[np.argsort(base[target])];q=np.linspace(0,1,k+2)[1:-1] if k>1 else [.5];chosen=[]
 for z in q:
  i=int(round(z*(len(order)-1)));chosen.append(int(order[i]))
 return np.asarray(chosen,int)
def adapted(method,base,support,query,y,X,D,bv,delta_pack):
 if method=="Base":return base[query]
 S=sims(query,support,bv);res=y[support]-base[support];k=len(support)
 if method=="Offset":return base[query]+k/(k+2)*res.mean()
 if method=="KernelResidual":
  W=np.maximum(S,1e-4)**3;return base[query]+(W@res)/(W.sum(1)+1.5)
 if method=="TanimotoGP":
  K=sims(support,support,bv)**2;alpha=np.linalg.solve(K+.35*np.eye(k),res);return base[query]+(S**2)@alpha
 if method in {"DeltaSAR","DeltaSARHybrid"}:
  m,sel,mean,std,_=delta_pack;aa=np.repeat(support,len(query));bb=np.tile(query,len(support));Z=pair_matrix(aa,bb,X,D,bv,sel,mean,std);dp=m.predict(Z).reshape(k,len(query));W=np.maximum(S.T,1e-4)**3;absolute=((y[support,None]+dp)*W).sum(0)/(W.sum(0)+1e-8)
  if method=="DeltaSAR":return absolute
  lam=k/(k+2);return (1-lam)*base[query]+lam*absolute
 raise ValueError(method)
def safe_rho(y,p):return float(spearmanr(y,p).statistic) if len(y)>=3 and np.std(p)>0 else 0.
def main():
 d=pd.read_csv(DATA).reset_index(drop=True);X,D,bv=features(d.canonical_smiles);y=d.pEC50.to_numpy(float);g=d.source_component.astype(str).to_numpy();groups=[x for x,n in d.source_component.value_counts().items() if n>=MIN_GROUP];methods=["Base","Offset","KernelResidual","TanimotoGP","DeltaSAR","DeltaSARHybrid"];rows=[]
 for oi,group in enumerate(groups):
  target=np.where(g==group)[0];train=np.where(g!=group)[0];base=np.full(len(d),np.nan);base[target]=base_model(X,D,y,train,target);delta=train_delta(train,X,D,y,g,bv,42+oi);print(group,"delta_pairs",delta[-1],flush=True)
  # Same nested support permutations are shared by every method and shot.
  perms=[np.random.default_rng(10000+oi*100+r).permutation(target) for r in range(REPEATS)]
  for shot in SHOTS:
   policies=[("none",[np.asarray([],int)])] if shot==0 else [("random",[p[:shot] for p in perms]),("diverse",[choose_diverse(target,shot,bv)]),("facility",[choose_facility(target,shot,bv)]),("joint_facility",[choose_facility(target,shot,bv,base)]),("predicted_span",[choose_predicted_span(target,shot,base)])]
   for policy,supports in policies:
    for rep,support in enumerate(supports):
     query=np.asarray([x for x in target if x not in set(support)],int)
     for method in methods:
      if shot==0 and method!="Base":continue
      pred=base[query] if shot==0 else adapted(method,base,support,query,y,X,D,bv,delta)
      rows.append({"group":group,"n_group":len(target),"shot":shot,"policy":policy,"repeat":rep,"method":method,"n_query":len(query),"Spearman":safe_rho(y[query],pred),"MAE":float(mean_absolute_error(y[query],pred)),"support_ids":"|".join(d.loc[support,"canonical_molecule_id"])})
 r=pd.DataFrame(rows);r.to_csv(OUT/"episodes.csv",index=False);per=r.groupby(["shot","policy","method","group"])[["Spearman","MAE"]].mean().reset_index();summary=per.groupby(["shot","policy","method"]).agg(macro_Spearman=("Spearman","mean"),median_Spearman=("Spearman","median"),worst_Spearman=("Spearman","min"),positive_series=("Spearman",lambda x:int((x>0).sum())),n_series=("Spearman","size"),macro_MAE=("MAE","mean")).reset_index();summary.to_csv(OUT/"summary.csv",index=False);payload=summary.to_dict("records");(OUT/"metrics.json").write_text(json.dumps(payload,indent=2),encoding="utf-8");print(summary.to_string(index=False))
if __name__=="__main__":main()
