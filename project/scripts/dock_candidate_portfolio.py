# -*- coding: utf-8 -*-
"""Dock the frozen PACER candidate portfolio into the ligand-free ACh M4 state.

Vina proposes a pose; it is not interpreted as PAM efficacy.  The promoted
structural readout is coverage of a predeclared allosteric-pocket residue set.
The run is resumable through a JSONL ledger.
"""
from __future__ import annotations
import argparse,hashlib,json,subprocess,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1]
DATA=P/"results"/"pacer_candidates_v01"/"predock_portfolio.csv"
STRUCT=P/"results"/"structure"/"ensemble"
RECEPTOR=STRUCT/"7TRS_R_meeko.pdbqt";RECEPTOR_PDB=STRUCT/"7TRS_R.pdb"
VINA=P/"tools"/"vina.exe";OUT=P/"results"/"pacer_candidates_v01"/"7trs_docking"
LEDGER=OUT/"features.jsonl";CSV=OUT/"features.csv";POSES=OUT/"poses"
CENTER=(110.2121366,107.8369051,68.1181557);SIZE=(22.,22.,22.)
POCKET=(89,92,93,96,184,186,190,423,432,433,435,436,439);CONTACT=4.5

def atoms(path):
 out=[];model=False;saw=False
 for line in path.read_text(encoding="utf-8",errors="ignore").splitlines():
  if line.startswith("MODEL"):saw=True;model=True;continue
  if line.startswith("ENDMDL"):break
  if not line.startswith(("ATOM","HETATM")) or (saw and not model):continue
  name=line[12:16].strip();element=line[76:79].strip().upper()
  if element in {"H","HD"} or name.startswith("H"):continue
  try:out.append((float(line[30:38]),float(line[38:46]),float(line[46:54]),int(line[22:26])))
  except ValueError:pass
 return out

REC=atoms(RECEPTOR_PDB);RECXYZ=np.asarray([x[:3] for x in REC],np.float32);RECNUM=np.asarray([x[3] for x in REC],np.int32)
MASK=np.isin(RECNUM,POCKET);PXYZ=RECXYZ[MASK];PNUM=RECNUM[MASK]

def pose_features(path):
 lig=np.asarray([x[:3] for x in atoms(path)],np.float32)
 if not len(lig):raise ValueError("empty pose")
 d=np.linalg.norm(lig[:,None,:]-PXYZ[None,:,:],axis=2);hit=d<CONTACT
 contacted=sorted({int(x) for x in PNUM[np.any(hit,axis=0)]});obs=set(contacted);ref=set(POCKET)
 feat={"n_heavy_atoms":len(lig),"pocket_atom_pair_contacts":int(hit.sum()),
       "ligand_atoms_contacting_pocket":int(np.any(hit,axis=1).sum()),
       "unique_pocket_residue_contacts":len(contacted),
       "pocket_residue_coverage":len(obs&ref)/len(ref),
       "pocket_jaccard":len(obs&ref)/len(obs|ref) if obs|ref else 0.,
       "box_centroid_distance_A":float(np.linalg.norm(lig.mean(0)-np.asarray(CENTER))),
       "contacted_residues":";".join(map(str,contacted))}
 for r in POCKET:
  feat[f"res{r}_min_A"]=float(d[:,PNUM==r].min());feat[f"res{r}_pairs"]=int(hit[:,PNUM==r].sum())
 return feat

def ligand(smiles,path,seed):
 from meeko import MoleculePreparation
 from meeko.preparation import PDBQTWriterLegacy
 m=Chem.AddHs(Chem.MolFromSmiles(smiles));p=AllChem.ETKDGv3();p.randomSeed=seed
 if AllChem.EmbedMolecule(m,p)!=0:raise ValueError("3D embedding failed")
 try:AllChem.MMFFOptimizeMolecule(m,maxIters=300)
 except Exception:pass
 setup=MoleculePreparation().prepare(m)[0];text,ok,error=PDBQTWriterLegacy.write_string(setup)
 if not ok:raise ValueError(str(error))
 path.write_text(text,encoding="utf-8")

def one(task):
 cid,smi,seed,exh=task;key=hashlib.sha1(f"{cid}|{seed}".encode()).hexdigest()[:14]
 lp=POSES/f"{key}.lig.pdbqt";pp=POSES/f"{key}.pose.pdbqt"
 try:
  if not pp.exists():
   ligand(smi,lp,seed)
   cmd=[str(VINA),"--receptor",str(RECEPTOR),"--ligand",str(lp),"--center_x",str(CENTER[0]),"--center_y",str(CENTER[1]),"--center_z",str(CENTER[2]),"--size_x",str(SIZE[0]),"--size_y",str(SIZE[1]),"--size_z",str(SIZE[2]),"--exhaustiveness",str(exh),"--num_modes","1","--cpu","1","--seed",str(seed),"--out",str(pp)]
   run=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
   if run.returncode:raise RuntimeError((run.stderr+run.stdout)[-500:])
  affinity=np.nan
  for line in pp.read_text(encoding="utf-8",errors="ignore").splitlines():
   if line.startswith("REMARK VINA RESULT:"):affinity=float(line.split()[3]);break
  return {"candidate_id":cid,"seed":seed,"vina_affinity":affinity,**pose_features(pp),"error":""}
 except Exception as e:return {"candidate_id":cid,"seed":seed,"error":repr(e)[:500]}

def main():
 ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=8);ap.add_argument("--exhaustiveness",type=int,default=4);ap.add_argument("--seeds",default="42");ap.add_argument("--limit",type=int,default=0);a=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);POSES.mkdir(exist_ok=True);df=pd.read_csv(DATA)
 if a.limit:df=df.head(a.limit)
 seeds=[int(x) for x in a.seeds.split(",")];old=[]
 if LEDGER.exists():old=[json.loads(x) for x in LEDGER.read_text(encoding="utf-8").splitlines() if x]
 done={(str(x["candidate_id"]),int(x["seed"])) for x in old};tasks=[(str(r.candidate_id),r.canonical_smiles,s,a.exhaustiveness) for r in df.itertuples() for s in seeds if (str(r.candidate_id),s) not in done]
 print(f"scheduled={len(tasks)} done={len(done)}",flush=True);t=time.time();n=0
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  for f in as_completed([pool.submit(one,x) for x in tasks]):
   row=f.result();n+=1
   with LEDGER.open("a",encoding="utf-8") as h:h.write(json.dumps(row,ensure_ascii=False)+"\n")
   if n%20==0 or n==len(tasks):print(f"completed={n}/{len(tasks)} elapsed_s={time.time()-t:.0f}",flush=True)
 rows=[json.loads(x) for x in LEDGER.read_text(encoding="utf-8").splitlines() if x];pd.DataFrame(rows).to_csv(CSV,index=False)
 meta={"state":"7TRS_ACh","n_candidates":len(df),"seeds":seeds,"exhaustiveness":a.exhaustiveness,"pocket_residues":POCKET,"interpretation":"Pose-proposal and pocket-compatibility gate only; not a PAM efficacy predictor."}
 (OUT/"metadata.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
if __name__=="__main__":main()
