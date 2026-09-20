# -*- coding: utf-8 -*-
"""Prepare a preregistered M4 receptor-state ensemble (7TRQ/7TRP/7TRS).

7TRQ: VU0467154 + iperoxo; 7TRP: LY2033298 + iperoxo;
7TRS: acetylcholine, no PAM. 7TRS receives the allosteric box by a backbone
Kabsch transform from 7TRQ, not by ligand-label optimization.
"""
from __future__ import annotations
import json, subprocess
from pathlib import Path
import numpy as np

PROJECT=Path(__file__).resolve().parents[1]
PDBDIR=PROJECT/"data"/"pdb"; OUT=PROJECT/"results"/"structure"/"ensemble"
STATES={"7TRQ":{"pam":"IUI","probe":"IXO"},"7TRP":{"pam":"IUE","probe":"IXO"},
        "7TRS":{"pam":None,"probe":"ACH"}}
REF_CENTER=np.asarray([107.948,85.737,70.412],float)

def lines(pdb,kind="ATOM",chain="R",resname=None):
    out=[]
    for l in pdb.read_text(errors="ignore").splitlines(True):
        if not l.startswith(kind): continue
        if chain is not None and l[21]!=chain: continue
        if resname is not None and l[17:20].strip()!=resname: continue
        out.append(l)
    return out

def ca_map(pdb):
    return {int(l[22:26]):np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            for l in lines(pdb) if l[12:16].strip()=="CA"}

def transform(ref,target):
    common=sorted(set(ref)&set(target)); x=np.vstack([ref[i] for i in common]); y=np.vstack([target[i] for i in common])
    xc,yc=x.mean(0),y.mean(0); u,_,vt=np.linalg.svd((x-xc).T@(y-yc)); r=u@vt
    if np.linalg.det(r)<0: u[:,-1]*=-1; r=u@vt
    t=yc-xc@r; rms=float(np.sqrt(np.mean(np.sum((x@r+t-y)**2,axis=1))))
    return r,t,rms,len(common)

def centroid(pdb,resname):
    a=lines(pdb,kind="HETATM",chain=None,resname=resname)
    xyz=np.asarray([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in a
                    if (l[76:78].strip() or l[12:16].strip()[0]).upper()!="H"])
    return xyz.mean(0),len(xyz)

def main():
    OUT.mkdir(parents=True,exist_ok=True); ref=ca_map(PDBDIR/"7TRQ.pdb"); report={}
    for state,meta in STATES.items():
        pdb=PDBDIR/f"{state}.pdb"; receptor=OUT/f"{state}_R.pdb"; receptor.write_text("".join(lines(pdb)),encoding="utf-8")
        r,t,rms,n=transform(ref,ca_map(pdb)); mapped=REF_CENTER@r+t
        direct=None
        if meta["pam"]:
            direct,nlig=centroid(pdb,meta["pam"])
            ligand=OUT/f"{state}_{meta['pam']}.pdb"; ligand.write_text("".join(lines(pdb,"HETATM",None,meta["pam"])),encoding="utf-8")
        run=subprocess.run(["mk_prepare_receptor.exe","--read_pdb",str(receptor),"-o",str(OUT/f"{state}_R_meeko"),
                            "-p","--charge_model","gasteiger","-a"],capture_output=True,text=True)
        report[state]={"probe":meta["probe"],"pam":meta["pam"],"common_CA":n,"backbone_fit_RMSD_A":rms,
                       "mapped_7TRQ_box_center":mapped.tolist(),"native_PAM_centroid":direct.tolist() if direct is not None else None,
                       "center_disagreement_A":float(np.linalg.norm(mapped-direct)) if direct is not None else None,
                       "meeko_ok":run.returncode==0}
    (OUT/"state_ensemble_audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
