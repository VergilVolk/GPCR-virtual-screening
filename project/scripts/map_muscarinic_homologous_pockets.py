#!/usr/bin/env python3
"""Map the frozen M4 extracellular allosteric residue panel to M1/M2/M3/M5."""
from __future__ import annotations
import argparse, json
from pathlib import Path

AA={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
PANEL=[89,92,93,96,184,186,190,423,432,433,435,436,439]

def residues(path,chain):
    seen=set(); out=[]
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[21].strip()!=chain: continue
        name=line[17:20].strip(); key=(line[22:26].strip(),line[26].strip())
        if name not in AA or key in seen: continue
        seen.add(key); out.append({"resid":int(key[0]),"icode":key[1],"aa":AA[name],"resname":name})
    return out

def align(a,b,match=2,mismatch=-1,gap=-2):
    n,m=len(a),len(b); score=[[0]*(m+1) for _ in range(n+1)]; trace=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): score[i][0]=i*gap; trace[i][0]=1
    for j in range(1,m+1): score[0][j]=j*gap; trace[0][j]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            choices=(score[i-1][j-1]+(match if a[i-1]==b[j-1] else mismatch),score[i-1][j]+gap,score[i][j-1]+gap)
            score[i][j]=max(choices); trace[i][j]=choices.index(score[i][j])
    pairs=[]; i,j=n,m
    while i or j:
        t=trace[i][j]
        if i and j and t==0: pairs.append((i-1,j-1)); i-=1; j-=1
        elif i and (j==0 or t==1): pairs.append((i-1,None)); i-=1
        else: pairs.append((None,j-1)); j-=1
    return list(reversed(pairs)),score[n][m]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--pdb-dir",type=Path,required=True); ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args(); specs={"M1":("5CXV","A"),"M2":("4MQS","A"),"M3":("4U15","A"),"M4":("7TRQ","R"),"M5":("6OL9","A")}
    seq={k:residues(args.pdb_dir/f"{p}.pdb",c) for k,(p,c) in specs.items()}; ref=seq["M4"]
    ref_index={r["resid"]:i for i,r in enumerate(ref)}; missing=set(PANEL)-set(ref_index)
    if missing: raise ValueError(f"M4 panel missing from structure: {sorted(missing)}")
    report={"reference":"M4 7TRQ chain R","reference_panel":PANEL,"targets":{}}
    for subtype,(pdb,chain) in specs.items():
        pairs,score=align("".join(r["aa"] for r in ref),"".join(r["aa"] for r in seq[subtype]))
        amap={i:j for i,j in pairs if i is not None and j is not None}; mapped=[]
        for resid in PANEL:
            i=ref_index[resid]; j=amap.get(i)
            if j is None: mapped.append({"m4_resid":resid,"target_resid":None,"identity":False})
            else: mapped.append({"m4_resid":resid,"m4_aa":ref[i]["aa"],"target_resid":seq[subtype][j]["resid"],
                                "target_aa":seq[subtype][j]["aa"],"identity":ref[i]["aa"]==seq[subtype][j]["aa"]})
        report["targets"][subtype]={"pdb":pdb,"chain":chain,"alignment_score":score,
            "resolved_residues":len(seq[subtype]),"mapped_panel":mapped,
            "target_resids":[r["target_resid"] for r in mapped if r["target_resid"] is not None],
            "identity_fraction":sum(r.get("identity",False) for r in mapped)/len(mapped)}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(report,indent=2),encoding="utf-8"); print(json.dumps(report,indent=2))
if __name__=="__main__": main()
