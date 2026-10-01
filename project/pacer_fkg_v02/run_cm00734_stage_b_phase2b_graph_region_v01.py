#!/usr/bin/env python
from __future__ import annotations
import argparse,hashlib,json,os,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

REPO=Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path: sys.path.insert(0,str(REPO))
from project.pacer_fkg_v02 import run_phase2_calibration as phase2
from project.pacer_fkg_v02.run_phase3_fkg import (
    REQUIRED_FREEZE_SHA256,GRAPH_PATH,CONTRASTS,HISTORICAL_AXIS,
    diffuse,graph_transition,historical_contrast,load_frozen_state,
    normalize_frozen,region_indices,rff_frozen,verify_frozen_anchor,
)

VERSION="PACER_DC_CM00734_STAGE_B_PHASE2B_GRAPH_REGION_v01"
ANALYSIS=REPO/"project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01"
PHASE2A=ANALYSIS/"track_b_fkg_v02/phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json"
CM_RFF_ROOT=REPO/"project/cache/pacer_dc_cm00734_stage_b_20ns_v01/track_b_fkg_v02/phase2a_frozen_apply"
OUTPUT=ANALYSIS/"track_b_fkg_v02/phase2b_graph_region"
RESULT=OUTPUT/"STAGE_B_MATCHED_20NS_RESULTS_v01.json"
SOURCE_REPO=Path(os.environ.get("PACER_STAGE_B_SOURCE_REPO",r"C:\projects\GPCR-virtual-screening"))
LY_RESULT=SOURCE_REPO/"project/results/pacer_dc_close_loop_20ns_v01/track_b_fkg_v02/phase2b_graph_region/TRACK_B_MATCHED_20NS_RESULTS_v01.json"

HIST={"A":"apo","P":"probe_only","C":"compound110__candidate_no_probe","CP":"compound110__candidate_probe"}
CM={"A":"apo","P":"probe_only","C":"CM00734__candidate_no_probe","CP":"CM00734__candidate_probe"}
REPLICAS=(1,2,3); BRANCHES=("STATE_MOTION","SIGNED_DRIFT")

def now(): return datetime.now(timezone.utc).astimezone().isoformat()
def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def load_json(p): return json.loads(Path(p).read_text(encoding="utf-8-sig"))
def atomic_json(p,obj):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); t=p.with_name(f".{p.name}.partial.{os.getpid()}")
    t.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8"); os.replace(t,p)
def atomic_npz(p,arrs):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); t=p.with_name(f".{p.name}.partial.{os.getpid()}.npz")
    np.savez_compressed(t,**arrs); os.replace(t,p)

def construct(frames):
    f=np.asarray(frames,dtype=np.float32)
    if f.shape!=(400,270,256): raise ValueError(f.shape)
    b=f.reshape(20,20,270,256)
    mu=b.mean(axis=1,dtype=np.float64).astype(np.float32)
    d=np.diff(b,axis=1)
    sd=d.mean(axis=1,dtype=np.float64).astype(np.float32)
    rms=np.sqrt(np.mean(np.square(d,dtype=np.float64),axis=1)).astype(np.float32)
    endpoint=(b[:,-1]-b[:,0])/np.float32(19)
    if not np.allclose(sd,endpoint,rtol=2e-5,atol=2e-6): raise RuntimeError("SIGNED_DRIFT endpoint identity failed")
    return {"STATE_MOTION":np.concatenate((mu,rms),axis=-1).astype(np.float32,copy=False),"SIGNED_DRIFT":sd}

def historical_rff(anchor,system,rep,branch):
    f=np.load(phase2.cache_path(system,rep),mmap_mode="r",allow_pickle=False)
    if f.shape!=(1000,270,256): raise RuntimeError(f"historical cache mismatch {system} R{rep}")
    d=construct(np.asarray(f[:400]))
    st=load_frozen_state(anchor,branch)
    out=rff_frozen(normalize_frozen(d[branch],st),st)
    if out.shape!=(20,270,512) or not np.isfinite(out).all(): raise RuntimeError("historical RFF contract failed")
    return out

def cm_rff(system,rep,branch):
    p=CM_RFF_ROOT/system/f"replica_{rep:02d}"/f"{branch}_RFF.npy"
    a=np.load(p,mmap_mode="r",allow_pickle=False)
    if a.shape!=(20,270,512) or not np.isfinite(a).all(): raise RuntimeError(f"CM RFF invalid {p}")
    return np.asarray(a)

def cosine(a,b):
    a=np.asarray(a,dtype=np.float64); b=np.asarray(b,dtype=np.float64)
    den=float(np.linalg.norm(a)*np.linalg.norm(b))
    return None if den<=1e-12 else float(np.clip(np.dot(a,b)/den,-1,1))

def summarize(v):
    v=np.asarray(v,dtype=np.float64)
    mean=v.mean(axis=0); mag=np.linalg.norm(v,axis=1)
    q10,q50,q90=np.quantile(mag,[0.10,0.50,0.90])
    return {"mean_vector_norm":float(np.linalg.norm(mean)),
            "block_magnitude_mean":float(mag.mean()),
            "block_magnitude_q10":float(q10),
            "block_magnitude_median":float(q50),
            "block_magnitude_q90":float(q90)}

def eval_panel(name,branch,rep,contexts,transition,regions):
    cr={}
    for alias,val in contexts.items():
        dv=diffuse(val,transition)
        cr[alias]={r:dv[:,idx,:].mean(axis=1,dtype=np.float64).astype(np.float32) for r,idx in regions.items()}
    ctr={}
    for c in CONTRASTS:
        dv=diffuse(historical_contrast(contexts,c),transition)
        ctr[c]={r:dv[:,idx,:].mean(axis=1,dtype=np.float64).astype(np.float32) for r,idx in regions.items()}
    arr={**{f"context__{a}__{r}":v for a,m in cr.items() for r,v in m.items()},
         **{f"contrast__{c}__{r}":v for c,m in ctr.items() for r,v in m.items()}}
    p=OUTPUT/branch.lower()/name/f"replica_{rep:02d}"/"BLOCK_REGIONAL_VECTORS.npz"
    atomic_npz(p,arr)
    return ctr,{"path":p.relative_to(REPO).as_posix(),"sha256":sha256(p),"bytes":p.stat().st_size}

def run():
    if RESULT.exists(): raise FileExistsError(f"refusing overwrite: {RESULT}")
    p2=load_json(PHASE2A)
    if p2.get("status")!="PHASE2A_FROZEN" or p2.get("historical_fkg_v02_freeze_sha256")!=REQUIRED_FREEZE_SHA256:
        raise RuntimeError("Phase2A anchor invalid")
    anchor=verify_frozen_anchor()
    graph=load_json(GRAPH_PATH); transition=graph_transition(graph); regions=region_indices(graph)
    allres={}; vec={}; arts=[]
    for branch in BRANCHES:
        allres[branch]={"compound110_reference":{},"CM00734":{}}
        vec[branch]={"compound110_reference":{},"CM00734":{}}
        for rep in REPLICAS:
            hist={a:historical_rff(anchor,s,rep,branch) for a,s in HIST.items()}
            cm={"A":hist["A"],"P":hist["P"],"C":cm_rff(CM["C"],rep,branch),"CP":cm_rff(CM["CP"],rep,branch)}
            for name,contexts in (("compound110_reference",hist),("CM00734",cm)):
                contrasts,art=eval_panel(name,branch,rep,contexts,transition,regions)
                vec[branch][name][rep]=contrasts
                arts.append({"branch":branch,"panel":name,"replica":rep,**art})
        for name in ("compound110_reference","CM00734"):
            panel={}
            for c in CONTRASTS:
                panel[c]={"historical_axis":HISTORICAL_AXIS[c],"regions":{}}
                for region in regions:
                    reps={f"R{rep}":summarize(vec[branch][name][rep][c][region]) for rep in REPLICAS}
                    r1=vec[branch][name][1][c][region].mean(axis=0)
                    r3=vec[branch][name][3][c][region].mean(axis=0)
                    panel[c]["regions"][region]={"replicas":reps,"R1_R3_direction_cosine":cosine(r1,r3)}
            allres[branch][name]=panel

    comp={}
    for branch in BRANCHES:
        comp[branch]={}
        for c in CONTRASTS:
            comp[branch][c]={}
            for region in regions:
                reps={}
                for rep in REPLICAS:
                    ref=vec[branch]["compound110_reference"][rep][c][region].mean(axis=0)
                    cm=vec[branch]["CM00734"][rep][c][region].mean(axis=0)
                    rn=float(np.linalg.norm(ref)); cn=float(np.linalg.norm(cm))
                    reps[f"R{rep}"]={"direction_cosine_CM_vs_compound110":cosine(cm,ref),
                                    "CM_mean_vector_norm":cn,"compound110_mean_vector_norm":rn,
                                    "magnitude_ratio_CM_over_compound110":cn/rn if rn>1e-12 else None}
                comp[branch][c][region]={"replicas":reps}

    ly_compare=None
    if LY_RESULT.is_file():
        ly=load_json(LY_RESULT); ly_compare={}
        for branch in BRANCHES:
            ly_compare[branch]={}
            for c in CONTRASTS:
                ly_compare[branch][c]={}
                for region in regions:
                    cmreg=allres[branch]["CM00734"][c]["regions"][region]
                    lyreg=ly["results"][branch]["LY2119620"][c]["regions"][region]
                    reps={}
                    for rep in REPLICAS:
                        k=f"R{rep}"; cn=cmreg["replicas"][k]["mean_vector_norm"]; ln=lyreg["replicas"][k]["mean_vector_norm"]
                        reps[k]={"CM_mean_vector_norm":cn,"LY_mean_vector_norm":ln,
                                 "magnitude_ratio_CM_over_LY":cn/ln if ln>1e-12 else None}
                    ly_compare[branch][c][region]={"replicas":reps,
                        "CM_R1_R3_direction_cosine":cmreg["R1_R3_direction_cosine"],
                        "LY_R1_R3_direction_cosine":lyreg["R1_R3_direction_cosine"]}

    report={
        "schema":"pacer_dc.cm00734.stage_b.matched20ns.graph_region.v1","version":VERSION,
        "created_at":now(),"status":"STAGE_B_MATCHED_20NS_COMPLETE",
        "upstream":{"phase2a_receipt_sha256":sha256(PHASE2A),"historical_fkg_v02_freeze_sha256":REQUIRED_FREEZE_SHA256},
        "matched_design":{"frames_per_trajectory":400,"frame_spacing_ps":50,"duration_ns":20,
            "block_frames":20,"blocks_per_trajectory":20,
            "shared_controls":{"A":"apo first 20 ns","P":"probe_only first 20 ns"},
            "compound110_reference":{"C":"compound110 no-probe first 20 ns","CP":"compound110 probe first 20 ns"},
            "CM00734":{"C":"CM00734 no-probe 20 ns","CP":"CM00734 probe 20 ns"}},
        "contrasts":{n:{"historical_axis":HISTORICAL_AXIS[n],"formula":f} for n,f in CONTRASTS.items()},
        "graph":{"path":GRAPH_PATH.relative_to(REPO).as_posix(),"sha256":sha256(GRAPH_PATH),
                 "diffusion_alpha":0.65,"diffusion_steps":20,
                 "regions":{n:int(len(i)) for n,i in regions.items()}},
        "results":allres,
        "CM_vs_compound110_descriptive_comparison":comp,
        "CM_vs_LY_descriptive_comparison":ly_compare,
        "LY_reference_path":str(LY_RESULT) if LY_RESULT.is_file() else None,
        "vector_artifacts":arts,
        "claim_boundary":"Matched-20ns descriptive hard-negative evaluation. No threshold tuning, classifier fitting, efficacy prediction, or independent-block inference.",
        "forbidden_operations_performed":{"calibration":False,"normalization_refit":False,"bandwidth_refit":False,
            "rff_refit":False,"graph_refit":False,"region_refit":False,"threshold_selection":False,"outcome_driven_tuning":False},
    }
    atomic_json(RESULT,report)
    return {"status":report["status"],"output":str(RESULT),"output_sha256":sha256(RESULT),
            "branches":list(BRANCHES),"panels":["compound110_reference","CM00734"],
            "CM_vs_LY_comparison_included":ly_compare is not None,"vector_artifacts":len(arts)}

def verify():
    r=load_json(RESULT)
    if r.get("status")!="STAGE_B_MATCHED_20NS_COMPLETE": raise RuntimeError("result status invalid")
    verify_frozen_anchor()
    for a in r["vector_artifacts"]:
        p=REPO/a["path"]
        if not p.is_file() or sha256(p)!=a["sha256"] or p.stat().st_size!=a["bytes"]:
            raise RuntimeError(f"vector artifact failed {p}")
    return {"status":"STAGE_B_MATCHED_20NS_VERIFIED","output_sha256":sha256(RESULT),
            "vector_artifacts":len(r["vector_artifacts"])}

def main():
    ap=argparse.ArgumentParser(); g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run",action="store_true"); g.add_argument("--verify",action="store_true")
    a=ap.parse_args(); print(json.dumps(run() if a.run else verify(),indent=2))
if __name__=="__main__": main()
