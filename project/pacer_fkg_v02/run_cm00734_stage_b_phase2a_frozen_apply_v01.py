#!/usr/bin/env python
from __future__ import annotations
import argparse,hashlib,json,os,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

REPO=Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path: sys.path.insert(0,str(REPO))
from project.pacer_fkg_v02.run_phase3_fkg import REQUIRED_FREEZE_SHA256,load_frozen_state,normalize_frozen,rff_frozen,verify_frozen_anchor

VERSION="PACER_DC_CM00734_STAGE_B_PHASE2A_FROZEN_APPLY_v01"
ANALYSIS_ROOT=REPO/"project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01"
PHASE1=ANALYSIS_ROOT/"track_b_fkg_v02/phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json"
CACHE_ROOT=REPO/"project/cache/pacer_dc_cm00734_stage_b_20ns_v01/track_b_fkg_v02/phase2a_frozen_apply"
REPORT_ROOT=ANALYSIS_ROOT/"track_b_fkg_v02/phase2a_frozen_apply"
RECEIPT=REPORT_ROOT/"PHASE2A_FREEZE_RECEIPT_v01.json"
SYSTEMS=("CM00734__candidate_no_probe","CM00734__candidate_probe")
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
def atomic_npy(p,a):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); t=p.with_name(f".{p.name}.partial.{os.getpid()}")
    with t.open("wb") as h: np.save(h,np.asarray(a),allow_pickle=False)
    os.replace(t,p)
def artifact(p):
    p=Path(p); return {"path":p.relative_to(REPO).as_posix(),"sha256":sha256(p),"bytes":p.stat().st_size}

def construct(frames):
    if frames.shape!=(400,270,256): raise ValueError(frames.shape)
    b=np.asarray(frames,dtype=np.float32).reshape(20,20,270,256)
    mu=b.mean(axis=1,dtype=np.float64).astype(np.float32)
    d=np.diff(b,axis=1)
    sd=d.mean(axis=1,dtype=np.float64).astype(np.float32)
    rms=np.sqrt(np.mean(np.square(d,dtype=np.float64),axis=1)).astype(np.float32)
    sm=np.concatenate((mu,rms),axis=-1).astype(np.float32,copy=False)
    endpoint=(b[:,-1]-b[:,0])/np.float32(19)
    if not np.allclose(sd,endpoint,rtol=2e-5,atol=2e-6): raise RuntimeError("SIGNED_DRIFT endpoint identity failed")
    return {"STATE_MOTION":sm,"SIGNED_DRIFT":sd}

def phase1_records():
    r=load_json(PHASE1)
    if r.get("status")!="PHASE1_FROZEN" or r.get("cache_count")!=6: raise RuntimeError("invalid Phase1 receipt")
    rec={(x["system"],int(x["replica"])):x for x in r["caches"]}
    out=[]
    for s in SYSTEMS:
        for rep in REPLICAS:
            x=rec[(s,rep)]; q=REPO/x["cache"]["path"]
            if sha256(q)!=x["cache"]["sha256"]: raise RuntimeError(f"Phase1 cache hash mismatch {s} R{rep}")
            out.append((s,rep,q,x["cache"]))
    return out

def run():
    if RECEIPT.exists(): raise FileExistsError(f"refusing overwrite: {RECEIPT}")
    anchor=verify_frozen_anchor(); outputs=[]
    for s,rep,q,crec in phase1_records():
        print(f"=== PHASE2A {s} R{rep} ===",flush=True)
        frames=np.load(q,mmap_mode="r",allow_pickle=False); desc=construct(frames)
        for branch in BRANCHES:
            state=load_frozen_state(anchor,branch)
            mapped=rff_frozen(normalize_frozen(desc[branch],state),state)
            if mapped.shape!=(20,270,512) or not np.isfinite(mapped).all(): raise RuntimeError("RFF contract failed")
            ddir=CACHE_ROOT/s/f"replica_{rep:02d}"; dp=ddir/f"{branch}.npy"; rp=ddir/f"{branch}_RFF.npy"
            atomic_npy(dp,desc[branch]); atomic_npy(rp,mapped)
            outputs.append({"system":s,"replica":rep,"branch":branch,"input_cache":crec,
                            "descriptor":artifact(dp),"rff":artifact(rp)})
    receipt={"schema":"pacer_dc.cm00734.stage_b.phase2a_freeze_receipt.v1","version":VERSION,"created_at":now(),
             "status":"PHASE2A_FROZEN","historical_fkg_v02_freeze_sha256":REQUIRED_FREEZE_SHA256,
             "phase1_receipt_sha256":sha256(PHASE1),"outputs":outputs,
             "blocks_per_trajectory":20,"block_frames":20,"frame_spacing_ps":50,"block_duration_ns":1.0,
             "forbidden_operations_performed":{"calibration":False,"normalization_refit":False,"bandwidth_refit":False,
                 "rff_refit":False,"graph_refit":False,"region_refit":False,"threshold_selection":False,"outcome_driven_tuning":False}}
    atomic_json(RECEIPT,receipt)
    return {"status":"PHASE2A_FROZEN","receipt":str(RECEIPT),"receipt_sha256":sha256(RECEIPT),"outputs":len(outputs)}

def verify():
    verify_frozen_anchor(); r=load_json(RECEIPT)
    if r.get("status")!="PHASE2A_FROZEN" or len(r.get("outputs",[]))!=12: raise RuntimeError("invalid Phase2A receipt")
    for x in r["outputs"]:
        for key in ("descriptor","rff"):
            q=REPO/x[key]["path"]
            if not q.is_file() or sha256(q)!=x[key]["sha256"] or q.stat().st_size!=x[key]["bytes"]: raise RuntimeError(f"artifact verification failed {q}")
    return {"status":"PHASE2A_VERIFIED","receipt_sha256":sha256(RECEIPT),"outputs":12}

def main():
    ap=argparse.ArgumentParser(); g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run",action="store_true"); g.add_argument("--verify",action="store_true")
    a=ap.parse_args(); print(json.dumps(run() if a.run else verify(),indent=2))
if __name__=="__main__": main()
