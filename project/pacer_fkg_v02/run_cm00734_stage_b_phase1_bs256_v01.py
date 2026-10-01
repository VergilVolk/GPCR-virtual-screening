#!/usr/bin/env python
from __future__ import annotations
import argparse,csv,hashlib,json,os,sys,time
from dataclasses import dataclass
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

REPO=Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path: sys.path.insert(0,str(REPO))

from project.pacer_fkg_v02.run_phase1_bs256 import (
    BACKBONE_ATOMS,BATCH_SIZE,BS_WIDTH,N_RESIDUES,RESTYPE_3TO1,
    infer_bs256,load_model,replay_metrics,topology_mapping,
    verify_phase0_and_frozen_candidate,
)
from project.pacer_fkg_v02.run_phase3_fkg import REQUIRED_FREEZE_SHA256,verify_frozen_anchor
from project.encoder_intermediate_v01.run_phase_0_2 import DEFAULT_CHECKPOINT

VERSION="PACER_DC_CM00734_STAGE_B_PHASE1_BS256_v01"
SOURCE_REPO=Path(os.environ.get("PACER_STAGE_B_SOURCE_REPO",r"C:\projects\GPCR-virtual-screening"))
SOURCE_RESULT=SOURCE_REPO/"project/results/pacer_dc_cm00734_stage_b_20ns_v01"
PRODUCTION_ROOT=SOURCE_RESULT/"production"
SOURCE_TOPOLOGY_ROOT=SOURCE_REPO/"project/results/pacer_dc_membrane_reference_v01"
ANALYSIS_ROOT=REPO/"project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01"
REPORT_ROOT=ANALYSIS_ROOT/"track_b_fkg_v02/phase1_bs256"
TOPOLOGY_ROOT=ANALYSIS_ROOT/"topology_sanitized"
CACHE_ROOT=REPO/"project/cache/pacer_dc_cm00734_stage_b_20ns_v01/track_b_fkg_v02/bs256"
RECEIPT=REPORT_ROOT/"PHASE1_FREEZE_RECEIPT_v01.json"
MAPPING=REPO/"project/results/pacer_dc_four_context_v01/compound110/G2_RESIDUE_MAPPING_v01.csv"
SYSTEMS=("CM00734__candidate_no_probe","CM00734__candidate_probe")
REPLICAS=(1,2,3); SEEDS={1:27101,2:38201,3:49301}
N_FRAMES=400; REPLAY_FRAME_IDS=(0,1,2,3)

def now(): return datetime.now(timezone.utc).astimezone().isoformat()
def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def load_json(path): return json.loads(Path(path).read_text(encoding="utf-8-sig"))
def atomic_json(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(f".{path.name}.partial.{os.getpid()}")
    tmp.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8"); os.replace(tmp,path)
def close_memmap(a):
    m=getattr(a,"_mmap",None) if a is not None else None
    if m is not None: m.close()

def frozen_sequence():
    with MAPPING.open(newline="",encoding="utf-8-sig") as h: rows=list(csv.DictReader(h))
    if len(rows)!=N_RESIDUES: raise RuntimeError("frozen residue mapping count mismatch")
    seq=[]
    for row in rows:
        aa=RESTYPE_3TO1.get(row["resname"].strip())
        if aa is None: raise RuntimeError(f"unsupported residue {row['resname']}")
        seq.append(aa)
    return "".join(seq)

def sanitize_topologies():
    TOPOLOGY_ROOT.mkdir(parents=True,exist_ok=True); rec={}
    for system in SYSTEMS:
        src=SOURCE_TOPOLOGY_ROOT/system/"minimized.pdb"
        if not src.is_file(): raise FileNotFoundError(src)
        dst=TOPOLOGY_ROOT/f"{system}_minimized_no_conect.pdb"
        content="".join(x for x in src.read_text(encoding="utf-8",errors="ignore").splitlines(True) if not x.startswith("CONECT"))
        if dst.exists():
            if dst.read_text(encoding="utf-8",errors="ignore")!=content: raise RuntimeError(f"existing sanitized topology differs: {dst}")
        else: dst.write_text(content,encoding="utf-8")
        rec[system]={"source_path":str(src),"source_sha256":sha256(src),"sanitized_path":str(dst),"sanitized_sha256":sha256(dst)}
    return rec

@dataclass(frozen=True)
class Job:
    system:str; replica:int; topology:Path; trajectory:Path; trajectory_sha256:str; topology_sha256:str; sequence:str
    @property
    def key(self): return f"{self.system}__replica_{self.replica:02d}"
    @property
    def frame_ids(self): return tuple(range(N_FRAMES))
    @property
    def output_shape(self): return (N_FRAMES,N_RESIDUES,BS_WIDTH)
    @property
    def cache_path(self): return CACHE_ROOT/self.system/f"replica_{self.replica:02d}"/"C1_BS256.npy"
    @property
    def manifest_path(self): return REPORT_ROOT/"manifests"/f"{self.key}.json"

def build_jobs(topos):
    seq=frozen_sequence(); jobs=[]
    for system in SYSTEMS:
        top=Path(topos[system]["sanitized_path"]); th=topos[system]["sanitized_sha256"]
        for rep in REPLICAS:
            d=PRODUCTION_ROOT/system/f"replica_{rep:02d}"
            traj=d/"trajectory.dcd"; prog=d/"progress.json"; state=d/"state.csv"
            for q in (traj,prog,state):
                if not q.is_file(): raise FileNotFoundError(q)
            pj=load_json(prog)
            if pj.get("status")!="complete" or float(pj.get("completed_ns",-1))!=20.0: raise RuntimeError(f"incomplete production {system} R{rep}")
            if int(pj.get("seed",-1))!=SEEDS[rep]: raise RuntimeError(f"seed mismatch {system} R{rep}")
            rows=sum(1 for _ in state.open(encoding="utf-8"))-1
            if rows!=N_FRAMES: raise RuntimeError(f"frame count mismatch {system} R{rep}: {rows}")
            jobs.append(Job(system,rep,top,traj,sha256(traj),th,seq))
    return jobs

def extract_one(job,model,device,model_record):
    if job.cache_path.exists() or job.manifest_path.exists(): raise FileExistsError(f"refusing overwrite {job.key}")
    job.cache_path.parent.mkdir(parents=True,exist_ok=True); job.manifest_path.parent.mkdir(parents=True,exist_ok=True)
    tmp=job.cache_path.with_name(f".{job.cache_path.name}.partial.{os.getpid()}")
    u=out=cache=None; started=time.perf_counter()
    try:
        u,atom_indices,atom_z,residue_index,bb,sc,records=topology_mapping(job)
        out=np.lib.format.open_memmap(tmp,mode="w+",dtype=np.float32,shape=job.output_shape)
        for start in range(0,N_FRAMES,BATCH_SIZE):
            ids=job.frame_ids[start:start+BATCH_SIZE]
            out[start:start+len(ids)]=infer_bs256(model,u,atom_indices,atom_z,residue_index,bb,sc,ids,device)
            out.flush(); done=start+len(ids)
            print(f"[{job.key}] {done:4d}/{N_FRAMES} ({100*done/N_FRAMES:6.2f}%)",flush=True)
        gly=np.asarray([aa=="G" for aa in job.sequence])
        if np.count_nonzero(out[:,gly,128:])!=0: raise RuntimeError("Gly side-chain output not zero")
        close_memmap(out); out=None
        cache=np.load(tmp,mmap_mode="r",allow_pickle=False)
        observed=np.asarray(cache[list(REPLAY_FRAME_IDS)])
        replayed=infer_bs256(model,u,atom_indices,atom_z,residue_index,bb,sc,REPLAY_FRAME_IDS,device)
        replay=replay_metrics(observed,replayed); replay["frame_ids"]=list(REPLAY_FRAME_IDS)
        if not replay["passed"]: raise RuntimeError(f"deterministic replay failed {job.key}")
        close_memmap(cache); cache=None
        digest=sha256(tmp); os.replace(tmp,job.cache_path)
        manifest={
            "schema":"pacer_dc.cm00734.stage_b.phase1_bs256.v1","version":VERSION,"created_at":now(),
            "system":job.system,"replica":job.replica,"shape":list(job.output_shape),"dtype":"float32",
            "frame_ids":list(job.frame_ids),"sequence":job.sequence,"selection":"chainID E and protein",
            "preprocessing":{"make_whole":False,"center":False,"align":False,"pbc_transform":False},
            "historical_fkg_v02_freeze_sha256":REQUIRED_FREEZE_SHA256,
            "trajectory":{"path":str(job.trajectory),"sha256":job.trajectory_sha256,"bytes":job.trajectory.stat().st_size},
            "topology":{"path":str(job.topology),"sha256":job.topology_sha256,"bytes":job.topology.stat().st_size},
            "cache":{"path":job.cache_path.relative_to(REPO).as_posix(),"sha256":digest,"bytes":job.cache_path.stat().st_size},
            "model_provenance":model_record,
            "atom_mapping":{"heavy_atoms":int(len(atom_indices)),"backbone_atoms":list(BACKBONE_ATOMS),
                            "sidechain_present_residues":int(sum(r["sidechain_present"] for r in records)),
                            "glycine_count":job.sequence.count("G")},
            "replay":replay,"elapsed_seconds":time.perf_counter()-started,
            "output_scope":"C1-BS256 only; no calibration or biological metrics",
        }
        atomic_json(job.manifest_path,manifest); return manifest
    finally:
        close_memmap(cache); close_memmap(out)
        if tmp.exists(): tmp.unlink()
        if u is not None: u.trajectory.close()

def run(checkpoint,device,geom2vec_source):
    if RECEIPT.exists(): raise FileExistsError(f"refusing overwrite: {RECEIPT}")
    frozen_candidate=verify_phase0_and_frozen_candidate(); verify_frozen_anchor()
    topos=sanitize_topologies(); jobs=build_jobs(topos)
    model,model_record=load_model(checkpoint,device,geom2vec_source)
    records=[]
    for i,j in enumerate(jobs,1):
        print(f"=== STAGE-B PHASE1 {i}/6 {j.key} ===",flush=True)
        records.append(extract_one(j,model,device,model_record))
    receipt={
        "schema":"pacer_dc.cm00734.stage_b.phase1_freeze_receipt.v1","version":VERSION,"created_at":now(),
        "status":"PHASE1_FROZEN","cache_count":6,"historical_fkg_v02_freeze_sha256":REQUIRED_FREEZE_SHA256,
        "source_result_root":str(SOURCE_RESULT),"topologies":topos,
        "caches":[{"system":m["system"],"replica":m["replica"],"trajectory":m["trajectory"],"topology":m["topology"],"cache":m["cache"],"replay":m["replay"]} for m in records],
        "forbidden_operations_performed":{"calibration":False,"normalization_refit":False,"bandwidth_refit":False,"rff_refit":False,
            "graph_refit":False,"region_refit":False,"threshold_selection":False,"outcome_driven_tuning":False},
        "frozen_candidate_anchor":frozen_candidate,
    }
    atomic_json(RECEIPT,receipt)
    return {"status":"PHASE1_FROZEN","receipt":str(RECEIPT),"receipt_sha256":sha256(RECEIPT),"cache_count":6}

def verify():
    verify_frozen_anchor(); r=load_json(RECEIPT)
    if r.get("status")!="PHASE1_FROZEN" or r.get("cache_count")!=6: raise RuntimeError("invalid Phase1 receipt")
    for item in r["caches"]:
        q=REPO/item["cache"]["path"]
        if not q.is_file() or sha256(q)!=item["cache"]["sha256"] or q.stat().st_size!=item["cache"]["bytes"]: raise RuntimeError(f"cache verification failed {q}")
        a=np.load(q,mmap_mode="r",allow_pickle=False)
        if a.shape!=(400,270,256) or a.dtype!=np.float32 or not np.isfinite(a).all(): raise RuntimeError(f"cache contract failed {q}")
    return {"status":"PHASE1_VERIFIED","receipt_sha256":sha256(RECEIPT),"cache_count":6}

def main():
    ap=argparse.ArgumentParser(); g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--run",action="store_true"); g.add_argument("--verify",action="store_true")
    ap.add_argument("--device",default="cuda"); ap.add_argument("--checkpoint",type=Path,default=DEFAULT_CHECKPOINT)
    ap.add_argument("--geom2vec-source",type=Path,default=Path(r"C:\projects\geom2vec-source"))
    a=ap.parse_args(); print(json.dumps(run(a.checkpoint,a.device,a.geom2vec_source) if a.run else verify(),indent=2))
if __name__=="__main__": main()
