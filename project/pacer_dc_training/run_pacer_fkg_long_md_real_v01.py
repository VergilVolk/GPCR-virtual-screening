from __future__ import annotations

from pathlib import Path
import argparse, csv, hashlib, json, sys
import numpy as np

from pacer_fkg_long_md_math_v02 import (
    AXIS_SIGNS,
    preprocess_features,
    rff_map,
    temporal_block_mean,
    signed_block_contrast,
    region_block_vectors,
    pooled_replica_vector,
    summarize_three_replicas,
)

EXPECTED_CONTRACT_SHA="0fc3121a6e105f3d4e4ddf2f09dfa99daf5152907b89f5b91b2cd0366bc0efc0"
EXPECTED_WEIGHT_SHA="9df1418a0825600e8890a0ce5f50c6886472bf9031f30d82fd34f496b94d9061"
EXPECTED_BIAS_SHA="c6a1e79365173d9ba3db0ecaa264d283c760d04755480e30ef0ea4183cbfeddc"

SYSTEM_CONTEXT={
    "apo":"apo",
    "probe_only":"probe",
    "compound110__candidate_no_probe":"candidate",
    "compound110__candidate_probe":"candidate_probe",
}
REPLICAS=(1,2,3)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def arrsha(a):
    a=np.ascontiguousarray(a)
    h=hashlib.sha256()
    h.update(str(a.dtype).encode("ascii"))
    h.update(str(tuple(a.shape)).encode("ascii"))
    h.update(a.tobytes(order="C"))
    return h.hexdigest()

def require(x,msg):
    if not x: raise RuntimeError(msg)

def basis(width,seed,dim):
    rng=np.random.default_rng(seed)
    w=rng.normal(0,1/width,size=(128,dim)).astype(np.float32)
    b=rng.uniform(0,2*np.pi,size=dim).astype(np.float32)
    return w,b

def graph_region(graph,name):
    return np.asarray(
        [int(m["embedding_index"]) for m in graph["regions"][name]],
        dtype=np.int64
    )

def write_csv(path,rows):
    if path.exists(): raise FileExistsError(f"REFUSE_OVERWRITE: {path}")
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def run(a):
    root=a.repo_root.resolve()
    out=a.output_root.resolve()
    if out.exists(): raise FileExistsError(f"REFUSE_OVERWRITE: {out}")
    out.mkdir(parents=True)

    contract=json.loads(a.contract.read_text(encoding="utf-8"))
    manifest=json.loads(a.feature_manifest.read_text(encoding="utf-8"))
    freeze=json.loads(a.math_freeze.read_text(encoding="utf-8"))
    graph=json.loads(a.graph.read_text(encoding="utf-8"))

    require(sha(a.contract)==EXPECTED_CONTRACT_SHA,"CONTRACT_HASH_MISMATCH")
    require(freeze["status"]=="PASS","MATH_FREEZE_NOT_PASS")
    require(freeze["contract_v03_sha256"]==EXPECTED_CONTRACT_SHA,"FREEZE_CONTRACT_MISMATCH")
    require(sha(Path(freeze["math_source"]))==freeze["math_source_sha256"],"MATH_SOURCE_CHANGED")
    require(sha(Path(freeze["synthetic_test_source"]))==freeze["synthetic_test_source_sha256"],"TEST_SOURCE_CHANGED")
    require(freeze["synthetic_test_replayed"]=="PASS","SYNTHETIC_NOT_FROZEN_PASS")
    require(freeze["kernel_u2_used"] is False,"KERNEL_U2_FORBIDDEN")
    require(freeze["qualification_gate"]=="NOT_DEFINED","GATE_FORBIDDEN")

    pre=contract["preprocessing"]
    ker=contract["common_kernel"]

    med=np.asarray(pre["channel_median"],dtype=np.float64)
    scale=np.asarray(pre["channel_scale"],dtype=np.float64)
    require(med.shape==(128,) and scale.shape==(128,),"BAD_PREPROCESS_SHAPE")
    require((scale>0).all(),"BAD_SCALE")
    require(pre["preprocessing_refit_on_long_md"] is False,"PREPROCESS_REFIT_FORBIDDEN")
    require(ker["bandwidth_refit_on_long_md"] is False,"WIDTH_REFIT_FORBIDDEN")
    require(ker["rff_basis_refit_on_long_md"] is False,"RFF_REFIT_FORBIDDEN")

    width=float(ker["shared_node_bandwidth"])
    dim=int(ker["rff_dimension"])
    seed=int(ker["actual_rff_rng_seed"])

    weights,bias=basis(width,seed,dim)
    require(arrsha(weights)==EXPECTED_WEIGHT_SHA,"RFF_WEIGHT_HASH_MISMATCH")
    require(arrsha(bias)==EXPECTED_BIAS_SHA,"RFF_BIAS_HASH_MISMATCH")

    entries={(e["system"],int(e["replica"])):e for e in manifest["entries"]}
    require(len(entries)==12,"EXPECTED_12_FEATURES")

    regions=sorted(graph["regions"])
    require("stable_core_control" in regions,"MISSING_STABLE_CORE")
    region_ix={r:graph_region(graph,r) for r in regions}
    stable_ix=region_ix["stable_core_control"]
    axes=("dAGO","dPAM","dINT")

    # Preprocess/RFF once per trajectory, reduce immediately to 50 temporal blocks.
    blocks={r:{} for r in REPLICAS}

    print("ENCODING_FROZEN_COMMON_RFF_BLOCKS")
    for rep in REPLICAS:
        for system,canonical in SYSTEM_CONTEXT.items():
            e=entries[(system,rep)]
            p=Path(e["path"])
            if not p.is_absolute(): p=(root/p).resolve()
            require(sha(p)==e["sha256"],f"FEATURE_HASH_MISMATCH {system} R{rep}")

            with np.load(p,allow_pickle=False) as z:
                raw=z["residue_features"]
                require(raw.shape==(1000,270,128),f"BAD_SHAPE {system} R{rep}")

                x=preprocess_features(raw,med,scale)
                mapped=rff_map(x,weights,bias)
                bmean=temporal_block_mean(mapped,20)

            require(bmean.shape==(50,270,256),"BAD_BLOCK_SHAPE")
            require(np.isfinite(bmean).all(),"NONFINITE_BLOCKS")
            blocks[rep][canonical]=bmean
            print(f"R{rep} {canonical:15s} PASS {bmean.shape}",flush=True)

    block_rows=[]
    summary_rows=[]

    block_store=np.empty(
        (len(axes),len(regions),3,50,256),
        dtype=np.float64
    )
    pooled_store=np.empty(
        (len(axes),len(regions),3,256),
        dtype=np.float64
    )

    for ai,axis in enumerate(axes):
        rep_node={}
        for rep in REPLICAS:
            rep_node[rep]=signed_block_contrast(blocks[rep],axis)
            require(rep_node[rep].shape==(50,270,256),"CONTRAST_SHAPE_FAIL")

        stable={
            rep:region_block_vectors(rep_node[rep],stable_ix)
            for rep in REPLICAS
        }

        for ri,region in enumerate(regions):
            ix=region_ix[region]
            vec={
                rep:region_block_vectors(rep_node[rep],ix)
                for rep in REPLICAS
            }
            pooled={
                rep:pooled_replica_vector(vec[rep])
                for rep in REPLICAS
            }

            for rep in REPLICAS:
                block_store[ai,ri,rep-1]=vec[rep]
                pooled_store[ai,ri,rep-1]=pooled[rep]

                n=np.linalg.norm(vec[rep],axis=1)
                s=np.linalg.norm(stable[rep],axis=1)
                gap=n-s

                for bi in range(50):
                    block_rows.append({
                        "axis":axis,
                        "region":region,
                        "replica":rep,
                        "block":bi,
                        "start_ns":float(bi),
                        "end_ns":float(bi+1),
                        "shared_rff_norm":float(n[bi]),
                        "stable_shared_rff_norm":float(s[bi]),
                        "shared_norm_minus_stable":float(gap[bi]),
                    })

            summ=summarize_three_replicas(pooled)

            rn=summ["replica_norms"]
            cos=summ["pairwise_direction_cosines"]

            row={
                "axis":axis,
                "region":region,
                "residues":len(ix),
                "R1_pooled_norm":rn[1],
                "R2_pooled_norm":rn[2],
                "R3_pooled_norm":rn[3],
                "mean_replica_vector_norm":summ["mean_replica_vector_norm"],
                "norm_of_mean_vector":summ["norm_of_mean_vector"],
                "R1_vs_R2_cosine":cos["R1_vs_R2"],
                "R1_vs_R3_cosine":cos["R1_vs_R3"],
                "R2_vs_R3_cosine":cos["R2_vs_R3"],
                "qualification_gate":"NOT_DEFINED",
            }
            summary_rows.append(row)

    write_csv(out/"PACER_FKG_LONG_MD_BLOCKS_v01.csv",block_rows)
    write_csv(out/"PACER_FKG_LONG_MD_REGION_AXIS_v01.csv",summary_rows)

    np.savez_compressed(
        out/"PACER_FKG_LONG_MD_VECTORS_v01.npz",
        axes=np.asarray(axes),
        regions=np.asarray(regions),
        replicas=np.asarray(REPLICAS,dtype=np.int64),
        block_vectors=block_store,
        pooled_vectors=pooled_store,
    )

    audit={
        "status":"PASS",
        "evidence_level":"descriptive_common_kernel_long_md",
        "contract_sha256":sha(a.contract),
        "math_freeze_sha256":sha(a.math_freeze),
        "feature_manifest_sha256":sha(a.feature_manifest),
        "graph_sha256":sha(a.graph),
        "shared_node_bandwidth":width,
        "rff_dimension":dim,
        "actual_rff_rng_seed":seed,
        "rff_weights_sha256":arrsha(weights),
        "rff_bias_sha256":arrsha(bias),
        "contexts":4,
        "paired_replicas":3,
        "blocks_per_trajectory":50,
        "block_duration_ns":1.0,
        "axes":list(axes),
        "regions":regions,
        "kernel_u2_gate":"REVOKED",
        "qualification_gate":"NOT_DEFINED",
        "p_values_computed":False,
        "automatic_biological_claims":False,
    }

    ap=out/"PACER_FKG_LONG_MD_RUN_AUDIT_v01.json"
    ap.write_text(json.dumps(audit,indent=2),encoding="utf-8")

    # Hash finished artifacts after all writes.
    artifacts={}
    for p in sorted(out.iterdir()):
        if p.is_file():
            artifacts[p.name]=sha(p)

    hp=out/"PACER_FKG_LONG_MD_ARTIFACT_HASHES_v01.json"
    hp.write_text(json.dumps(artifacts,indent=2),encoding="utf-8")

    print()
    print("PACER_FKG_LONG_MD_REAL_DESCRIPTIVE_PASS")
    print("axes",len(axes))
    print("regions",len(regions))
    print("independent_replicas",3)
    print("block_rows",len(block_rows))
    print("summary_rows",len(summary_rows))
    print("kernel_u2_gate REVOKED")
    print("qualification_gate NOT_DEFINED")
    print("p_values_computed False")
    print("automatic_biological_claims False")
    print("output",out)
    print()
    print("axis region R1norm R2norm R3norm cos12 cos13 cos23 normMean")
    for r in summary_rows:
        print(
            f'{r["axis"]:4s} {r["region"]:28s} '
            f'{r["R1_pooled_norm"]:.6f} '
            f'{r["R2_pooled_norm"]:.6f} '
            f'{r["R3_pooled_norm"]:.6f} '
            f'{r["R1_vs_R2_cosine"]:.4f} '
            f'{r["R1_vs_R3_cosine"]:.4f} '
            f'{r["R2_vs_R3_cosine"]:.4f} '
            f'{r["norm_of_mean_vector"]:.6f}'
        )

def cli():
    p=argparse.ArgumentParser()
    p.add_argument("--repo-root",type=Path,default=Path.cwd())
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--math-freeze",type=Path,required=True)
    p.add_argument("--feature-manifest",type=Path,required=True)
    p.add_argument("--graph",type=Path,required=True)
    p.add_argument("--output-root",type=Path,required=True)
    run(p.parse_args())

if __name__=="__main__":
    try: cli()
    except Exception as e:
        print(f"PACER_FKG_LONG_MD_REAL_FAIL: {type(e).__name__}: {e}")
        sys.exit(1)
