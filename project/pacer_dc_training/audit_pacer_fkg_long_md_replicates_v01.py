from pathlib import Path
import csv, hashlib, json
import numpy as np

ROOT=Path(r"C:\projects\GPCR-virtual-screening")
RUN=ROOT/"project/results/pacer_dc_long_md_v02/fkg_common_kernel_v01"
VEC=RUN/"PACER_FKG_LONG_MD_VECTORS_v01.npz"
BLOCKCSV=RUN/"PACER_FKG_LONG_MD_BLOCKS_v01.csv"
SUMMARYCSV=RUN/"PACER_FKG_LONG_MD_REGION_AXIS_v01.csv"
HASHES=RUN/"PACER_FKG_LONG_MD_ARTIFACT_HASHES_v01.json"

OUTCSV=RUN/"PACER_FKG_LONG_MD_REPLICA_AUDIT_v01.csv"
OUTJSON=RUN/"PACER_FKG_LONG_MD_REPLICA_AUDIT_v01.json"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()

def cosine(a,b):
    a=np.asarray(a,dtype=np.float64)
    b=np.asarray(b,dtype=np.float64)
    den=np.linalg.norm(a)*np.linalg.norm(b)
    if den<=0: return np.nan
    return float(np.dot(a,b)/den)

def q(x,p):
    x=np.asarray(x,dtype=np.float64)
    x=x[np.isfinite(x)]
    return float(np.percentile(x,p)) if len(x) else np.nan

if OUTCSV.exists() or OUTJSON.exists():
    raise SystemExit("REFUSE_OVERWRITE")

# Verify the existing run before analysing it.
stored=json.loads(HASHES.read_text(encoding="utf-8"))
for p in (VEC,BLOCKCSV,SUMMARYCSV,RUN/"PACER_FKG_LONG_MD_RUN_AUDIT_v01.json"):
    if p.name not in stored:
        raise SystemExit(f"MISSING_STORED_HASH {p.name}")
    if sha(p)!=stored[p.name]:
        raise SystemExit(f"ARTIFACT_HASH_MISMATCH {p.name}")

with np.load(VEC,allow_pickle=False) as z:
    axes=[str(x) for x in z["axes"]]
    regions=[str(x) for x in z["regions"]]
    reps=[int(x) for x in z["replicas"]]
    block=np.asarray(z["block_vectors"],dtype=np.float64)
    pooled=np.asarray(z["pooled_vectors"],dtype=np.float64)

if block.shape!=(3,9,3,50,256):
    raise SystemExit(f"BAD_BLOCK_VECTOR_SHAPE {block.shape}")
if pooled.shape!=(3,9,3,256):
    raise SystemExit(f"BAD_POOLED_VECTOR_SHAPE {pooled.shape}")
if reps != [1,2,3]:
    raise SystemExit(f"BAD_REPLICAS {reps}")

stable_i=regions.index("stable_core_control")
rows=[]

for ai,axis in enumerate(axes):
    stable_block=block[ai,stable_i]
    stable_pool=pooled[ai,stable_i]

    for ri,region in enumerate(regions):
        rv=pooled[ai,ri]
        rb=block[ai,ri]

        norms=np.linalg.norm(rv,axis=1)
        stable_norms=np.linalg.norm(stable_pool,axis=1)
        pooled_gap=norms-stable_norms

        meanrep=float(norms.mean())
        meanvec=rv.mean(axis=0)
        normmean=float(np.linalg.norm(meanvec))
        coherence=float(normmean/meanrep) if meanrep>0 else np.nan

        c12=cosine(rv[0],rv[1])
        c13=cosine(rv[0],rv[2])
        c23=cosine(rv[1],rv[2])

        block_norm=np.linalg.norm(rb,axis=-1)
        stable_block_norm=np.linalg.norm(stable_block,axis=-1)
        block_gap=block_norm-stable_block_norm

        matched={}
        for name,i,j in (
            ("12",0,1),
            ("13",0,2),
            ("23",1,2),
        ):
            vals=[
                cosine(rb[i,b],rb[j,b])
                for b in range(50)
            ]
            matched[name]=np.asarray(vals,dtype=np.float64)

        row={
            "axis":axis,
            "region":region,
            "residues":None,

            "R1_pooled_norm":float(norms[0]),
            "R2_pooled_norm":float(norms[1]),
            "R3_pooled_norm":float(norms[2]),

            "R1_pooled_norm_minus_stable":float(pooled_gap[0]),
            "R2_pooled_norm_minus_stable":float(pooled_gap[1]),
            "R3_pooled_norm_minus_stable":float(pooled_gap[2]),

            "mean_replica_vector_norm":meanrep,
            "norm_of_mean_vector":normmean,
            "coherence_ratio":coherence,

            "R1_vs_R2_pooled_cosine":c12,
            "R1_vs_R3_pooled_cosine":c13,
            "R2_vs_R3_pooled_cosine":c23,
            "mean_pooled_pairwise_cosine":float(np.nanmean([c12,c13,c23])),
            "min_pooled_pairwise_cosine":float(np.nanmin([c12,c13,c23])),

            "R1_block_gap_median":q(block_gap[0],50),
            "R2_block_gap_median":q(block_gap[1],50),
            "R3_block_gap_median":q(block_gap[2],50),

            "R1_block_gap_positive_fraction":float(np.mean(block_gap[0]>0)),
            "R2_block_gap_positive_fraction":float(np.mean(block_gap[1]>0)),
            "R3_block_gap_positive_fraction":float(np.mean(block_gap[2]>0)),

            "R1R2_matched_block_cosine_median":q(matched["12"],50),
            "R1R2_matched_block_cosine_q025":q(matched["12"],2.5),
            "R1R2_matched_block_cosine_q975":q(matched["12"],97.5),

            "R1R3_matched_block_cosine_median":q(matched["13"],50),
            "R1R3_matched_block_cosine_q025":q(matched["13"],2.5),
            "R1R3_matched_block_cosine_q975":q(matched["13"],97.5),

            "R2R3_matched_block_cosine_median":q(matched["23"],50),
            "R2R3_matched_block_cosine_q025":q(matched["23"],2.5),
            "R2R3_matched_block_cosine_q975":q(matched["23"],97.5),

            "qualification_gate":"NOT_DEFINED",
        }
        rows.append(row)

with OUTCSV.open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

audit={
    "status":"PASS",
    "method":"descriptive three-replica common-RFF consistency audit",
    "input_vectors_sha256":sha(VEC),
    "input_blocks_sha256":sha(BLOCKCSV),
    "rows":len(rows),
    "axes":axes,
    "regions":regions,
    "independent_replication_unit":"replica_seed",
    "n_independent_replicas":3,
    "matched_temporal_blocks_per_replica":50,
    "block_duration_ns":1.0,
    "p_values_computed":False,
    "qualification_gate":"NOT_DEFINED",
    "automatic_biological_claims":False,
}
OUTJSON.write_text(json.dumps(audit,indent=2),encoding="utf-8")

print("PACER_FKG_REPLICA_CONSISTENCY_AUDIT_PASS")
print("rows",len(rows))
print("p_values_computed False")
print("qualification_gate NOT_DEFINED")
print()
print("axis region coherence meanCos minCos gapR1 gapR2 gapR3 posFracR1 posFracR2 posFracR3")

for r in rows:
    if r["region"]=="stable_core_control":
        continue
    print(
        f'{r["axis"]:4s} {r["region"]:28s} '
        f'{r["coherence_ratio"]:+.3f} '
        f'{r["mean_pooled_pairwise_cosine"]:+.3f} '
        f'{r["min_pooled_pairwise_cosine"]:+.3f} '
        f'{r["R1_pooled_norm_minus_stable"]:+.4f} '
        f'{r["R2_pooled_norm_minus_stable"]:+.4f} '
        f'{r["R3_pooled_norm_minus_stable"]:+.4f} '
        f'{r["R1_block_gap_positive_fraction"]:.2f} '
        f'{r["R2_block_gap_positive_fraction"]:.2f} '
        f'{r["R3_block_gap_positive_fraction"]:.2f}'
    )

print()
print("audit_csv_sha256",sha(OUTCSV))
print("audit_json_sha256",sha(OUTJSON))
