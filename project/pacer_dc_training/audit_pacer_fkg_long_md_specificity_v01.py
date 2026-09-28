from pathlib import Path
import csv, hashlib, json
import numpy as np

ROOT = Path(r"C:\projects\GPCR-virtual-screening")
RUN = ROOT / "project/results/pacer_dc_long_md_v02/fkg_common_kernel_v01"

VEC = RUN / "PACER_FKG_LONG_MD_VECTORS_v01.npz"
HASHES = RUN / "PACER_FKG_LONG_MD_ARTIFACT_HASHES_v01.json"

OUTCSV = RUN / "PACER_FKG_LONG_MD_SPECIFICITY_AUDIT_v01.csv"
OUTJSON = RUN / "PACER_FKG_LONG_MD_SPECIFICITY_AUDIT_v01.json"

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def cosine(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    den = np.linalg.norm(a) * np.linalg.norm(b)
    if den <= 0:
        return np.nan
    return float(np.dot(a, b) / den)

def med(x):
    return float(np.median(np.asarray(x, dtype=np.float64)))

if OUTCSV.exists() or OUTJSON.exists():
    raise SystemExit("REFUSE_OVERWRITE")

stored = json.loads(HASHES.read_text(encoding="utf-8"))
if stored.get(VEC.name) != sha(VEC):
    raise SystemExit("VECTOR_ARTIFACT_HASH_MISMATCH")

with np.load(VEC, allow_pickle=False) as z:
    axes = [str(x) for x in z["axes"]]
    regions = [str(x) for x in z["regions"]]
    reps = [int(x) for x in z["replicas"]]
    block = np.asarray(z["block_vectors"], dtype=np.float64)
    pooled = np.asarray(z["pooled_vectors"], dtype=np.float64)

if block.shape != (3, 9, 3, 50, 256):
    raise SystemExit(f"BAD_BLOCK_SHAPE {block.shape}")
if pooled.shape != (3, 9, 3, 256):
    raise SystemExit(f"BAD_POOLED_SHAPE {pooled.shape}")
if reps != [1, 2, 3]:
    raise SystemExit(f"BAD_REPLICAS {reps}")

stable_i = regions.index("stable_core_control")
distal_i = regions.index("distal_control")

targets = [
    r for r in regions
    if r not in ("stable_core_control", "distal_control")
]

rows = []

for ai, axis in enumerate(axes):
    stable_pool = pooled[ai, stable_i]
    distal_pool = pooled[ai, distal_i]

    stable_block = block[ai, stable_i]
    distal_block = block[ai, distal_i]

    stable_pool_norm = np.linalg.norm(stable_pool, axis=1)
    distal_pool_norm = np.linalg.norm(distal_pool, axis=1)

    stable_block_norm = np.linalg.norm(stable_block, axis=-1)
    distal_block_norm = np.linalg.norm(distal_block, axis=-1)

    for region in targets:
        ri = regions.index(region)

        rv = pooled[ai, ri]
        rb = block[ai, ri]

        pool_norm = np.linalg.norm(rv, axis=1)
        block_norm = np.linalg.norm(rb, axis=-1)

        excess_stable = pool_norm - stable_pool_norm
        excess_distal = pool_norm - distal_pool_norm

        gap_stable = block_norm - stable_block_norm
        gap_distal = block_norm - distal_block_norm

        cos_stable = [
            cosine(rv[i], stable_pool[i])
            for i in range(3)
        ]
        cos_distal = [
            cosine(rv[i], distal_pool[i])
            for i in range(3)
        ]

        row = {
            "axis": axis,
            "region": region,

            "R1_excess_vs_stable": float(excess_stable[0]),
            "R2_excess_vs_stable": float(excess_stable[1]),
            "R3_excess_vs_stable": float(excess_stable[2]),

            "R1_excess_vs_distal": float(excess_distal[0]),
            "R2_excess_vs_distal": float(excess_distal[1]),
            "R3_excess_vs_distal": float(excess_distal[2]),

            "replicas_positive_vs_stable": int(np.sum(excess_stable > 0)),
            "replicas_positive_vs_distal": int(np.sum(excess_distal > 0)),

            "R1_block_fraction_gt_stable": float(np.mean(gap_stable[0] > 0)),
            "R2_block_fraction_gt_stable": float(np.mean(gap_stable[1] > 0)),
            "R3_block_fraction_gt_stable": float(np.mean(gap_stable[2] > 0)),

            "R1_block_fraction_gt_distal": float(np.mean(gap_distal[0] > 0)),
            "R2_block_fraction_gt_distal": float(np.mean(gap_distal[1] > 0)),
            "R3_block_fraction_gt_distal": float(np.mean(gap_distal[2] > 0)),

            "R1_block_gap_vs_stable_median": med(gap_stable[0]),
            "R2_block_gap_vs_stable_median": med(gap_stable[1]),
            "R3_block_gap_vs_stable_median": med(gap_stable[2]),

            "R1_block_gap_vs_distal_median": med(gap_distal[0]),
            "R2_block_gap_vs_distal_median": med(gap_distal[1]),
            "R3_block_gap_vs_distal_median": med(gap_distal[2]),

            "R1_target_vs_stable_cosine": cos_stable[0],
            "R2_target_vs_stable_cosine": cos_stable[1],
            "R3_target_vs_stable_cosine": cos_stable[2],

            "R1_target_vs_distal_cosine": cos_distal[0],
            "R2_target_vs_distal_cosine": cos_distal[1],
            "R3_target_vs_distal_cosine": cos_distal[2],

            "mean_target_vs_stable_cosine": float(np.nanmean(cos_stable)),
            "mean_target_vs_distal_cosine": float(np.nanmean(cos_distal)),

            "qualification_gate": "NOT_DEFINED",
        }
        rows.append(row)

with OUTCSV.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

audit = {
    "status": "PASS",
    "method": "descriptive target-vs-stable-and-distal common-RFF specificity audit",
    "input_vectors_sha256": sha(VEC),
    "axes": axes,
    "target_regions": targets,
    "controls": ["stable_core_control", "distal_control"],
    "independent_replication_unit": "replica_seed",
    "n_independent_replicas": 3,
    "matched_blocks_per_replica": 50,
    "p_values_computed": False,
    "qualification_gate": "NOT_DEFINED",
    "automatic_biological_claims": False,
    "rows": len(rows),
}
OUTJSON.write_text(json.dumps(audit, indent=2), encoding="utf-8")

print("PACER_FKG_SPECIFICITY_AUDIT_PASS")
print("rows", len(rows))
print("controls stable_core_control distal_control")
print("p_values_computed False")
print("qualification_gate NOT_DEFINED")
print()
print(
    "axis region "
    "rep>stable rep>distal "
    "exS_R1 exS_R2 exS_R3 "
    "exD_R1 exD_R2 exD_R3 "
    "fracD_R1 fracD_R2 fracD_R3 "
    "cosStable cosDistal"
)

for r in rows:
    print(
        f'{r["axis"]:4s} {r["region"]:28s} '
        f'{r["replicas_positive_vs_stable"]:1d} '
        f'{r["replicas_positive_vs_distal"]:1d} '
        f'{r["R1_excess_vs_stable"]:+.4f} '
        f'{r["R2_excess_vs_stable"]:+.4f} '
        f'{r["R3_excess_vs_stable"]:+.4f} '
        f'{r["R1_excess_vs_distal"]:+.4f} '
        f'{r["R2_excess_vs_distal"]:+.4f} '
        f'{r["R3_excess_vs_distal"]:+.4f} '
        f'{r["R1_block_fraction_gt_distal"]:.2f} '
        f'{r["R2_block_fraction_gt_distal"]:.2f} '
        f'{r["R3_block_fraction_gt_distal"]:.2f} '
        f'{r["mean_target_vs_stable_cosine"]:+.3f} '
        f'{r["mean_target_vs_distal_cosine"]:+.3f}'
    )

print()
print("specificity_csv_sha256", sha(OUTCSV))
print("specificity_json_sha256", sha(OUTJSON))
