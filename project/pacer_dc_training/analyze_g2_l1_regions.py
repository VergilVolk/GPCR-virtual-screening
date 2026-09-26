
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("project/results/pacer_dc_four_context_v01/compound110")
OUTPUT = ROOT / "G2_L1_region_cross_replica_v01.csv"

CONTEXTS = (
    "candidate_probe",
    "candidate_no_probe",
    "probe_only",
    "apo",
)

regions = json.loads(
    (ROOT / "G2_REGION_MAP_v02.json").read_text(encoding="utf-8")
)["regions"]

g1 = pd.read_csv(ROOT / "five_layer_cross_replica_G1_v01.csv")
g1 = g1[g1["layer"] == "L1"]

assert len(g1) == 15
assert not OUTPUT.exists(), "Output exists; refusing overwrite"

rows = []

for w in range(5):
    data = {}

    for r in (2, 3):
        folder = ROOT / f"five_layer_diagnostic_R{r}_W{w}_v01"

        for c in CONTEXTS:
            path = folder / f"L1_{c}_R{r}_W{w}.npy"
            arr = np.load(path, mmap_mode="r")

            assert arr.shape == (100, 270, 21)
            assert np.isfinite(arr).all()
            data[r, c] = arr

    for region, residues in regions.items():
        ids = [x["embedding_index"] for x in residues]
        z = {}

        for r in (2, 3):
            for c in CONTEXTS:
                z[r, c] = np.asarray(
                    data[r, c][:, ids, :],
                    dtype=np.float64,
                ).mean(axis=(0, 1))

        contrasts = {}

        for r in (2, 3):
            pam = z[r, "candidate_probe"] - z[r, "probe_only"]
            ago = z[r, "candidate_no_probe"] - z[r, "apo"]

            contrasts[r] = {
                "dPAM": pam,
                "dAGO": ago,
                "dINT": pam - ago,
            }

        for metric in ("dPAM", "dAGO", "dINT"):
            a = contrasts[2][metric]
            b = contrasts[3][metric]

            n2 = float(np.linalg.norm(a))
            n3 = float(np.linalg.norm(b))

            assert n2 > 1e-12 and n3 > 1e-12

            distance = float(np.linalg.norm(a - b))
            cosine = float(np.dot(a, b) / (n2 * n3))

            rows.append({
                "window": w,
                "region": region,
                "metric": metric,
                "cosine": cosine,
                "norm_R2": n2,
                "norm_R3": n3,
                "distance": distance,
                "relative_distance": 2 * distance / (n2 + n3),
            })

result = pd.DataFrame(rows)

assert len(result) == 105
assert np.isfinite(
    result[
        ["cosine", "norm_R2", "norm_R3",
         "distance", "relative_distance"]
    ].to_numpy()
).all()

# Verify all 15 global L1 records against G1.
merged = result[result["region"] == "global"].merge(
    g1,
    on=["window", "metric"],
    suffixes=("_G2", "_G1"),
    validate="one_to_one",
)

assert len(merged) == 15

columns = [
    "cosine", "norm_R2", "norm_R3",
    "distance", "relative_distance",
]

errors = {
    col: float(np.max(np.abs(
        merged[f"{col}_G2"] - merged[f"{col}_G1"]
    )))
    for col in columns
}

print("ROWS:", len(result))
print("REGIONS:", result["region"].nunique())
print("WINDOWS:", result["window"].nunique())
print("MAX_ERRORS:", errors)

assert all(v < 1e-4 for v in errors.values()), errors

result.to_csv(OUTPUT, index=False)

print("G1_GLOBAL_MATCH: True")
print("OUTPUT:", OUTPUT)
print("ALL_PASS: True")
