
from pathlib import Path
import numpy as np

ROOT = Path(
    "project/results/pacer_dc_four_context_v01/compound110"
)

HISTORICAL_DPAM_COS = [
    0.502744, 0.712309, -0.353644, 0.262047, 0.085430
]

def load_vector(replica, window, context):
    p = (
        ROOT
        / f"replica_{replica:02d}"
        / f"window_{window:03d}"
        / "embeddings"
        / f"{context}_w{window:03d}.npy"
    )
    x = np.load(p, allow_pickle=False).astype(np.float64)
    if x.shape != (1024,) or not np.isfinite(x).all():
        raise ValueError(f"Invalid embedding: {p}")
    return x

def cosine(a, b):
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if min(na, nb) < 1e-10:
        return float("nan")
    return float(np.clip(np.dot(a, b) / (na * nb), -1, 1))

def metrics(a, b):
    c = cosine(a, b)
    return {
        "cosine": c,
        "angle": float(np.degrees(np.arccos(c))),
        "distance": float(np.linalg.norm(a - b)),
        "norm_R2": float(np.linalg.norm(a)),
        "norm_R3": float(np.linalg.norm(b)),
    }

data = {}
for r in (2, 3):
    for w in range(5):
        ca = load_vector(r, w, "candidate_probe")
        c = load_vector(r, w, "candidate_no_probe")
        a = load_vector(r, w, "probe_only")
        apo = load_vector(r, w, "apo")

        dpam = ca - a
        dago = c - apo
        dint = ca - a - c + apo

        residual = np.linalg.norm(dint - (dpam - dago))
        if residual > 1e-10:
            raise RuntimeError("Difference identity failed")

        data[r, w] = {
            "dPAM": dpam,
            "dAGO": dago,
            "dINT": dint,
        }

print("\nG0-3: FIVE-WINDOW COMPARISON")
print("Window | Metric | Cosine | Angle_deg | L2_distance")

all_valid = True

for w in range(5):
    for name in ("dPAM", "dAGO", "dINT"):
        m = metrics(data[2, w][name], data[3, w][name])
        print(
            f"W{w} | {name} | "
            f"{m['cosine']:.6f} | "
            f"{m['angle']:.3f} | "
            f"{m['distance']:.6f}"
        )
        if not np.isfinite(m["cosine"]):
            all_valid = False

    observed = cosine(
        data[2, w]["dPAM"],
        data[3, w]["dPAM"],
    )
    if abs(observed - HISTORICAL_DPAM_COS[w]) > 0.001:
        print(f"WARNING: W{w} historical dPAM mismatch")
        all_valid = False

print("\nFIVE-WINDOW MEAN VECTOR COMPARISON")

for name in ("dPAM", "dAGO", "dINT"):
    mean2 = np.mean([data[2, w][name] for w in range(5)], axis=0)
    mean3 = np.mean([data[3, w][name] for w in range(5)], axis=0)
    m = metrics(mean2, mean3)
    print(
        f"{name}: cosine={m['cosine']:.6f}, "
        f"angle={m['angle']:.3f}, "
        f"distance={m['distance']:.6f}, "
        f"norm_R2={m['norm_R2']:.6f}, "
        f"norm_R3={m['norm_R3']:.6f}"
    )
    if not np.isfinite(m["cosine"]):
        all_valid = False

print("\nG0_3_PASS:", all_valid)


# G0-4: Export validated audit results
import csv
import json
import hashlib
from datetime import datetime, timezone

if not all_valid:
    raise RuntimeError("G0-3 validation failed; refusing export")

OUT = ROOT / "interaction_audit_R2_R3_v01"
OUT.mkdir(parents=True, exist_ok=True)

records = []
mean_records = []

for w in range(5):
    for name in ("dPAM", "dAGO", "dINT"):
        m = metrics(data[2, w][name], data[3, w][name])
        records.append({
            "window": w,
            "metric": name,
            "cosine": m["cosine"],
            "angle_deg": m["angle"],
            "l2_distance": m["distance"],
            "norm_R2": m["norm_R2"],
            "norm_R3": m["norm_R3"],
        })

for name in ("dPAM", "dAGO", "dINT"):
    a = np.mean(
        [data[2, w][name] for w in range(5)], axis=0
    )
    b = np.mean(
        [data[3, w][name] for w in range(5)], axis=0
    )
    m = metrics(a, b)
    mean_records.append({
        "metric": name,
        "cosine": m["cosine"],
        "angle_deg": m["angle"],
        "l2_distance": m["distance"],
        "norm_R2": m["norm_R2"],
        "norm_R3": m["norm_R3"],
    })

# Preserve input provenance for reproducibility.
input_hashes = {}

for r in (2, 3):
    for w in range(5):
        for name in (
            "candidate_probe",
            "candidate_no_probe",
            "probe_only",
            "apo",
        ):
            p = (
                ROOT
                / f"replica_{r:02d}"
                / f"window_{w:03d}"
                / "embeddings"
                / f"{name}_w{w:03d}.npy"
            )
            input_hashes[str(p)] = hashlib.sha256(
                p.read_bytes()
            ).hexdigest()

report = {
    "analysis": "compound110_R2_R3_interaction_audit_v01",
    "created_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "embedding_dimension": 1024,
    "replicas": [2, 3],
    "windows": [0, 1, 2, 3, 4],
    "definitions": {
        "dPAM": "z_CA - z_A",
        "dAGO": "z_C - z_0",
        "dINT": "z_CA - z_A - z_C + z_0",
    },
    "window_metrics": records,
    "mean_vector_metrics": mean_records,
    "input_sha256": input_hashes,
    "validation_pass": True,
    "limitations": [
        "compound110 only",
        "two independent replicas",
        "five windows are not independent biological samples",
        "dINT is an embedding interaction contrast, not proof of pharmacological synergy",
    ],
}

csv_path = OUT / "interaction_comparison_R2_R3_v01.csv"
json_path = OUT / "interaction_comparison_R2_R3_v01.json"

with csv_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "window",
        "metric",
        "cosine",
        "angle_deg",
        "l2_distance",
        "norm_R2",
        "norm_R3",
    ])
    writer.writeheader()
    writer.writerows(records)

with json_path.open("w", encoding="utf-8") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)

print("\nG0-4 EXPORT")
print("WINDOW_RECORDS:", len(records))
print("MEAN_RECORDS:", len(mean_records))
print("INPUT_HASHES:", len(input_hashes))
print("CSV:", csv_path)
print("JSON:", json_path)
print("G0_4_EXPORT_PASS:", (
    csv_path.is_file()
    and json_path.is_file()
    and len(records) == 15
    and len(mean_records) == 3
    and len(input_hashes) == 40
))
