from pathlib import Path
import hashlib
import json
import csv
import numpy as np

ROOT = Path(
    "project/results/pacer_dc_four_context_v01/compound110"
)
OUT = ROOT / "direction_audit_R2_R3_v01"

if OUT.exists():
    raise SystemExit("OUTPUT_EXISTS: refusing to overwrite")

def load(replica, window, name):
    path = (
        ROOT
        / f"replica_{replica:02d}"
        / f"window_{window:03d}"
        / "embeddings"
        / f"{name}_w{window:03d}.npy"
    )

    if not path.is_file():
        raise FileNotFoundError(path)

    raw = np.load(path, allow_pickle=False)

    if raw.shape != (1024,):
        raise ValueError(f"INVALID_SHAPE: {path}: {raw.shape}")

    if not np.isfinite(raw).all():
        raise ValueError(f"NONFINITE: {path}")

    sha = hashlib.sha256(path.read_bytes()).hexdigest()

    return raw.astype(np.float64), {
        "path": str(path),
        "sha256": sha,
        "shape": list(raw.shape),
        "dtype": str(raw.dtype),
    }

def norm(x):
    return float(np.linalg.norm(x))

def cosine(x, y):
    nx = norm(x)
    ny = norm(y)

    if nx <= 1e-12 or ny <= 1e-12:
        raise ValueError("ZERO_NORM_VECTOR")

    return float(np.dot(x, y) / (nx * ny))

def angle(x, y):
    return float(
        np.degrees(
            np.arccos(np.clip(cosine(x, y), -1.0, 1.0))
        )
    )

rows = []
provenance = []

for w in range(5):
    z = {}
    files = {}

    for r in (2, 3):
        z[r] = {}
        files[str(r)] = {}

        for key, name in (
            ("ca", "candidate_probe"),
            ("a", "probe_only"),
            ("d", "dPAM"),
        ):
            vector, info = load(r, w, name)
            z[r][key] = vector
            files[str(r)][key] = info

        residual = norm(
            z[r]["d"] - (z[r]["ca"] - z[r]["a"])
        )

        if residual >= 1e-5:
            raise ValueError(
                f"DIFFERENTIAL_MISMATCH R{r} W{w}: "
                f"{residual}"
            )

    delta_ca = z[3]["ca"] - z[2]["ca"]
    delta_a = z[3]["a"] - z[2]["a"]
    delta_d = z[3]["d"] - z[2]["d"]

    ca_sq = norm(delta_ca) ** 2
    a_sq = norm(delta_a) ** 2
    cross = float(-2 * np.dot(delta_ca, delta_a))

    predicted_sq = ca_sq + a_sq + cross
    observed_sq = norm(delta_d) ** 2
    residual = abs(predicted_sq - observed_sq)

    if residual >= 1e-5:
        raise ValueError(
            f"DECOMPOSITION_MISMATCH W{w}: {residual}"
        )

    row = {
        "window": w,
        "ca_cross_replica_cosine":
            cosine(z[2]["ca"], z[3]["ca"]),
        "a_cross_replica_cosine":
            cosine(z[2]["a"], z[3]["a"]),
        "dpam_cross_replica_cosine":
            cosine(z[2]["d"], z[3]["d"]),
        "dpam_cross_replica_angle_deg":
            angle(z[2]["d"], z[3]["d"]),
        "r2_dpam_norm": norm(z[2]["d"]),
        "r3_dpam_norm": norm(z[3]["d"]),
        "delta_ca_norm": norm(delta_ca),
        "delta_a_norm": norm(delta_a),
        "delta_dpam_norm": norm(delta_d),
        "delta_ca_a_cosine": cosine(delta_ca, delta_a),
        "ca_squared_contribution": ca_sq,
        "a_squared_contribution": a_sq,
        "cross_term": cross,
        "predicted_delta_dpam_sq": predicted_sq,
        "observed_delta_dpam_sq": observed_sq,
        "decomposition_residual": residual,
        "r2_differential_residual":
            norm(z[2]["d"] - (z[2]["ca"] - z[2]["a"])),
        "r3_differential_residual":
            norm(z[3]["d"] - (z[3]["ca"] - z[3]["a"])),
    }

    rows.append(row)
    provenance.append({
        "window": w,
        "files": files,
    })

    print(f"\n=== WINDOW {w} ===")
    print(
        "RAW_COSINES:",
        f"CA={row['ca_cross_replica_cosine']:.6f}",
        f"A={row['a_cross_replica_cosine']:.6f}"
    )
    print(
        "DELTA_CA_A_COSINE:",
        f"{row['delta_ca_a_cosine']:.6f}"
    )
    print("CROSS_TERM:", f"{cross:.6f}")
    print(
        "DPAM_COSINE:",
        f"{row['dpam_cross_replica_cosine']:.6f}"
    )
    print(
        "DPAM_ANGLE_DEG:",
        f"{row['dpam_cross_replica_angle_deg']:.3f}"
    )
    print(
        "DELTA_DPAM_NORM:",
        f"{row['delta_dpam_norm']:.6f}"
    )
    print(
        "DECOMPOSITION_RESIDUAL:",
        f"{residual:.10f}"
    )

# 全部计算和验证通过后才创建输出目录
OUT.mkdir(parents=True)

report = {
    "candidate": "compound110",
    "replicas": [2, 3],
    "windows": list(range(5)),
    "embedding_dim": 1024,
    "definition": "dPAM = candidate_probe - probe_only",
    "delta_definition": "R3 - R2",
    "results": rows,
    "provenance": provenance,
}

json_path = OUT / "dpam_direction_decomposition_R2_R3_v01.json"
csv_path = OUT / "dpam_direction_decomposition_R2_R3_v01.csv"

json_path.write_text(
    json.dumps(report, indent=2, ensure_ascii=False),
    encoding="utf-8"
)

with csv_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

print("\nALL_FIVE_WINDOWS_PASS:", len(rows) == 5)
print("JSON:", json_path)
print("CSV:", csv_path)
