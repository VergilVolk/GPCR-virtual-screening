"""G4: frozen-region, four-stage mechanism analysis."""
import json
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

ROOT = Path(
    "project/results/pacer_dc_four_context_v01/compound110"
)
OUTPUT = ROOT / "G4_four_stage_regions_v01.csv"

assert not OUTPUT.exists(), f"Output exists: {OUTPUT}"

base = runpy.run_path(
    "project/pacer_dc_training/"
    "g4_l1_reconstruction_smoke_v01.py"
)

shift = base["shift"].detach().numpy().reshape(384).astype(np.float64)
scale = base["scale"].detach().numpy().reshape(384).astype(np.float64)
weight = base["lw"].detach().numpy().astype(np.float64)
bias = base["lb"].detach().numpy().astype(np.float64)
del base

all_regions = json.loads(
    (ROOT / "G2_REGION_MAP_v02.json").read_text(
        encoding="utf-8"
    )
)["regions"]

region_counts = {
    "ACh_pocket": 12,
    "compound110_pocket": 12,
    "distal_control": 12,
    "global": 270,
}

regions = {}
for name, count in region_counts.items():
    ids = [x["embedding_index"] for x in all_regions[name]]
    assert len(ids) == count
    assert len(set(ids)) == count
    assert all(0 <= i < 270 for i in ids)
    regions[name] = ids

contexts = (
    "candidate_probe",
    "probe_only",
    "candidate_no_probe",
    "apo",
)
stages = ("S0", "S1", "S2", "S3")
metrics = ("dPAM", "dAGO", "dINT")
columns = (
    "cosine",
    "norm_R2",
    "norm_R3",
    "distance",
    "relative_distance",
)

historical = {
    layer: pd.read_csv(
        ROOT / f"G2_{layer}_region_cross_replica_v01.csv"
    ).set_index(["window", "region", "metric"], verify_integrity=True)
    for layer in ("L0", "L1")
}

def vector_metrics(a, b):
    n2 = float(np.linalg.norm(a))
    n3 = float(np.linalg.norm(b))
    assert n2 > 1e-12 and n3 > 1e-12
    distance = float(np.linalg.norm(a - b))
    return {
        "cosine": float(np.dot(a, b) / (n2 * n3)),
        "norm_R2": n2,
        "norm_R3": n3,
        "distance": distance,
        "relative_distance": 2 * distance / (n2 + n3),
    }

rows = []
history_checks = 0
maximum_l1_error = 0.0
maximum_g2_error = 0.0

for w in range(5):
    for region, ids in regions.items():
        z = {}
        region_max_error = 0.0

        for r in (2, 3):
            folder = ROOT / f"five_layer_diagnostic_R{r}_W{w}_v01"

            for context in contexts:
                stem = f"{context}_R{r}_W{w}"
                l0 = np.load(
                    folder / f"L0_{stem}.npy",
                    mmap_mode="r",
                    allow_pickle=False,
                )
                l1 = np.load(
                    folder / f"L1_{stem}.npy",
                    mmap_mode="r",
                    allow_pickle=False,
                )

                assert l0.shape == (100, 270, 384)
                assert l1.shape == (100, 270, 21)

                raw = np.array(
                    l0[:, ids, :],
                    dtype=np.float32,
                    copy=True,
                )
                reference = np.asarray(
                    l1[:, ids, :],
                    dtype=np.float64,
                ).mean(axis=(0, 1))

                assert np.isfinite(raw).all()
                assert np.isfinite(reference).all()

                s0 = raw.astype(np.float64).mean(axis=(0, 1))

                with torch.no_grad():
                    normalized = F.layer_norm(
                        torch.from_numpy(raw),
                        (384,),
                        eps=1e-6,
                    )

                s1 = normalized.numpy().astype(
                    np.float64
                ).mean(axis=(0, 1))
                s2 = s1 * (1.0 + scale)
                s3 = (s2 + shift) @ weight.T + bias

                reconstruction_error = float(np.max(np.abs(s3 - reference)))
                maximum_l1_error = max(maximum_l1_error, reconstruction_error)
                region_max_error = max(region_max_error, reconstruction_error)
                assert reconstruction_error < 1e-4, (
                    w, region, r, context, reconstruction_error
                )

                for stage, vector in zip(
                    stages, (s0, s1, s2, s3)
                ):
                    assert np.isfinite(vector).all()
                    z[stage, r, context] = vector

        for stage in stages:
            contrasts = {}

            for r in (2, 3):
                pam = (
                    z[stage, r, "candidate_probe"]
                    - z[stage, r, "probe_only"]
                )
                ago = (
                    z[stage, r, "candidate_no_probe"]
                    - z[stage, r, "apo"]
                )
                contrasts[r] = {
                    "dPAM": pam,
                    "dAGO": ago,
                    "dINT": pam - ago,
                }

            for metric in metrics:
                values = vector_metrics(
                    contrasts[2][metric],
                    contrasts[3][metric],
                )

                if stage in ("S0", "S3"):
                    layer = "L0" if stage == "S0" else "L1"
                    old = historical[layer].loc[
                        (w, region, metric)
                    ]

                    for column in columns:
                        error = abs(
                            float(old[column]) - values[column]
                        )
                        maximum_g2_error = max(
                            maximum_g2_error, error
                        )
                        assert error < 1e-4, (
                            layer, w, region, metric,
                            column, error
                        )
                    history_checks += 1

                if metric == "dINT":
                    p_bias = (
                        contrasts[2]["dPAM"]
                        - contrasts[3]["dPAM"]
                    )
                    a_bias = (
                        contrasts[2]["dAGO"]
                        - contrasts[3]["dAGO"]
                    )
                    bias_metrics = vector_metrics(
                        p_bias, a_bias
                    )
                    denominator = (
                        np.linalg.norm(p_bias)
                        + np.linalg.norm(a_bias)
                    )
                    assert denominator > 1e-12
                    cancel = float(
                        np.linalg.norm(p_bias - a_bias)
                        / denominator
                    )

                    rows.append({
                        "window": w,
                        "region": region,
                        "stage": stage,
                        **values,
                        "bias_cosine": bias_metrics["cosine"],
                        "cancel_ratio": cancel,
                        "l1_reconstruction_max_error": region_max_error
                        if stage == "S3" else np.nan,
                    })

        print(f"W{w} {region}: QC PASS", flush=True)

result = pd.DataFrame(rows)

assert len(result) == 80
assert history_checks == 120
assert not result.duplicated(
    ["window", "region", "stage"]
).any()

numeric = [
    "cosine", "norm_R2", "norm_R3",
    "distance", "relative_distance",
    "bias_cosine", "cancel_ratio",
]
assert np.isfinite(result[numeric].to_numpy()).all()
assert result["l1_reconstruction_max_error"].isna().sum() == 60

print("RESULT_ROWS:", len(result))
print("HISTORICAL_G2_CHECKS:", history_checks)
print("MAX_G2_ERROR:", maximum_g2_error)
print("MAX_L1_RECONSTRUCTION_ERROR:", maximum_l1_error)
print("G4_SUMMARY_QC: PASS")

result.to_csv(OUTPUT, index=False, mode="x")
print("OUTPUT:", OUTPUT)
