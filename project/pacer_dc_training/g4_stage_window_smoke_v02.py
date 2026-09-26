import runpy
import sys
w = int(sys.argv[1])
assert 0 <= w <= 4
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

ROOT = Path(
    "project/results/pacer_dc_four_context_v01/compound110"
)

# Reuse the reconstruction validated across all 40 samples.
base = runpy.run_path(
    "project/pacer_dc_training/"
    "g4_l1_reconstruction_smoke_v01.py"
)

shift = base["shift"].detach().numpy().reshape(384).astype(np.float64)
scale = base["scale"].detach().numpy().reshape(384).astype(np.float64)
weight = base["lw"].detach().numpy().astype(np.float64)
bias = base["lb"].detach().numpy().astype(np.float64)
del base

import json

regions = json.loads(
    (ROOT / "G2_REGION_MAP_v02.json").read_text(
        encoding="utf-8"
    )
)["regions"]

ids = [
    x["embedding_index"]
    for x in regions["ACh_pocket"]
]

assert len(ids) == 12

contexts = (
    "candidate_probe",
    "probe_only",
    "candidate_no_probe",
    "apo",
)

z = {}
max_error = 0.0

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

        selected = np.array(
            l0[:, ids, :],
            dtype=np.float32,
            copy=True,
        )

        assert np.isfinite(selected).all()

        # S0: Original L0.
        s0 = selected.astype(np.float64).mean(
            axis=(0, 1)
        )

        # S1: Per-token LayerNorm, then regional pooling.
        with torch.no_grad():
            normalized = F.layer_norm(
                torch.from_numpy(selected),
                (384,),
                eps=1e-6,
            )

        s1 = normalized.numpy().astype(
            np.float64
        ).mean(axis=(0, 1))

        # S2: Shared conditional scaling.
        s2 = s1 * (1.0 + scale)

        # S3: Full affine reconstruction of pooled L1.
        s3 = (s2 + shift) @ weight.T + bias

        reference = np.asarray(
            l1[:, ids, :],
            dtype=np.float64,
        ).mean(axis=(0, 1))

        error = float(np.max(np.abs(s3 - reference)))
        max_error = max(max_error, error)

        assert error < 1e-4, (r, context, error)

        for name, vector in (
            ("S0", s0),
            ("S1", s1),
            ("S2", s2),
            ("S3", s3),
        ):
            assert np.isfinite(vector).all()
            z[name, r, context] = vector

def cosine(a, b):
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    assert na > 1e-12 and nb > 1e-12
    return float(np.dot(a, b) / (na * nb))

print(
    "STAGE DINT_COS NORM_R2 NORM_R3 "
    "REL_DISTANCE BIAS_COS CANCEL_RATIO"
)

for stage in ("S0", "S1", "S2", "S3"):
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
            "P": pam,
            "A": ago,
            "I": pam - ago,
        }

    a = contrasts[2]["I"]
    b = contrasts[3]["I"]

    n2 = float(np.linalg.norm(a))
    n3 = float(np.linalg.norm(b))

    distance = float(np.linalg.norm(a - b))
    relative = 2 * distance / (n2 + n3)
    dint_cos = cosine(a, b)

    p_bias = contrasts[2]["P"] - contrasts[3]["P"]
    a_bias = contrasts[2]["A"] - contrasts[3]["A"]
    i_bias = p_bias - a_bias

    bias_cos = cosine(p_bias, a_bias)
    cancel = float(
        np.linalg.norm(i_bias)
        / (
            np.linalg.norm(p_bias)
            + np.linalg.norm(a_bias)
        )
    )

    print(
        stage,
        *[
            round(x, 6)
            for x in (
                dint_cos, n2, n3,
                relative, bias_cos, cancel
            )
        ],
    )

    # Historical G2 consistency check.
    if stage in ("S0", "S3"):
        layer = "L0" if stage == "S0" else "L1"

        historical = pd.read_csv(
            ROOT / f"G2_{layer}_region_cross_replica_v01.csv"
        )

        match = historical[
            (historical["window"] == w)
            & (historical["region"] == "ACh_pocket")
            & (historical["metric"] == "dINT")
        ]

        assert len(match) == 1

        old = match.iloc[0]

        for column, value in (
            ("cosine", dint_cos),
            ("norm_R2", n2),
            ("norm_R3", n3),
            ("distance", distance),
            ("relative_distance", relative),
        ):
            assert abs(float(old[column]) - value) < 1e-4, (
                layer, column, old[column], value
            )

        print(f"G2_{layer}_MATCH: PASS")

print("MAX_POOLED_L1_ERROR:", max_error)
print(f"G4_STAGE_W{w}_QC: PASS")
