import runpy
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(
    "project/results/pacer_dc_four_context_v01/compound110"
)

SMOKE = Path(
    "project/pacer_dc_training/"
    "g4_l1_reconstruction_smoke_v01.py"
)

# Reuse the previously validated reconstruction.
base = runpy.run_path(str(SMOKE))

shift = base["shift"].reshape(1, 1, 384)
scale = base["scale"].reshape(1, 1, 384)
weight = base["lw"]
bias = base["lb"]

del base

contexts = (
    "candidate_probe",
    "probe_only",
    "candidate_no_probe",
    "apo",
)

checked = 0
global_max_error = 0.0

for r in (2, 3):
    for w in range(5):
        folder = ROOT / (
            f"five_layer_diagnostic_R{r}_W{w}_v01"
        )

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

            sample_max = 0.0

            # Process 10 frames at a time.
            for start in range(0, 100, 10):
                raw = np.array(
                    l0[start:start + 10],
                    dtype=np.float32,
                    copy=True,
                )
                reference = torch.from_numpy(
                    np.array(
                        l1[start:start + 10],
                        dtype=np.float32,
                        copy=True,
                    )
                )

                assert np.isfinite(raw).all()
                assert torch.isfinite(reference).all()

                with torch.no_grad():
                    x = torch.from_numpy(raw)
                    x = F.layer_norm(
                        x, (384,), eps=1e-6
                    )
                    x = (
                        x * (1 + scale)
                        + shift
                    )
                    reconstructed = F.linear(
                        x, weight, bias
                    )

                    error = torch.abs(
                        reconstructed - reference
                    )

                    assert torch.isfinite(error).all()

                    maximum = float(error.max())
                    sample_max = max(
                        sample_max, maximum
                    )

                    assert maximum < 1e-4, (
                        r, w, context, start, maximum
                    )

                    assert torch.allclose(
                        reconstructed,
                        reference,
                        rtol=1e-5,
                        atol=1e-4,
                    )

            checked += 1
            global_max_error = max(
                global_max_error, sample_max
            )

            print(
                f"R{r} W{w} {context}: "
                f"MAX_ERROR={sample_max:.9g}",
                flush=True,
            )

assert checked == 40

print("CHECKED:", checked)
print("GLOBAL_MAX_ERROR:", global_max_error)
print("G4_L1_MATRIX_QC: PASS")
