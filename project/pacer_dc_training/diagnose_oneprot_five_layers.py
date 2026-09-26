
"""G1: single-window, five-layer OneProt-MD diagnostic."""

import gc
import json
from pathlib import Path
import importlib.util

import numpy as np
import torch

ROOT = Path("project/results/pacer_dc_four_context_v01/compound110")
UNIT = ROOT / "replica_02/window_000"
CONTEXT = "candidate_probe"
CKPT = Path(
    "project/tools/oneprot-embeddings/artifacts/"
    "Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt"
)
G0_SCRIPT = Path(
    "project/pacer_dc_training/oneprot_g0_audit/"
    "scripts/run_g0_embedding.py"
)
OUT = ROOT / "five_layer_diagnostic_R2_W0_v01"

EXPECTED_SHA = (
    "a79bce2e8f91d9cd59840965e145a119c3202e30516e1e0f69452dee6f211997"
)

# 1. Check inputs and historical provenance.
summary = json.loads(
    (UNIT / "unit_summary.json").read_text(encoding="utf-8")
)
assert summary["checkpoint_sha256"] == EXPECTED_SHA
assert summary["window_frames"] == 100
assert summary["completeness"] == "complete"

atom_path = UNIT / "atom14" / f"{CONTEXT}_w000.npy"
csv_path = atom_path.with_suffix(".csv")
history_path = UNIT / "embeddings" / f"{CONTEXT}_w000.npy"

for path in (G0_SCRIPT, CKPT, atom_path, csv_path, history_path):
    if not path.is_file():
        raise FileNotFoundError(path)

# 2. Import exactly the same G0 preprocessing module.
spec = importlib.util.spec_from_file_location(
    "g0_five_layer_diagnostic", G0_SCRIPT
)
g0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g0)

arr = np.load(atom_path, allow_pickle=False)
seqres = csv_path.read_text(
    encoding="utf-8"
).strip().splitlines()[1].split(",")[1]

assert arr.shape == (100, 270, 14, 3)
assert len(seqres) == 270

latents, kwargs = g0.build_batch(arr, seqres)
assert tuple(latents.shape) == (1, 100, 270, 21)

print("INPUT_SHAPE:", tuple(latents.shape), flush=True)

# 3. Reconstruct the audited production encoder.
enc = g0.TrajectoryEncoder(
    output_dim=1024,
    hidden_size=21,
    num_layers=4,
    num_heads=8,
    pretrained=False,
    frozen=True,
    proj_type="mlp",
    num_frames=100,
    suffix="_i1",
)

checkpoint = torch.load(
    CKPT, map_location="cpu", weights_only=False
)
md = {
    k[len("network.md."):]: v
    for k, v in checkpoint["state_dict"].items()
    if k.startswith("network.md.")
}

incompat = enc.load_state_dict(md, strict=False)
missing = list(incompat.missing_keys)
unexpected = list(incompat.unexpected_keys)

print("MISSING_KEYS:", missing, flush=True)
print("UNEXPECTED_KEYS:", unexpected, flush=True)

assert not missing
assert set(unexpected) == {"norm.1.log_logit_scale"}

del checkpoint, md
gc.collect()
enc.eval()

# 4. Non-invasive hooks.
# L0 is the exact tensor entering MDGen's FinalLayer,
# not an assumed output of the last Transformer block.
captured = {}

def capture_pre(module, inputs):
    captured["L0"] = inputs[0].detach().cpu().clone()

def capture(name):
    def hook(module, inputs, output):
        if not isinstance(output, torch.Tensor):
            raise TypeError(f"{name}: unexpected output type")
        captured[name] = output.detach().cpu().clone()
    return hook

hooks = [
    enc.transformer.emb_to_latent.register_forward_pre_hook(
        capture_pre
    ),
    enc.transformer.emb_to_latent.register_forward_hook(
        capture("L1")
    ),
    enc.pooling.register_forward_hook(capture("L2")),
    enc.proj.register_forward_hook(capture("L3")),
    enc.norm.register_forward_hook(capture("L4")),
]

try:
    with torch.no_grad():
        first = enc(latents, 0, **kwargs)
finally:
    for hook in hooks:
        hook.remove()

# 5. Shape and finite-value checks.
expected = {
    "L0": (1, 100, 270, 384),
    "L1": (1, 100, 270, 21),
    "L2": (1, 21),
    "L3": (1, 1024),
    "L4": (1, 1024),
}

assert set(captured) == set(expected)

for name, shape in expected.items():
    x = captured[name]
    print(
        f"{name}: shape={tuple(x.shape)} "
        f"finite={bool(torch.isfinite(x).all())}",
        flush=True,
    )
    assert tuple(x.shape) == shape
    assert bool(torch.isfinite(x).all())

assert torch.equal(first.cpu(), captured["L4"])

# 6. Historical compatibility and repeatability.
historical = np.load(history_path, allow_pickle=False)
new = first[0].detach().cpu().numpy()

max_error = float(np.max(np.abs(new - historical)))
l2_error = float(np.linalg.norm(
    new.astype(np.float64) -
    historical.astype(np.float64)
))
compatible = bool(
    np.allclose(new, historical, rtol=1e-6, atol=1e-6)
)

with torch.no_grad():
    second = enc(latents, 0, **kwargs)

repeat_identical = bool(torch.equal(first, second))

print("HISTORICAL_MAX_ABS_ERROR:", max_error)
print("HISTORICAL_L2_ERROR:", l2_error)
print("HISTORICAL_COMPATIBLE:", compatible)
print("REPEAT_IDENTICAL:", repeat_identical)

passed = compatible and repeat_identical
print("G1_2F_PASS:", passed)

if not passed:
    raise RuntimeError(
        "Historical compatibility or repeatability failed. "
        "No layer files will be exported."
    )

# 7. Export only after validation.
OUT.mkdir(parents=True, exist_ok=True)

for name, tensor in captured.items():
    np.save(
        OUT / f"{name}_{CONTEXT}_R2_W0.npy",
        tensor.numpy()[0],
    )

report = {
    "analysis": "G1_five_layer_single_window_v01",
    "candidate": "compound110",
    "replica": 2,
    "window": 0,
    "context": CONTEXT,
    "checkpoint_sha256": EXPECTED_SHA,
    "layers": {
        name: {
            "shape": list(tensor.shape),
            "finite": bool(torch.isfinite(tensor).all()),
        }
        for name, tensor in captured.items()
    },
    "historical_max_abs_error": max_error,
    "historical_l2_error": l2_error,
    "historical_compatible": compatible,
    "repeat_identical": repeat_identical,
    "validation_pass": passed,
}

(OUT / "five_layer_report.json").write_text(
    json.dumps(report, indent=2),
    encoding="utf-8",
)

print("OUTPUT_DIR:", OUT)
print("EXPORTED_LAYERS:", len(captured))
print("G1_2F_EXPORT_PASS:", all(
    (OUT / f"{name}_{CONTEXT}_R2_W0.npy").is_file()
    for name in expected
))
