"""Step-by-step diagnostic for the CM00734 payload regression crash."""
from __future__ import annotations

import faulthandler
import sys
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _shim_helper import _Shim  # noqa: E402

faulthandler.enable()
REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))


def step(name):
    print(f"STEP {name}", flush=True)


try:
    step("import engines")
    from project.pacer_fkg_v02.frozen_runner import blocks, engines
    from project.pacer_fkg_v02.frozen_runner import fkg_phase2b as p2b

    step("frozen numerics")
    print("  ok =", engines.frozen_numerics_report()["ok"], flush=True)

    step("anchor (cached)")
    anchor = engines.cached_anchor()
    print("  phase2 =", anchor["phase2_verification"], flush=True)

    step("blocks equivalence small")
    rng = np.random.default_rng(13)
    frames400 = rng.standard_normal((400, 270, 256), dtype=np.float32)
    sm, sd = blocks.construct_blocks(frames400, block_frames=20, n_blocks=20)
    print("  shapes", sm.shape, sd.shape, flush=True)

    step("blocks equivalence 1000")
    frames1000 = rng.standard_normal((1000, 270, 256), dtype=np.float32)
    sm2, sd2 = blocks.construct_blocks(frames1000, block_frames=20, n_blocks=50)
    print("  shapes", sm2.shape, sd2.shape, flush=True)
    del frames1000, sm2, sd2

    step("graph + transition + regions")
    graph = engines.load_json(engines.GRAPH_PATH)
    transition = engines.graph_transition(graph)
    regions = engines.region_indices(graph)
    print("  regions", len(regions), flush=True)

    step("load frozen states")
    sm_state = engines.load_frozen_state(anchor, "STATE_MOTION")
    sd_state = engines.load_frozen_state(anchor, "SIGNED_DRIFT")
    print("  ok", flush=True)

    step("reference_branches one replica")
    desc = p2b.reference_branches(_Shim(), "apo", 1)
    print("  shapes", {k: v.shape for k, v in desc.items()}, flush=True)

    step("candidate rff load")
    arr = p2b.candidate_rff(_Shim(), "CM00734__candidate_no_probe", 1, "STATE_MOTION")
    print("  ", arr.shape, flush=True)

    step("diffuse one contrast")
    ctx = {"A": desc["STATE_MOTION"], "P": desc["STATE_MOTION"],
           "C": arr, "CP": arr}
    ctr = engines.historical_contrast(ctx, "Delta_INT")
    out = engines.diffuse(ctr, transition)
    print("  ", out.shape, np.isfinite(out).all(), flush=True)

    print("ALL STEPS OK", flush=True)
except Exception:
    traceback.print_exc()
    raise SystemExit(1)
