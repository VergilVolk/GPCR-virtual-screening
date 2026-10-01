"""Tests for the manifest-driven frozen PACER-FKG runner.

Runnable under pytest or directly: `python project/tests/test_frozen_runner.py`.
No scientific computation is performed; every check is a contract, identity or
policy assertion against already-frozen artifacts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02.frozen_runner import blocks, engines, paths as path_policy
from project.pacer_fkg_v02.frozen_runner.spec import RunSpecError, authenticate_inputs, load_run_spec

LEDGER_PATH = REPO / "project/integration_freeze/FROZEN_MODULE_LEDGER_v01.json"
REGRESSION_SPEC = REPO / "project/integration_freeze/run_specs/CM00734_REGRESSION_v01.json"


def _expect_raises(exc, fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc as error:
        return error
    raise AssertionError(f"expected {exc.__name__} from {fn.__name__}")


# ---------------------------------------------------------------- frozen state
def test_frozen_numerics_match_freeze_manifest():
    report = engines.frozen_numerics_report()
    assert report["ok"], f"frozen numeric drift: {report['mismatches']}"


def test_frozen_module_ledger_matches_recorded():
    assert LEDGER_PATH.is_file(), f"missing ledger: {LEDGER_PATH}"
    recorded = json.loads(LEDGER_PATH.read_text(encoding="utf-8"))
    live = engines.frozen_module_ledger(REPO)
    for rel, entry in recorded["modules"].items():
        assert live[rel]["lf_sha256"] == entry["lf_sha256"], f"frozen module changed: {rel}"


def test_frozen_anchor_verifies():
    anchor = engines.verify_frozen_anchor()
    assert anchor["phase2_verification"]["status"] == "VERIFIED_FROZEN"
    assert anchor["phase2_verification"]["fitting_replicas"] == [2]
    assert anchor["phase2_verification"]["excluded_replicas"] == [1, 3]


# ------------------------------------------------------------------- run specs
def test_regression_run_spec_is_valid():
    spec = load_run_spec(REGRESSION_SPEC, REPO)
    assert spec.molecule.candidate_id == "CM00734"
    assert spec.replicas == (1, 2, 3)
    assert spec.seeds == {1: 27101, 2: 38201, 3: 49301}
    assert spec.blocks_per_trajectory * spec.block_frames == spec.n_frames
    assert spec.calibration["r2_fitting_allowed"] is False


def test_regression_spec_inputs_authenticate():
    spec = load_run_spec(REGRESSION_SPEC, REPO)
    report = authenticate_inputs(spec, REPO)
    assert report["all_match"] and report["count"] >= 7


def _mutate(tmp_path: Path, **changes) -> Path:
    payload = json.loads(REGRESSION_SPEC.read_text(encoding="utf-8"))
    for key, value in changes.items():
        if isinstance(value, dict) and isinstance(payload.get(key), dict):
            payload[key].update(value)
        else:
            payload[key] = value
    target = tmp_path / "mutated_spec.json"
    target.write_text(json.dumps(payload), encoding="utf-8")
    return target


def test_spec_rejects_wrong_seed_group(tmp_path):
    bad = _mutate(tmp_path, seeds={"1": 1, "2": 2, "3": 3})
    _expect_raises(RunSpecError, load_run_spec, bad, REPO)


def test_spec_rejects_forbidden_claim(tmp_path):
    bad = _mutate(tmp_path, scope={"pam_classifier_claim": True})
    _expect_raises(RunSpecError, load_run_spec, bad, REPO)


def test_spec_rejects_calibration_fitting(tmp_path):
    bad = _mutate(tmp_path, calibration={"usage": "fit", "r2_fitting_allowed": True})
    _expect_raises(RunSpecError, load_run_spec, bad, REPO)


def test_spec_rejects_wrong_block_geometry(tmp_path):
    bad = _mutate(tmp_path, frames={"n_frames": 400, "block_frames": 25, "blocks_per_trajectory": 16})
    _expect_raises(RunSpecError, load_run_spec, bad, REPO)


# --------------------------------------------------------------- path policy
def test_path_policy_rejects_protected_root():
    _expect_raises(
        Exception, path_policy.guard_new_output,
        REPO / "project/results/pacer_dc_membrane_reference_v01/x.npy",
        REPO, (REPO / "project/results",))


def test_path_policy_rejects_outside_allowed_roots():
    _expect_raises(
        Exception, path_policy.guard_new_output,
        REPO / "project/results/pacer_final_newmol_v01/x.npy",
        REPO, (REPO / "project/cache/pacer_final_newmol_v01",))


def test_path_policy_allows_declared_root():
    target = path_policy.guard_new_output(
        REPO / "project/results/pacer_final_newmol_v01/x.npy",
        REPO, (REPO / "project/results/pacer_final_newmol_v01",))
    assert target.name == "x.npy"


# ----------------------------------------------------------------- block maths
def test_blocks_bit_identical_to_stage_b_driver():
    rng = np.random.default_rng(7)
    frames = rng.standard_normal((400, engines.N_RESIDUES, engines.P2_BS_WIDTH), dtype=np.float32)
    sm, sd = blocks.construct_blocks(frames, block_frames=20, n_blocks=20)
    from project.pacer_fkg_v02.run_cm00734_stage_b_phase2a_frozen_apply_v01 import construct as stage_b
    reference = stage_b(frames)
    assert np.array_equal(sm, reference["STATE_MOTION"])
    assert np.array_equal(sd, reference["SIGNED_DRIFT"])


def test_blocks_endpoint_identity_stage_b():
    rng = np.random.default_rng(11)
    frames = rng.standard_normal((400, engines.N_RESIDUES, engines.P2_BS_WIDTH), dtype=np.float32)
    _, sd = blocks.construct_blocks(frames, block_frames=20, n_blocks=20)
    assert blocks.verify_endpoint_identity(frames, sd, block_frames=20, n_blocks=20)


def test_blocks_bit_identical_to_frozen_phase2():
    rng = np.random.default_rng(13)
    frames = rng.standard_normal(
        (engines.P2_N_FRAMES, engines.N_RESIDUES, engines.P2_BS_WIDTH), dtype=np.float32)
    sm, sd = blocks.construct_blocks(frames, block_frames=engines.P2_BLOCK_FRAMES, n_blocks=engines.P2_N_BLOCKS)
    ref_sm, ref_sd = engines.construct_blocks(frames)
    assert np.array_equal(sm, ref_sm)
    assert np.array_equal(sd, ref_sd)


def test_blocks_reject_wrong_shape():
    _expect_raises(ValueError, blocks.construct_blocks,
                   np.zeros((10, engines.N_RESIDUES, engines.P2_BS_WIDTH), dtype=np.float32),
                   block_frames=20, n_blocks=20)


def _run_all() -> int:
    namespace = {k: v for k, v in sorted(globals().items()) if k.startswith("test_")}
    import tempfile, traceback
    failures = 0
    for name, fn in namespace.items():
        kwargs = {}
        if "tmp_path" in fn.__code__.co_varnames[: fn.__code__.co_argcount]:
            kwargs["tmp_path"] = Path(tempfile.mkdtemp(prefix="frozen_runner_test_"))
        try:
            fn(**kwargs)
            print(f"PASS  {name}")
        except Exception:
            failures += 1
            print(f"FAIL  {name}")
            traceback.print_exc()
    print(f"\n{len(namespace) - failures}/{len(namespace)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(_run_all())
