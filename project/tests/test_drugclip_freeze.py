"""Tests for the frozen DrugCLIP apply-only path.

Runnable under pytest or directly: `python project/tests/test_drugclip_freeze.py`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.drugclip_freeze import dual_model, model_a, model_b, projection

SHORTLIST = REPO / "project/results/project_wide_integration_benchmark_v02/candidate_dual_model_shortlist_v02.csv"


def _expect_raises(exc, fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc as error:
        return error
    raise AssertionError(f"expected {exc.__name__} from {fn.__name__}")


def test_frozen_identities_are_pinned():
    assert model_a.DRUGCLIP_COMMIT == "7a3a3fa33673f8668c811790f2e4681c98af44ef"
    assert model_a.UNICORE_COMMIT == "44f6386f4dcd7137fc1e5d5e768117d635d64a26"
    assert model_a.CHECKPOINT_BYTES == 1183713459
    assert model_a.CHECKPOINT_SHA256 == "dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e"
    assert model_a.EMBEDDING_DIM == 128 and model_a.REPRESENTATION_DIM == 512
    assert model_a.M4_POCKET_IDS == ("7TRQ_M4_allosteric", "7TRP_M4_allosteric", "7TRS_M4_allosteric")


def test_pinned_pdbs_match():
    from project.pacer_fkg_v02.frozen_runner.spec import sha256_file
    for pdb, expected in model_a.EXPECTED_PDB_SHA256.items():
        path = REPO / f"project/data/pdb/{pdb}.pdb"
        assert path.is_file(), f"missing {path}"
        assert sha256_file(path) == expected, f"{pdb} SHA256 mismatch"


def test_model_a_validate_rejects_missing_archive(tmp_path):
    _expect_raises(model_a.ModelAError, model_a.validate_archive, tmp_path / "absent.npz")


def test_projection_rejects_bad_shapes():
    head = {
        "linear1.weight": np.zeros((512, 512), dtype=np.float32),
        "linear1.bias": np.zeros((512,), dtype=np.float32),
        "linear2.weight": np.zeros((128, 512), dtype=np.float32),
        "linear2.bias": np.zeros((128,), dtype=np.float32),
    }
    _expect_raises(projection.ProjectionError, projection.forward, np.zeros((3, 256), dtype=np.float32), head)


def test_projection_forward_is_l2_normalised():
    rng = np.random.default_rng(3)
    head = {
        "linear1.weight": rng.standard_normal((512, 512), dtype=np.float32),
        "linear1.bias": np.zeros((512,), dtype=np.float32),
        "linear2.weight": rng.standard_normal((128, 512), dtype=np.float32),
        "linear2.bias": np.zeros((128,), dtype=np.float32),
    }
    out = projection.forward(rng.standard_normal((5, 512), dtype=np.float32), head)
    assert out.shape == (5, 128)
    assert np.allclose(np.linalg.norm(out, axis=1), 1.0, atol=1e-5)


def test_model_b_requires_three_deployed_seeds():
    _expect_raises(model_b.ModelBError, model_b.score,
                   REPO / "project/results/nonexistent.npz", [Path("a.pt")])


def test_dual_model_rule_intersection():
    records = [
        {"candidate_id": "A", "old_rank": 10, "new_rank": 5},
        {"candidate_id": "B", "old_rank": 10, "new_rank": 30},
        {"candidate_id": "C", "old_rank": 60, "new_rank": 5},
        {"candidate_id": "D", "old_rank": 50, "new_rank": 25},
    ]
    decision = dual_model.decide(records)
    assert decision["rule"] == "old_top50 INTERSECT new_top25"
    assert decision["dual_top_candidates"] == [
        {"candidate_id": "A", "old_rank": 10, "new_rank": 5},
        {"candidate_id": "D", "old_rank": 50, "new_rank": 25},
    ]
    assert decision["fusion_used"] is False


def test_dual_model_rejects_fused_columns():
    _expect_raises(dual_model.DualModelError, dual_model.assert_no_fusion, ["candidate_id", "fused_score"])


def test_committed_shortlist_reproduces_rule():
    assert SHORTLIST.is_file(), f"missing committed shortlist: {SHORTLIST}"
    report = dual_model.verify_committed_shortlist(SHORTLIST)
    assert report["n_rows"] >= 5
    assert report["rule_reproduced"], f"rule mismatches: {report['violations']}"


def _run_all() -> int:
    namespace = {k: v for k, v in sorted(globals().items()) if k.startswith("test_")}
    import tempfile, traceback
    failures = 0
    for name, fn in namespace.items():
        kwargs = {}
        if "tmp_path" in fn.__code__.co_varnames[: fn.__code__.co_argcount]:
            kwargs["tmp_path"] = Path(tempfile.mkdtemp(prefix="drugclip_test_"))
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
