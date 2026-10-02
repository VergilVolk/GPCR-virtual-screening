from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "select_pareto_candidates.py"
SPEC = importlib.util.spec_from_file_location("pareto_selector", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def test_pareto_mask_keeps_only_non_dominated_rows() -> None:
    values = np.asarray([[1.0, 1.0], [0.5, 0.5], [1.1, 0.8], [0.8, 1.1]])
    np.testing.assert_array_equal(MOD.pareto_mask(values), [True, False, True, True])


def test_assign_pareto_fronts_is_complete() -> None:
    frame = pd.DataFrame({"a": [1.0, 0.5, 0.1], "b": [1.0, 0.5, 0.1]})
    result = MOD.assign_pareto_fronts(frame, ["a", "b"])
    assert result.pareto_front.tolist() == [0, 1, 2]


def test_choose_diverse_returns_unique_candidates() -> None:
    frame = pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c"],
            "canonical_smiles": ["c1ccccc1", "Cc1ccccc1", "CCO"],
            "pareto_front": [0, 0, 1],
            "pocket_residue_coverage": [0.9, 0.8, 0.7],
            "QED": [0.7, 0.8, 0.6],
            "SA_score": [2.0, 2.1, 1.5],
        }
    )
    result = MOD.choose_diverse(frame, 2)
    assert len(result) == 2
    assert result.candidate_id.is_unique
    assert "diversity_cluster" in result
