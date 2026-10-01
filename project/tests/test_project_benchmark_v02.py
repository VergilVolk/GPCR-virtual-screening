from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit.ML.Scoring.Scoring import CalcBEDROC


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/recompute_binding_benchmark_v02.py"
SPEC = importlib.util.spec_from_file_location("binding_v02", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def test_bedroc_and_ef_match_reference_definition():
    y = np.array([1, 0, 1, 0, 0, 1, 0, 0, 0, 0])
    score = np.array([.9, .8, .7, .6, .5, .4, .3, .2, .1, 0.0])
    got = MOD.metric(y, score)
    order = np.argsort(-score, kind="mergesort")
    ranked = [[float(score[i]), int(y[i])] for i in order]
    assert got["bedroc_alpha20"] == pytest.approx(CalcBEDROC(ranked, 1, 20.0))
    assert got["ef1pct"] == pytest.approx(y[order[:1]].mean() / y.mean())
    assert got["ef5pct"] == pytest.approx(y[order[:1]].mean() / y.mean())


def test_ensemble_rejects_row_order_mismatch(tmp_path):
    base = pd.DataFrame({"target": ["A", "A"], "label": [1, 0],
                         "canonical_smiles": ["C", "CC"], "murcko_scaffold": ["", "C"],
                         "score": [.8, .2]})
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    base.to_csv(a, index=False)
    base.iloc[::-1].to_csv(b, index=False)
    with pytest.raises(ValueError, match="row identity/order mismatch"):
        MOD.ensemble([a, b], "score")
