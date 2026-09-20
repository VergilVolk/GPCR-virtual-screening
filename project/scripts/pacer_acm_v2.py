"""Frozen PACER-ACM v2 calibration and refusal engine.

This module deliberately exposes a small, auditable algorithm. It must not be
modified after a fourth-campaign source is registered without incrementing the
method version and declaring the external test developmental.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


PROJECT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT / "config" / "pacer_acm_v2_frozen.json"
ALLOWED_ENDPOINTS = {"log_alpha_beta", "log_tauB"}


@dataclass(frozen=True)
class Prediction:
    status: str
    endpoint: str
    prediction: float | None
    lower: float | None
    upper: float | None
    max_anchor_similarity: float | None
    reason: str


def load_config(path: Path = CONFIG_PATH) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _fingerprints(smiles: Iterable[str]):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    result = []
    for item in smiles:
        molecule = Chem.MolFromSmiles(str(item))
        if molecule is None:
            raise ValueError(f"Invalid SMILES: {item}")
        result.append(generator.GetFingerprint(molecule))
    return result


def select_three_anchors(smiles: list[str]) -> list[int]:
    """Label-blind medoid + deterministic max-min selection."""
    if len(smiles) < 4:
        raise ValueError("At least four molecules are required: three anchors plus one query.")
    fps = _fingerprints(smiles)
    matrix = np.asarray([[DataStructs.TanimotoSimilarity(a, b) for b in fps] for a in fps])
    selected = [int(np.argmax(matrix.mean(axis=1)))]
    while len(selected) < 3:
        remaining = [i for i in range(len(smiles)) if i not in selected]
        novelty = [min(1.0 - matrix[i, j] for j in selected) for i in remaining]
        selected.append(remaining[int(np.argmax(novelty))])
    return selected


def calibrate_one(
    *, endpoint: str, query_smiles: str, query_base: float,
    anchor_smiles: list[str], anchor_base: list[float], anchor_observed: list[float],
    anchor_qc: list[bool], context_matches: bool, same_series: bool,
    config: dict | None = None,
) -> Prediction:
    cfg = config or load_config()
    if endpoint not in ALLOWED_ENDPOINTS:
        return Prediction("REFUSE", endpoint, None, None, None, None, "endpoint_not_supported")
    if not (len(anchor_smiles) == len(anchor_base) == len(anchor_observed) == len(anchor_qc) == 3):
        return Prediction("REFUSE", endpoint, None, None, None, None, "exactly_three_anchors_required")
    if not context_matches:
        return Prediction("REFUSE", endpoint, None, None, None, None, "assay_context_mismatch")
    if not same_series:
        return Prediction("REFUSE", endpoint, None, None, None, None, "outside_declared_series")
    if not all(anchor_qc):
        return Prediction("REFUSE", endpoint, None, None, None, None, "anchor_qc_failed")
    values = np.asarray([query_base, *anchor_base, *anchor_observed], dtype=float)
    if not np.isfinite(values).all():
        return Prediction("REFUSE", endpoint, None, None, None, None, "nonfinite_input")

    fps = _fingerprints([query_smiles, *anchor_smiles])
    max_similarity = max(DataStructs.TanimotoSimilarity(fps[0], fp) for fp in fps[1:])
    min_similarity = float(cfg["applicability_domain"]["minimum_query_to_anchor_tanimoto"])
    if max_similarity < min_similarity:
        return Prediction("REFUSE", endpoint, None, None, None, float(max_similarity), "outside_chemical_domain")

    residual = np.asarray(anchor_observed, float) - np.asarray(anchor_base, float)
    residual_sd = float(np.std(residual, ddof=1))
    max_sd = float(cfg["applicability_domain"]["maximum_anchor_residual_sd_log_units"])
    if residual_sd > max_sd:
        return Prediction("REFUSE", endpoint, None, None, None, float(max_similarity), "anchor_residuals_incoherent")

    prediction = float(query_base + residual.mean())
    # t(2, .975)=4.303; retain a floor because n=3 can otherwise look falsely precise.
    half_width = max(4.303 * residual_sd / np.sqrt(3),
                     float(cfg["uncertainty"]["minimum_interval_half_width_log_units"]))
    return Prediction("PREDICT", endpoint, prediction, prediction - half_width,
                      prediction + half_width, float(max_similarity), "within_frozen_domain")


def smoke_test() -> None:
    accepted = calibrate_one(
        endpoint="log_alpha_beta", query_smiles="CCOc1ccccc1", query_base=0.5,
        anchor_smiles=["COc1ccccc1", "CCOc1ccccc1C", "CCOc1ccc(F)cc1"],
        anchor_base=[0.2, 0.4, 0.6], anchor_observed=[0.4, 0.7, 0.8],
        anchor_qc=[True, True, True], context_matches=True, same_series=True,
    )
    refused = calibrate_one(
        endpoint="pKB", query_smiles="CCO", query_base=1.0,
        anchor_smiles=["CCO", "CCCO", "CCCCO"], anchor_base=[1, 1, 1],
        anchor_observed=[1, 1, 1], anchor_qc=[True] * 3,
        context_matches=True, same_series=True,
    )
    assert accepted.status == "PREDICT" and refused.status == "REFUSE"
    print(json.dumps({"accepted": accepted.__dict__, "refused": refused.__dict__}, indent=2))


if __name__ == "__main__":
    smoke_test()
