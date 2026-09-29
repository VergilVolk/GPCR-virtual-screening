#!/usr/bin/env python3
"""Audit saved DrugCLIP LIT-PCBA outputs against the official scoring protocol."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit.ML.Scoring.Scoring import CalcAUC, CalcBEDROC, CalcEnrichment


FRACTIONS = [0.005, 0.01, 0.02, 0.05]


def official_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, object]:
    ranked = np.column_stack((scores, labels))
    ranked = ranked[np.argsort(ranked[:, 0], kind="mergesort")[::-1]]
    enrichments = CalcEnrichment(ranked, 1, FRACTIONS)
    return {
        "roc_auc": float(CalcAUC(ranked, 1)),
        "bedroc_alpha80_5": float(CalcBEDROC(ranked, 1, 80.5)),
        "enrichment": {str(fraction): float(value) for fraction, value in zip(FRACTIONS, enrichments)},
    }


def audit_target(npz_path: Path) -> dict[str, object]:
    metric_path = npz_path.with_suffix(".metrics.json")
    saved_metric = json.loads(metric_path.read_text(encoding="utf-8"))
    with np.load(npz_path, allow_pickle=False) as payload:
        molecule_ids = payload["molecule_ids"].astype(str)
        pocket_ids = payload["pocket_ids"].astype(str)
        molecule_embeddings = payload["molecule_embeddings"].astype(np.float64)
        pocket_embeddings = payload["pocket_embeddings"].astype(np.float64)
        labels = payload["labels"].astype(np.int64)
        saved_scores = payload["scores"].astype(np.float64)

    recomputed_scores = (pocket_embeddings @ molecule_embeddings.T).max(axis=0)
    official = official_metrics(labels, saved_scores)
    stored_ef = {str(fraction): float(saved_metric[f"ef{fraction:g}"]) for fraction in FRACTIONS}
    metric_deltas = {
        "roc_auc": abs(float(saved_metric["roc_auc"]) - float(official["roc_auc"])),
        "bedroc_alpha80_5": abs(
            float(saved_metric["bedroc_alpha80_5"]) - float(official["bedroc_alpha80_5"])
        ),
        "enrichment_max": max(
            abs(stored_ef[key] - official["enrichment"][key]) for key in official["enrichment"]
        ),
    }
    return {
        "target": saved_metric["target"],
        "n": int(labels.size),
        "actives": int(labels.sum()),
        "molecule_id_unique": int(np.unique(molecule_ids).size),
        "molecule_id_duplicates": int(labels.size - np.unique(molecule_ids).size),
        "pockets": int(pocket_embeddings.shape[0]),
        "pocket_id_unique": int(np.unique(pocket_ids).size),
        "finite": bool(
            np.isfinite(molecule_embeddings).all()
            and np.isfinite(pocket_embeddings).all()
            and np.isfinite(saved_scores).all()
        ),
        "molecule_norm_max_abs_error": float(
            np.max(np.abs(np.linalg.norm(molecule_embeddings, axis=1) - 1.0))
        ),
        "pocket_norm_max_abs_error": float(
            np.max(np.abs(np.linalg.norm(pocket_embeddings, axis=1) - 1.0))
        ),
        "score_recompute_max_abs_error": float(np.max(np.abs(saved_scores - recomputed_scores))),
        "official_metric_deltas": metric_deltas,
        "stored": {
            "roc_auc": float(saved_metric["roc_auc"]),
            "bedroc_alpha80_5": float(saved_metric["bedroc_alpha80_5"]),
            "ef0.01": float(saved_metric["ef0.01"]),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows = [audit_target(path) for path in sorted(args.output.glob("*.npz"))]
    report = {
        "n_targets": len(rows),
        "targets": rows,
        "pass": bool(rows)
        and all(
            row["finite"]
            and row["score_recompute_max_abs_error"] < 1e-6
            # Saved embeddings/scores are float32, so parity is assessed at a
            # float32-appropriate tolerance rather than bitwise equality.
            and max(row["official_metric_deltas"].values()) < 1e-6
            for row in rows
        ),
        "claim_boundary": (
            "This verifies implementation parity with the repository scoring protocol only; "
            "it does not establish paper-table parity, statistical significance, or SOTA."
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
