"""Post-hoc mechanism few-shot adaptation after cross-publication failure.

This is a developmental analysis, not an untouched external result. Three
label-blind, chemically diverse query-domain anchors are revealed jointly for
five operational-allostery endpoints. The remaining molecules are evaluated.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "mechanistic_chembl" / "mechanistic_multidomain.csv"
BASE = P / "results" / "pacer_mechanistic_multidomain_v01" / "cross_document_predictions.csv"
OUT = P / "results" / "pacer_mechanistic_fewshot_v01"
ENDPOINTS = ("pKB", "cAMP_log_tauB", "cAMP_log_alpha_beta",
             "arrestin_log_tauB", "arrestin_log_alpha_beta")
SEED = 20260830


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def fps(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in smiles]


def diverse_anchors(smiles, k=3):
    bitvectors = fps(smiles)
    similarity = np.asarray([[DataStructs.TanimotoSimilarity(left, right)
                              for right in bitvectors] for left in bitvectors])
    selected = [int(np.argmax(similarity.mean(axis=1)))]
    while len(selected) < k:
        remaining = [i for i in range(len(smiles)) if i not in selected]
        score = [min(1.0 - similarity[i, j] for j in selected) for i in remaining]
        selected.append(remaining[int(np.argmax(score))])
    return np.asarray(selected), similarity


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    molecules = pd.read_csv(DATA)
    base = pd.read_csv(BASE)
    rows, audit_folds = [], []
    for query_document in ("Jorg2023", "Liu2024"):
        block = base[(base.query_document == query_document) & base.endpoint.isin(ENDPOINTS)].copy()
        observed = block.pivot(index="molecule_chembl_id", columns="endpoint", values="observed")
        predicted = block.pivot(index="molecule_chembl_id", columns="endpoint", values="Pooled-Ridge")
        complete = observed.dropna(subset=list(ENDPOINTS)).index.intersection(predicted.dropna(subset=list(ENDPOINTS)).index)
        metadata = molecules[molecules.molecule_chembl_id.isin(complete)].set_index("molecule_chembl_id").loc[complete]
        anchors_position, similarity = diverse_anchors(metadata.canonical_smiles.tolist(), 3)
        anchor_ids = metadata.index.to_numpy()[anchors_position]
        query_ids = metadata.index.difference(anchor_ids)
        residual = observed.loc[anchor_ids, ENDPOINTS].to_numpy(float) - predicted.loc[anchor_ids, ENDPOINTS].to_numpy(float)
        center = residual.mean(axis=0)
        u, singular, vt = np.linalg.svd(residual - center, full_matrices=False)
        rank1 = center + np.outer(u[:, 0] * singular[0], vt[0])
        id_to_position = {item: i for i, item in enumerate(metadata.index)}
        for molecule_id in query_ids:
            qpos = id_to_position[molecule_id]
            sims = similarity[qpos, anchors_position]
            weights = np.maximum(sims, 1e-4) ** 3
            weights /= weights.sum()
            for endpoint_index, endpoint in enumerate(ENDPOINTS):
                baseline = float(predicted.loc[molecule_id, endpoint])
                endpoint_residual = float(weights @ residual[:, endpoint_index])
                mechanism_residual = float(weights @ rank1[:, endpoint_index])
                rows.append({
                    "query_document": query_document, "molecule_chembl_id": molecule_id,
                    "endpoint": endpoint, "observed": float(observed.loc[molecule_id, endpoint]),
                    "CrossDoc-Pooled": baseline,
                    "ThreeShot-Offset": baseline + float(center[endpoint_index]),
                    "ThreeShot-KernelResidual": baseline + endpoint_residual,
                    "PACER-MechanismFS": baseline + mechanism_residual,
                })
        audit_folds.append({"query_document": query_document, "anchor_ids": anchor_ids.tolist(),
                            "n_complete_before_anchors": int(len(complete)),
                            "n_query_after_anchors": int(len(query_ids)),
                            "anchor_selection": "label-blind ECFP4 max-min diversity",
                            "rank1_explained_fraction": float(singular[0] ** 2 / np.sum(singular ** 2))})

    predictions = pd.DataFrame(rows)
    methods = ("CrossDoc-Pooled", "ThreeShot-Offset", "ThreeShot-KernelResidual", "PACER-MechanismFS")
    metrics = []
    for endpoint, block in predictions.groupby("endpoint"):
        for method in methods:
            directions = [rho(part.observed, part[method]) for _, part in block.groupby("query_document") if len(part) >= 3]
            metrics.append({"endpoint": endpoint, "method": method, "n": len(block),
                            "Spearman": rho(block.observed, block[method]),
                            "MAE": float(mean_absolute_error(block.observed, block[method])),
                            "worst_direction_Spearman": min(directions) if directions else np.nan})
    metrics = pd.DataFrame(metrics)
    predictions.to_csv(OUT / "predictions.csv", index=False)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    audit = {"status": "post_hoc_mechanistic_fewshot_development",
             "anchors_per_document": 3, "folds": audit_folds,
             "query_labels_used": "only the three declared anchors per document",
             "claim_boundary": "Developmental retrospective assay calibration; not independent external validation or PAM confirmation."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER-M4 mechanism few-shot adaptation", "",
              "Three label-blind diversity anchors per query publication calibrate five operational-allostery endpoints.",
              "This was developed after observing the zero-shot cross-publication failure and is therefore post-hoc.", "",
              "| Endpoint | Method | N | Spearman | MAE | Worst direction |", "|---|---|---:|---:|---:|---:|"]
    for row in metrics.sort_values(["endpoint", "Spearman"], ascending=[True, False]).itertuples():
        report.append(f"| {row.endpoint} | {row.method} | {int(row.n)} | {row.Spearman:.3f} | {row.MAE:.3f} | {row.worst_direction_Spearman:.3f} |")
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
