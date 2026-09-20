"""Endpoint-complete and random-anchor robustness for assay calibration."""

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
OUT = P / "results" / "pacer_assay_anchor_robustness_v01"
ENDPOINTS = ("pKB", "cAMP_log_tauB", "cAMP_log_alpha_beta",
             "arrestin_log_tauB", "arrestin_log_alpha_beta")
SEED = 20260830
N_RANDOM = 1000


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def similarity_matrix(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in smiles]
    return np.asarray([[DataStructs.TanimotoSimilarity(left, right) for right in fps] for left in fps])


def diverse(S, k):
    selected = [int(np.argmax(S.mean(axis=1)))]
    while len(selected) < k:
        remaining = [i for i in range(len(S)) if i not in selected]
        scores = [min(1 - S[i, j] for j in selected) for i in remaining]
        selected.append(remaining[int(np.argmax(scores))])
    return np.asarray(selected)


def adapt(block, anchor_positions, S):
    query_positions = np.asarray([i for i in range(len(block)) if i not in set(anchor_positions)])
    residual = block.observed.to_numpy()[anchor_positions] - block["Pooled-Ridge"].to_numpy()[anchor_positions]
    offset = float(np.mean(residual))
    base = block["Pooled-Ridge"].to_numpy()[query_positions]
    kernel = []
    for position in query_positions:
        weights = np.maximum(S[position, anchor_positions], 1e-4) ** 3
        kernel.append(base[len(kernel)] + float(weights @ residual / weights.sum()))
    return query_positions, base, base + offset, np.asarray(kernel)


def summarize(predictions):
    methods = ("CrossDoc-Pooled", "ThreeShot-Offset", "ThreeShot-KernelResidual")
    rows = []
    for endpoint, block in predictions.groupby("endpoint"):
        for method in methods:
            directions = [rho(part.observed, part[method]) for _, part in block.groupby("query_document") if len(part) >= 3]
            rows.append({"endpoint": endpoint, "method": method, "n": len(block),
                         "Spearman": rho(block.observed, block[method]),
                         "MAE": float(mean_absolute_error(block.observed, block[method])),
                         "worst_direction_Spearman": min(directions) if directions else np.nan})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    molecules = pd.read_csv(DATA).set_index("molecule_chembl_id")
    base = pd.read_csv(BASE)
    deterministic_rows, selections = [], []
    random_records = []
    rng = np.random.default_rng(SEED)
    for endpoint in ENDPOINTS:
        direction_blocks = {}
        for document in ("Jorg2023", "Liu2024"):
            block = base[(base.query_document == document) & (base.endpoint == endpoint)].dropna(
                subset=["observed", "Pooled-Ridge"]).reset_index(drop=True)
            smiles = molecules.loc[block.molecule_chembl_id, "canonical_smiles"].tolist()
            S = similarity_matrix(smiles)
            anchors = diverse(S, min(3, len(block) - 1))
            query, raw, offset, kernel = adapt(block, anchors, S)
            for position, b, o, k in zip(query, raw, offset, kernel):
                deterministic_rows.append({"endpoint": endpoint, "query_document": document,
                                           "molecule_chembl_id": block.loc[position, "molecule_chembl_id"],
                                           "observed": block.loc[position, "observed"],
                                           "CrossDoc-Pooled": b, "ThreeShot-Offset": o,
                                           "ThreeShot-KernelResidual": k})
            selections.append({"endpoint": endpoint, "query_document": document,
                               "anchor_ids": block.loc[anchors, "molecule_chembl_id"].tolist(),
                               "n_available": len(block), "n_query": len(query)})
            direction_blocks[document] = (block, S)

        for repeat in range(N_RANDOM):
            collected = []
            for document, (block, S) in direction_blocks.items():
                anchors = rng.choice(len(block), min(3, len(block) - 1), replace=False)
                query, raw, offset, kernel = adapt(block, anchors, S)
                for position, b, o, k in zip(query, raw, offset, kernel):
                    collected.append({"document": document, "observed": block.loc[position, "observed"],
                                      "CrossDoc-Pooled": b, "ThreeShot-Offset": o,
                                      "ThreeShot-KernelResidual": k})
            sample = pd.DataFrame(collected)
            for method in ("ThreeShot-Offset", "ThreeShot-KernelResidual"):
                directions = [rho(part.observed, part[method]) for _, part in sample.groupby("document") if len(part) >= 3]
                random_records.append({"endpoint": endpoint, "repeat": repeat, "method": method,
                                       "Spearman": rho(sample.observed, sample[method]),
                                       "MAE": float(mean_absolute_error(sample.observed, sample[method])),
                                       "worst_direction_Spearman": min(directions) if directions else np.nan})

    predictions = pd.DataFrame(deterministic_rows)
    metrics = summarize(predictions)
    random = pd.DataFrame(random_records)
    random_summary = random.groupby(["endpoint", "method"]).agg(
        median_Spearman=("Spearman", "median"), q025_Spearman=("Spearman", lambda x: x.quantile(.025)),
        q975_Spearman=("Spearman", lambda x: x.quantile(.975)),
        median_MAE=("MAE", "median"), median_worst=("worst_direction_Spearman", "median")
    ).reset_index()
    predictions.to_csv(OUT / "diverse_anchor_predictions.csv", index=False)
    metrics.to_csv(OUT / "diverse_anchor_metrics.csv", index=False)
    random.to_csv(OUT / "random_anchor_repeats.csv", index=False)
    random_summary.to_csv(OUT / "random_anchor_summary.csv", index=False)
    audit = {"status": "post_hoc_assay_anchor_robustness", "random_repeats": N_RANDOM,
             "selections": selections,
             "claim_boundary": "Retrospective anchor sensitivity analysis; no independent SOTA or PAM confirmation."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metrics.to_string(index=False)); print(random_summary.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
