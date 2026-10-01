#!/usr/bin/env python3
"""Frozen evaluation of DrugCLIP variants on heterogeneous external M4 data.

Assays are never pooled. Binary PAM calls, continuous potency, and patent
ordinal bins are scored with endpoint-appropriate metrics. Exact M4 training
overlaps are excluded before any metric is computed.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata, spearmanr
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F


class Projection(nn.Module):
    def __init__(self, state: dict):
        super().__init__()
        self.linear1 = nn.Linear(state["linear1.weight"].shape[1], state["linear1.weight"].shape[0])
        self.linear2 = nn.Linear(state["linear2.weight"].shape[1], state["linear2.weight"].shape[0])
        self.load_state_dict(state)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return F.normalize(self.linear2(F.relu(self.linear1(x))), dim=-1)


def retrieval(mol_rep: np.ndarray, pocket_rep: np.ndarray, state: dict) -> np.ndarray:
    mol = Projection(copy.deepcopy(state["mol_project"])).eval()
    pocket = Projection(copy.deepcopy(state["pocket_project"])).eval()
    with torch.inference_mode():
        zm = mol(torch.as_tensor(mol_rep, dtype=torch.float32))
        zp = pocket(torch.as_tensor(pocket_rep, dtype=torch.float32))
        state_scores = zm @ zp.T
        return (0.1 * torch.logsumexp(state_scores / 0.1, dim=1)).numpy()


def fingerprints(smiles: list[str]) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    output = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, value in enumerate(smiles):
        fp = generator.GetFingerprint(Chem.MolFromSmiles(value))
        DataStructs.ConvertToNumpyArray(fp, output[i])
    return output


def ci_bootstrap(y: np.ndarray, score: np.ndarray, metric: str, seed: int, n: int) -> list[float]:
    rng = np.random.default_rng(seed); values = []
    for _ in range(n):
        index = rng.integers(0, len(y), len(y)); yy, ss = y[index], score[index]
        if metric == "auc":
            if np.unique(yy).size < 2:
                continue
            value = roc_auc_score(yy, ss)
        else:
            value = spearmanr(yy, ss).statistic
            if not np.isfinite(value):
                continue
        values.append(value)
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975]))) if values else [float("nan")] * 3


def paired_delta_ci(y: np.ndarray, a: np.ndarray, b: np.ndarray, metric: str,
                    seed: int, n: int) -> list[float]:
    rng = np.random.default_rng(seed); values = []
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    for _ in range(n):
        if metric == "auc":
            index = np.concatenate([rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)])
            values.append(roc_auc_score(y[index], a[index]) - roc_auc_score(y[index], b[index]))
        else:
            index = rng.integers(0, len(y), len(y))
            da = spearmanr(y[index], a[index]).statistic
            db = spearmanr(y[index], b[index]).statistic
            if np.isfinite(da) and np.isfinite(db):
                values.append(da - db)
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975])))


def metric_block(frame: pd.DataFrame, methods: list[str], bootstrap: int) -> dict:
    endpoint = frame.endpoint.iloc[0]
    binary = frame.label_binary.notna().all()
    y_column = "label_binary" if binary else ("potency_value" if frame.potency_value.notna().all() else "ordinal_class")
    y = frame[y_column].to_numpy(float)
    result = {"endpoint": endpoint, "n": int(len(frame)), "outcome": y_column, "methods": {}}
    if binary:
        result.update({"positive": int(y.sum()), "negative": int(len(y) - y.sum())})
    for j, method in enumerate(methods):
        score = frame[method].to_numpy(float)
        if binary:
            result["methods"][method] = {
                "roc_auc": float(roc_auc_score(y, score)),
                "average_precision": float(average_precision_score(y, score)),
                "roc_auc_bootstrap_95ci": ci_bootstrap(y, score, "auc", 20260924 + j, bootstrap),
            }
        else:
            rho = spearmanr(y, score).statistic
            result["methods"][method] = {
                "spearman": float(rho),
                "spearman_bootstrap_95ci": ci_bootstrap(y, score, "spearman", 20260924 + j, bootstrap),
            }
    comparison_metric = "auc" if binary else "spearman"
    result["paired_differences"] = {
        "gpcr_triplet_minus_official_95ci": paired_delta_ci(
            y, frame.drugclip_gpcr_triplet.to_numpy(float), frame.drugclip_official.to_numpy(float),
            comparison_metric, 20261001, bootstrap),
        "m4_matched_minus_gpcr_triplet_95ci": paired_delta_ci(
            y, frame.drugclip_m4_matched.to_numpy(float), frame.drugclip_gpcr_triplet.to_numpy(float),
            comparison_metric, 20261002, bootstrap),
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--training-classification", type=Path, required=True)
    parser.add_argument("--training-potency", type=Path, required=True)
    parser.add_argument("--official-projection", type=Path, required=True)
    parser.add_argument("--gpcr-projection", type=Path, required=True)
    parser.add_argument("--matched-projections", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=3000)
    args = parser.parse_args()

    archive = np.load(args.representations, allow_pickle=False)
    molecule_ids = list(map(str, archive["molecule_ids"])); index = {v: i for i, v in enumerate(molecule_ids)}
    manifest = pd.read_csv(args.manifest)
    # DrugCLIP's retrieval dataset exposes the canonical SMILES as ``smi_name``
    # even when the LMDB also carries a human-readable label.
    keys = [str(mid) if str(mid) in index else str(smi)
            for mid, smi in zip(manifest.external_molecule_id, manifest.canonical_smiles)]
    order = np.asarray([index[v] for v in keys])
    molecule_rep = archive["molecule_representations"].astype(np.float32)[order]
    pocket_rep = archive["pocket_representations"].astype(np.float32)
    official = torch.load(args.official_projection, map_location="cpu")
    gpcr = torch.load(args.gpcr_projection, map_location="cpu")
    matched = [torch.load(path, map_location="cpu") for path in args.matched_projections]
    manifest["drugclip_official"] = retrieval(molecule_rep, pocket_rep, official)
    manifest["drugclip_gpcr_triplet"] = retrieval(molecule_rep, pocket_rep, gpcr)
    manifest["drugclip_m4_matched"] = np.mean(
        [retrieval(molecule_rep, pocket_rep, state) for state in matched], axis=0)

    classification = pd.read_csv(args.training_classification)
    potency = pd.read_csv(args.training_potency)
    train_class_x = fingerprints(classification.canonical_smiles.astype(str).tolist())
    train_pot_x = fingerprints(potency.canonical_smiles.astype(str).tolist())
    external_x = fingerprints(manifest.canonical_smiles.astype(str).tolist())
    classifier = RandomForestClassifier(
        n_estimators=600, class_weight="balanced_subsample", min_samples_leaf=2,
        random_state=42, n_jobs=-1).fit(train_class_x, classification.target.to_numpy(int))
    regressor = RandomForestRegressor(
        n_estimators=600, min_samples_leaf=2, random_state=42, n_jobs=-1
    ).fit(train_pot_x, potency.pEC50.to_numpy(float))
    manifest["ecfp_rf_pam"] = classifier.predict_proba(external_x)[:, 1]
    manifest["ecfp_rf_potency"] = regressor.predict(external_x)
    # Fixed 50:50 rank fusion, transferred unchanged from the retrospective M4
    # experiment. It is not fitted or reweighted on any external endpoint.
    manifest["drugclip_2d_rank_fusion"] = (
        rankdata(manifest.drugclip_m4_matched) + rankdata(manifest.ecfp_rf_pam)
    ) / (2.0 * len(manifest))

    eligible = manifest[manifest.eligible.astype(bool) & ~manifest.exact_training_overlap_computed.astype(bool)].copy()
    # Identical molecules within one endpoint are one statistical unit.
    score_columns = ["drugclip_official", "drugclip_gpcr_triplet", "drugclip_m4_matched",
                     "ecfp_rf_pam", "ecfp_rf_potency", "drugclip_2d_rank_fusion"]
    value_columns = ["label_binary", "potency_value", "ordinal_class", *score_columns]
    grouped = eligible.groupby(["dataset", "endpoint", "external_molecule_id"], as_index=False)[value_columns].mean()
    reports = {}
    for dataset, frame in grouped.groupby("dataset", sort=False):
        methods = ["drugclip_official", "drugclip_gpcr_triplet", "drugclip_m4_matched"]
        if frame.label_binary.notna().all():
            methods += ["ecfp_rf_pam", "drugclip_2d_rank_fusion"]
        else:
            methods += ["ecfp_rf_potency"]
        reports[dataset] = metric_block(frame.reset_index(drop=True), methods, args.bootstrap)
    report = {
        "method_freeze": "No external label was used for checkpoint, seed, weight, or hyperparameter selection.",
        "n_endpoint_rows": int(len(manifest)),
        "n_zero_shot_rows": int(len(eligible)),
        "n_unique_zero_shot_units": int(len(grouped)),
        "matched_checkpoints": list(map(str, args.matched_projections)),
        "datasets": reports,
        "claim_boundary": (
            "Cross-source retrospective stress test. Some sources were used in earlier PACER analyses; "
            "this is not a prospectively untouched or wet-lab PAM validation."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    manifest.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
