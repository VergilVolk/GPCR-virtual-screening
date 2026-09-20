"""Strict molecule-held-out development test of a context-conditioned kernel."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "pacer_external_us20260055116_v01" / "external_unique_active_moieties.csv"
OUT = PROJECT / "results" / "pacer_context_kernel_v01"
ORDER = {"A": 3.0, "B": 2.0, "C": 1.0, "D": 0.0}
ENDPOINTS = {"human_m4_perk": ("human", "pERK"),
             "rat_m4_perk": ("rat", "pERK"),
             "human_m4_gtpgs": ("human", "GTPgammaS")}
ALPHAS = [0.01, 0.1, 1.0, 10.0]
CONTEXTS = ["shared", "same_endpoint", "shared_species_readout", "weak_shared_species_readout"]
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def concordance(y, p):
    y, p = np.asarray(y), np.asarray(p)
    credit = []
    for i in range(len(y)):
        for j in range(i + 1, len(y)):
            if y[i] == y[j]:
                continue
            direction = np.sign(y[i] - y[j]) * np.sign(p[i] - p[j])
            credit.append(1.0 if direction > 0 else .5 if direction == 0 else 0.0)
    return float(np.mean(credit)) if credit else np.nan


def context_similarity(a, b, kind):
    same_endpoint = a["endpoint"] == b["endpoint"]
    same_species = a["species"] == b["species"]
    same_readout = a["readout"] == b["readout"]
    if kind == "shared":
        return 1.0
    if kind == "same_endpoint":
        return 1.0 if same_endpoint else 0.0
    if kind == "shared_species_readout":
        return 1.0 + float(same_species) + float(same_readout)
    if kind == "weak_shared_species_readout":
        return .25 + float(same_species) + float(same_readout)
    raise ValueError(kind)


def kernel(a: pd.DataFrame, b: pd.DataFrame, fps, kind: str):
    matrix = np.zeros((len(a), len(b)), float)
    for i, row_a in enumerate(a.to_dict("records")):
        for j, row_b in enumerate(b.to_dict("records")):
            chem = DataStructs.TanimotoSimilarity(fps[row_a["molecule_index"]],
                                                   fps[row_b["molecule_index"]])
            matrix[i, j] = chem * context_similarity(row_a, row_b, kind)
    return matrix


def predict(train, query, kernels, kind, alpha):
    train_ids = train.row_id.to_numpy(int)
    query_ids = query.row_id.to_numpy(int)
    full = kernels[kind]
    K = full[np.ix_(train_ids, train_ids)]
    Kq = full[np.ix_(query_ids, train_ids)]
    y = train.y.to_numpy(float)
    mean = y.mean()
    coef = np.linalg.solve(K + alpha * np.eye(len(K)), y - mean)
    return mean + Kq @ coef


def inner_select(train, kernels, allowed_kinds):
    best = None
    for kind in allowed_kinds:
        for alpha in ALPHAS:
            truth, prediction = [], []
            for held in train.compound_id.unique():
                fit = train[train.compound_id != held]
                query = train[train.compound_id == held]
                if not len(fit) or not len(query):
                    continue
                truth.extend(query.y)
                prediction.extend(predict(fit, query, kernels, kind, alpha))
            score = concordance(np.asarray(truth), np.asarray(prediction))
            candidate = (score, -alpha, kind, alpha)
            if best is None or candidate > best:
                best = candidate
    return best[2], best[3], best[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    molecules = pd.read_csv(SOURCE).reset_index(drop=True)
    fps = [FPGEN.GetFingerprint(Chem.MolFromSmiles(x)) for x in molecules.canonical_smiles]
    rows = []
    for i, mol in molecules.iterrows():
        for endpoint, (species, readout) in ENDPOINTS.items():
            label = mol[endpoint]
            if label in ORDER:
                rows.append({"molecule_index": i, "compound_id": mol.compound_id,
                             "endpoint": endpoint, "species": species, "readout": readout,
                             "label": label, "y": ORDER[label]})
    data = pd.DataFrame(rows)
    data["row_id"] = np.arange(len(data))
    kernels = {kind: kernel(data, data, fps, kind) for kind in CONTEXTS}
    predictions, selections = [], []
    for held in data.compound_id.unique():
        train = data[data.compound_id != held].reset_index(drop=True)
        query = data[data.compound_id == held].reset_index(drop=True)
        for method, kinds in {"EndpointOnly-KRR": ["same_endpoint"],
                              "PooledNoContext-KRR": ["shared"],
                              "PACER-ContextKRR": CONTEXTS}.items():
            kind, alpha, inner = inner_select(train, kernels, kinds)
            pred = predict(train, query, kernels, kind, alpha)
            selections.append({"held_compound": held, "method": method,
                               "selected_context_kernel": kind, "selected_alpha": alpha,
                               "inner_ordinal_concordance": inner})
            for (_, row), value in zip(query.iterrows(), pred):
                predictions.append({**row.to_dict(), "method": method,
                                    "prediction": float(value)})
    pred = pd.DataFrame(predictions)
    pred.to_csv(OUT / "outer_predictions.csv", index=False)
    pd.DataFrame(selections).to_csv(OUT / "inner_selections.csv", index=False)

    metric_rows = []
    for method, block in pred.groupby("method"):
        endpoint_scores = []
        for endpoint, part in block.groupby("endpoint"):
            c = concordance(part.y, part.prediction)
            endpoint_scores.append(c)
            metric_rows.append({"method": method, "scope": endpoint, "n": len(part),
                                "ordinal_concordance": c,
                                "spearman": float(spearmanr(part.y, part.prediction).statistic)})
        metric_rows.append({"method": method, "scope": "macro", "n": len(block),
                            "ordinal_concordance": float(np.mean(endpoint_scores)),
                            "spearman": np.nan})
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    macro = metrics[metrics.scope == "macro"].set_index("method").ordinal_concordance
    context_parts = metrics[(metrics.method == "PACER-ContextKRR") & (metrics.scope != "macro")]
    passed = bool(macro["PACER-ContextKRR"] > macro["EndpointOnly-KRR"] and
                  macro["PACER-ContextKRR"] > macro["PooledNoContext-KRR"] and
                  context_parts.ordinal_concordance.min() >= .50)
    audit = {"outer_unit": "entire_molecule_all_endpoints",
             "outer_endpoint_rows": int(len(data)), "outer_molecules": int(data.compound_id.nunique()),
             "all_hyperparameters_selected_inside_outer_training": True,
             "primary_macro_ordinal_concordance": macro.to_dict(),
             "development_pass": passed,
             "claim_boundary": "Viewed-source development; not external validation or SOTA."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER-ContextKernel development result", "",
              metrics.to_markdown(index=False, floatfmt=".3f"), "",
              f"Frozen development pass: **{passed}**.", "", audit["claim_boundary"]]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
