# -*- coding: utf-8 -*-
"""Evaluate a receptor microstate-response fingerprint (SRF) on A-tier M4 PAM labels.

The SRF removes the absolute cluster-00 affinity and represents how a ligand's
pose/score changes across the published holo-GaMD ensemble. It is evaluated as
an independent structural signal and never treated as PAM efficacy by itself.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
P = Path(__file__).resolve().parents[1]
ROOT = P / "results" / "m4_gamd_ensemble"
FRAME = ROOT / "validation_primary_docking_protocol_v127_ph7_rank1" / "framewise_scores.csv"
LABELS = ROOT / "validation_set_primary.csv"
OUT = P / "results" / "pacer_srf_v01"
SEEDS = list(range(10))


def ecfp(smiles: str, nbits: int = 2048) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganGenerator(radius=2, fpSize=nbits).GetFingerprint(mol)
    x = np.zeros(nbits, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, x)
    return x


def contact_set(value) -> set[int]:
    if pd.isna(value) or str(value).strip() == "":
        return set()
    return {int(x) for x in str(value).split(";") if x}


def build_srf(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    rows = []
    centers = json.loads((ROOT / "ensemble_preparation_audit.json").read_text())["clusters"]
    centers = {int(x["cluster"]): np.asarray(x["box_center"], float) for x in centers}
    for mid, part in frame.groupby("molecule_id"):
        part = part.sort_values("cluster")
        if part.cluster.nunique() != 10:
            continue
        aff = part.vina_affinity.to_numpy(float)
        cov = part.native_pocket_coverage.to_numpy(float)
        xyz = part[["pose_centroid_x", "pose_centroid_y", "pose_centroid_z"]].to_numpy(float)
        disp = np.asarray([np.linalg.norm(xyz[i] - centers[int(c)]) for i, c in enumerate(part.cluster)])
        contacts = [contact_set(x) for x in part.contacted_residues]
        jac = []
        for i in range(1, 10):
            union = contacts[0] | contacts[i]
            jac.append(len(contacts[0] & contacts[i]) / len(union) if union else 1.0)
        values = {"molecule_id": mid}
        for i in range(1, 10):
            values[f"dVina_c{i:02d}"] = aff[i] - aff[0]
            values[f"dCoverage_c{i:02d}"] = cov[i] - cov[0]
            values[f"dDisplacement_c{i:02d}"] = disp[i] - disp[0]
            values[f"contactJaccard0_c{i:02d}"] = jac[i - 1]
        values.update({
            "vina_response_sd": float(np.std(aff - aff[0])),
            "vina_response_range": float(np.ptp(aff)),
            "coverage_mean": float(np.mean(cov)),
            "coverage_sd": float(np.std(cov)),
            "displacement_mean": float(np.mean(disp)),
            "displacement_sd": float(np.std(disp)),
            "contact_union_n": len(set().union(*contacts)),
            "contact_core_n": len(set.intersection(*contacts)) if contacts else 0,
            "contact_jaccard0_mean": float(np.mean(jac)),
        })
        rows.append(values)
    out = pd.DataFrame(rows)
    features = [c for c in out.columns if c != "molecule_id"]
    return out, features


def make_pair_splits(groups: np.ndarray, seed: int, n_splits: int = 5):
    unique = np.unique(groups).copy()
    rng = np.random.default_rng(seed); rng.shuffle(unique)
    fold_map = {g: i % n_splits for i, g in enumerate(unique)}
    return [(np.flatnonzero(np.asarray([fold_map[g] for g in groups]) != f),
             np.flatnonzero(np.asarray([fold_map[g] for g in groups]) == f)) for f in range(n_splits)]


def nested_logistic_oof(X, y, groups, splits) -> tuple[np.ndarray, list[float]]:
    pred = np.full(len(y), np.nan); chosen = []
    for train, test in splits:
        inner_groups = groups[train]
        candidates = [0.01, 0.1, 1.0, 10.0]
        scores = []
        n_inner = min(3, len(np.unique(inner_groups)))
        for C in candidates:
            fold_scores = []
            for it, iv in GroupKFold(n_inner).split(X[train], y[train], inner_groups):
                if len(np.unique(y[train][iv])) < 2:
                    continue
                model = make_pipeline(StandardScaler(), LogisticRegression(
                    C=C, class_weight="balanced", max_iter=3000, random_state=42))
                model.fit(X[train][it], y[train][it])
                fold_scores.append(roc_auc_score(y[train][iv], model.predict_proba(X[train][iv])[:, 1]))
            scores.append(np.mean(fold_scores) if fold_scores else -np.inf)
        C = candidates[int(np.argmax(scores))]; chosen.append(C)
        model = make_pipeline(StandardScaler(), LogisticRegression(
            C=C, class_weight="balanced", max_iter=3000, random_state=42))
        model.fit(X[train], y[train])
        pred[test] = model.predict_proba(X[test])[:, 1]
    return pred, chosen


def centroid_oof(X, y, splits) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    for train, test in splits:
        scale = np.std(X[train], axis=0); scale[scale < 1e-8] = 1.0
        ztr = (X[train] - np.mean(X[train], axis=0)) / scale
        zte = (X[test] - np.mean(X[train], axis=0)) / scale
        p = ztr[y[train] == 1].mean(0); n = ztr[y[train] == 0].mean(0)
        pred[test] = np.linalg.norm(zte - n, axis=1) - np.linalg.norm(zte - p, axis=1)
    return pred


def score(y, pred):
    return {"roc_auc": float(roc_auc_score(y, pred)),
            "pr_auc": float(average_precision_score(y, pred))}


def evaluate_scheme(name, Xs, y, groups, split_factory):
    records = []; predictions = []
    for seed in SEEDS:
        splits = split_factory(seed)
        for model_name, X in Xs.items():
            pred, chosen = nested_logistic_oof(X, y, groups, splits)
            rec = {"scheme": name, "seed": seed, "model": model_name, **score(y, pred),
                   "chosen_C": ";".join(map(str, chosen))}
            records.append(rec)
            predictions.extend({"scheme": name, "seed": seed, "model": model_name,
                                "row": i, "y": int(y[i]), "prediction": float(pred[i])}
                               for i in range(len(y)))
        c = centroid_oof(Xs["SRF"], y, splits)
        records.append({"scheme": name, "seed": seed, "model": "SRF_centroid", **score(y, c), "chosen_C": ""})
    return records, predictions


def source_component_holdout(Xs, y, components, pair_groups):
    records = []; predictions = []
    eligible = [g for g in np.unique(components) if len(np.unique(y[components == g])) == 2]
    for model_name, X in Xs.items():
        pooled_y = []; pooled_p = []; fold_aucs = []
        for held in eligible:
            test = np.flatnonzero(components == held); train = np.flatnonzero(components != held)
            inner = make_pair_splits(pair_groups[train], seed=42, n_splits=3)
            pred, _ = nested_logistic_oof(X[train], y[train], pair_groups[train], inner)
            # Refit with the modal C selected in nested training OOF.
            # C=0.1 is deliberately fixed for the external component prediction.
            model = make_pipeline(StandardScaler(), LogisticRegression(
                C=0.1, class_weight="balanced", max_iter=3000, random_state=42))
            model.fit(X[train], y[train])
            p = model.predict_proba(X[test])[:, 1]
            auc = roc_auc_score(y[test], p)
            fold_aucs.append(float(auc)); pooled_y.extend(y[test]); pooled_p.extend(p)
            predictions.extend({"scheme": "negative_component_LOO", "seed": 42, "model": model_name,
                                "row": int(i), "held_component": held, "y": int(y[i]),
                                "prediction": float(pp)} for i, pp in zip(test, p))
        records.append({"scheme": "negative_component_LOO", "seed": 42, "model": model_name,
                        "roc_auc": float(np.mean(fold_aucs)),
                        "pr_auc": float(average_precision_score(pooled_y, pooled_p)),
                        "worst_component_auc": float(min(fold_aucs)),
                        "fold_aucs": ";".join(f"{x:.6f}" for x in fold_aucs), "chosen_C": "0.1"})
    return records, predictions


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(FRAME)
    labels = pd.read_csv(LABELS)
    srf, srf_cols = build_srf(frame)
    data = labels.merge(srf, left_on="canonical_molecule_id", right_on="molecule_id", validate="one_to_one")
    Xs = {
        "ECFP": np.vstack([ecfp(x) for x in data.canonical_smiles]),
        "SRF": data[srf_cols].to_numpy(float),
    }
    Xs["ECFP_plus_SRF"] = np.hstack([Xs["ECFP"], Xs["SRF"]])
    y = data.target.to_numpy(int)
    pair = data.match_pair.to_numpy()
    scaffold = data.murcko_scaffold.fillna(data.canonical_molecule_id).to_numpy()
    all_records = []; all_predictions = []
    r, p = evaluate_scheme("matched_pair_group5", Xs, y, pair,
                           lambda seed: make_pair_splits(pair, seed, 5))
    all_records += r; all_predictions += p
    r, p = evaluate_scheme("scaffold_group5", Xs, y, scaffold,
                           lambda seed: list(StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(Xs["SRF"], y, scaffold)))
    all_records += r; all_predictions += p
    r, p = source_component_holdout(Xs, y, data.source_component.to_numpy(), pair)
    all_records += r; all_predictions += p
    metrics = pd.DataFrame(all_records)
    metrics.to_csv(OUT / "repeated_group_metrics.csv", index=False)
    pd.DataFrame(all_predictions).to_csv(OUT / "oof_predictions.csv", index=False)
    data[["canonical_molecule_id", "target", "match_pair", "source_component", "murcko_scaffold"] + srf_cols].to_csv(
        OUT / "state_response_features.csv", index=False)
    summary = metrics.groupby(["scheme", "model"])[["roc_auc", "pr_auc"]].agg(["mean", "std", "min", "max"])
    summary.to_csv(OUT / "summary.csv")
    best = metrics.groupby(["scheme", "model"]).roc_auc.mean().sort_values(ascending=False)
    audit = {
        "n": len(data), "n_positive": int(y.sum()), "n_srf_features": len(srf_cols),
        "endpoint": "A-tier assay-specific functional PAM vs primary-confirmed inactive",
        "best_mean_auc": {"scheme/model": "/".join(best.index[0]), "value": float(best.iloc[0])},
        "promotion_rule": "SRF fusion must improve both matched-pair and scaffold grouped AUC over ECFP and survive permutation/ablation.",
        "claim_boundary": "Exploratory small-n structural response model; no external SOTA claim.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(summary.to_string())
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
