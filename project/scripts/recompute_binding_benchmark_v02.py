#!/usr/bin/env python3
"""Protocol-safe recomputation of binding-line metrics from raw predictions.

Unlike v01 this script recomputes every reported metric with one implementation,
checks row identity before ensembling, and records hashes for external artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score

KEYS = ["target", "label", "canonical_smiles", "murcko_scaffold"]
METRICS = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def metric(y: np.ndarray, score: np.ndarray) -> dict:
    order = np.argsort(-score, kind="mergesort")
    ranked = [[float(score[i]), int(y[i])] for i in order]
    out = {
        "roc_auc": float(roc_auc_score(y, score)),
        "pr_auc": float(average_precision_score(y, score)),
        "bedroc_alpha20": float(CalcBEDROC(ranked, 1, 20.0)),
    }
    for fraction, name in ((0.01, "ef1pct"), (0.05, "ef5pct")):
        n = max(1, int(np.ceil(len(y) * fraction)))
        out[name] = float(y[order[:n]].mean() / y.mean())
    return out


def evaluate(frame: pd.DataFrame, score_col: str) -> dict:
    per = {}
    for target, group in frame.groupby("target", sort=True):
        per[str(target)] = metric(group.label.to_numpy(int), group[score_col].to_numpy(float))
    macro = {name: float(np.mean([row[name] for row in per.values()])) for name in METRICS}
    return {"macro": macro, "per_target": per, "n_rows": int(len(frame)), "n_targets": len(per)}


def read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    missing = set(KEYS) - set(frame)
    if missing:
        raise ValueError(f"{path}: missing identity columns {sorted(missing)}")
    return frame


def ensemble(paths: list[Path], score_col: str) -> pd.DataFrame:
    frames = [read(path) for path in paths]
    base = frames[0][KEYS].reset_index(drop=True)
    scores = []
    for path, frame in zip(paths, frames):
        identity = frame[KEYS].reset_index(drop=True)
        if not identity.equals(base):
            raise ValueError(f"row identity/order mismatch: {path}")
        scores.append(frame[score_col].to_numpy(float))
    base = base.copy()
    base["ensemble_score"] = np.mean(scores, axis=0)
    return base


def target_bootstrap_delta(candidate: dict, baseline: dict, draws: int, seed: int) -> dict:
    targets = sorted(set(candidate["per_target"]) & set(baseline["per_target"]))
    rng = np.random.default_rng(seed)
    output = {}
    for name in METRICS:
        delta = np.asarray([candidate["per_target"][t][name] - baseline["per_target"][t][name] for t in targets])
        boot = delta[rng.integers(0, len(delta), size=(draws, len(delta)))].mean(axis=1)
        output[name] = {
            "delta": float(delta.mean()),
            "target_bootstrap_95ci": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
        }
    return output


def locate(relative: str, roots: list[Path]) -> Path:
    for root in roots:
        path = root / relative
        if path.is_file():
            return path
    raise FileNotFoundError(f"artifact absent under all roots: {relative}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument("--artifact-root", type=Path, help="Optional checkout containing git-ignored raw predictions")
    p.add_argument("--output", type=Path)
    p.add_argument("--bootstrap", type=int, default=20000)
    p.add_argument("--seed", type=int, default=20261001)
    args = p.parse_args()
    roots = [args.repo_root]
    if args.artifact_root and args.artifact_root.resolve() != args.repo_root.resolve():
        roots.append(args.artifact_root)
    base = "project/results/drugclip_science2026"
    methods, provenance = {}, {}

    def register(protocol: str, method: str, path: Path, score: str) -> None:
        frame = read(path)
        methods[f"{protocol}/{method}"] = evaluate(frame, score)
        provenance[f"{protocol}/{method}"] = {"files": [str(path)], "sha256": [sha256(path)], "score_column": score}

    # 13-target fixed protocol.
    seed13 = [locate(f"{base}/migration_13target_v01/loso_preserve0.2_seed{s}.predictions.csv", roots)
              for s in (20260925, 20260926, 20260927)]
    register("13T", "official_drugclip", seed13[0], "official")
    for s, path in zip((20260925, 20260926, 20260927), seed13):
        register("13T", f"ep80_seed{s}", path, "tuned")
    ep80 = ensemble(seed13, "tuned")
    methods["13T/ep80_score_ensemble"] = evaluate(ep80, "ensemble_score")
    provenance["13T/ep80_score_ensemble"] = {"files": [str(x) for x in seed13], "sha256": [sha256(x) for x in seed13]}
    ep40p = locate(f"{base}/migration_13target_v01/summary_preserve0_3seed.predictions.csv", roots)
    register("13T", "ep40_score_ensemble", ep40p, "tuned")
    ecfp13 = locate(f"{base}/migration_13target_v01/ecfp_loso.predictions.csv", roots)
    register("13T", "ecfp4_logistic", ecfp13, "ecfp4_logistic")

    fam13 = [locate(f"{base}/family_aug_v01/famaug_ep80_seed{s}.predictions.csv", roots)
             for s in (20260925, 20260926, 20260927)]
    for s, path in zip((20260925, 20260926, 20260927), fam13):
        register("13T", f"famaug_seed{s}", path, "tuned")
    fam_ens = ensemble(fam13, "tuned")
    methods["13T/famaug_score_ensemble"] = evaluate(fam_ens, "ensemble_score")
    provenance["13T/famaug_score_ensemble"] = {"files": [str(x) for x in fam13], "sha256": [sha256(x) for x in fam13]}

    # 20-target fixed protocol; never compare absolute values as if it were 13T.
    seed20 = [locate(f"{base}/extended20_finetune_v01/loso20_ep80_seed{s}.predictions.csv", roots)
              for s in (20260925, 20260926, 20260927)]
    register("20T", "official_drugclip", seed20[0], "official")
    for s, path in zip((20260925, 20260926, 20260927), seed20):
        register("20T", f"ep80_seed{s}", path, "tuned")
    ep20 = ensemble(seed20, "tuned")
    methods["20T/ep80_score_ensemble"] = evaluate(ep20, "ensemble_score")
    provenance["20T/ep80_score_ensemble"] = {"files": [str(x) for x in seed20], "sha256": [sha256(x) for x in seed20]}
    ecfp20 = locate(f"{base}/extended20_finetune_v01/ecfp20_loso.predictions.csv", roots)
    register("20T", "ecfp4_logistic", ecfp20, "ecfp4_logistic")
    fam20 = [locate(f"{base}/family_aug_v01/famaug20_ep80_seed{s}.predictions.csv", roots)
             for s in (20260925, 20260926, 20260927)]
    for s, path in zip((20260925, 20260926, 20260927), fam20):
        register("20T", f"famaug_seed{s}", path, "tuned")
    fam20e = ensemble(fam20, "tuned")
    methods["20T/famaug_score_ensemble"] = evaluate(fam20e, "ensemble_score")
    provenance["20T/famaug_score_ensemble"] = {"files": [str(x) for x in fam20], "sha256": [sha256(x) for x in fam20]}

    comparisons = {
        "13T/famaug_vs_official": target_bootstrap_delta(methods["13T/famaug_score_ensemble"], methods["13T/official_drugclip"], args.bootstrap, args.seed),
        "13T/famaug_vs_ecfp": target_bootstrap_delta(methods["13T/famaug_score_ensemble"], methods["13T/ecfp4_logistic"], args.bootstrap, args.seed + 1),
        "13T/famaug_vs_ep80": target_bootstrap_delta(methods["13T/famaug_score_ensemble"], methods["13T/ep80_score_ensemble"], args.bootstrap, args.seed + 2),
        "20T/ep80_vs_official": target_bootstrap_delta(methods["20T/ep80_score_ensemble"], methods["20T/official_drugclip"], args.bootstrap, args.seed + 3),
        "20T/famaug_vs_ep80": target_bootstrap_delta(methods["20T/famaug_score_ensemble"], methods["20T/ep80_score_ensemble"], args.bootstrap, args.seed + 4),
    }
    result = {
        "version": "v02", "metric_implementation": "sklearn ROC/AP + RDKit CalcBEDROC(alpha=20) + ceil EF",
        "protocol_rule": "13T and 20T are separate benchmarks; no cross-protocol winner selection",
        "methods": methods, "paired_target_bootstrap": comparisons, "provenance": provenance,
        "supersedes": "v01 ensemble BEDROC/EF values computed by ensemble_eval_13target.py",
    }
    output = args.output or args.repo_root / "project/results/project_wide_integration_benchmark_v02/binding_recomputed_v02.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    rows = []
    for key, node in methods.items():
        protocol, method = key.split("/", 1)
        rows.append({"protocol": protocol, "method": method, **node["macro"]})
    pd.DataFrame(rows).to_csv(output.with_suffix(".csv"), index=False)
    print(json.dumps({"output": str(output), "methods": len(methods),
                      "famaug13": methods["13T/famaug_score_ensemble"]["macro"],
                      "ep80_20": methods["20T/ep80_score_ensemble"]["macro"]}, indent=2))


if __name__ == "__main__":
    main()
