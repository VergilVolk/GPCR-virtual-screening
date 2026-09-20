# -*- coding: utf-8 -*-
"""Analyze M4R GaMD ensemble docking without claiming functional PAM efficacy."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

P = Path(__file__).resolve().parents[1]
ROOT = P / "results" / "m4_gamd_ensemble"
PMF_FILE = P / "data" / "miao2026_m4r" / "M4R_ensemble_PMF.xvg"


def load_pmf(path: Path) -> dict[int, float]:
    values = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith(("#", "@")):
            continue
        x, y, *_ = line.split()
        cluster = int(float(x))
        if cluster >= 0:
            values[cluster] = float(y)
    return values


def ef_at_fraction(y: np.ndarray, score: np.ndarray, fraction: float) -> dict:
    n = len(y); k = max(1, int(np.ceil(n * fraction)))
    order = np.argsort(-score, kind="stable")
    hits = int(y[order[:k]].sum())
    base = float(y.mean())
    return {"fraction": fraction, "n_selected": k, "hits": hits,
            "ef": float((hits / k) / base) if base else float("nan")}


def metrics(y: np.ndarray, lower_is_better: np.ndarray) -> dict:
    rank_score = -np.asarray(lower_is_better, float)
    return {
        "roc_auc": float(roc_auc_score(y, rank_score)),
        "pr_auc": float(average_precision_score(y, rank_score)),
        "ef_10pct": ef_at_fraction(y, rank_score, 0.10),
        "ef_05pct": ef_at_fraction(y, rank_score, 0.05),
    }


def stratified_paired_bootstrap(df: pd.DataFrame, a: str, b: str,
                                n_boot: int = 4000, seed: int = 20260829) -> dict:
    """Paired bootstrap CI for AUC(b)-AUC(a), preserving class balance."""
    rng = np.random.default_rng(seed)
    pos = np.flatnonzero(df.target.to_numpy() == 1)
    neg = np.flatnonzero(df.target.to_numpy() == 0)
    diffs = []
    for _ in range(n_boot):
        idx = np.r_[rng.choice(pos, len(pos), replace=True),
                    rng.choice(neg, len(neg), replace=True)]
        y = df.target.to_numpy()[idx]
        diffs.append(roc_auc_score(y, -df[b].to_numpy()[idx]) -
                     roc_auc_score(y, -df[a].to_numpy()[idx]))
    lo, hi = np.quantile(diffs, [0.025, 0.975])
    return {"delta_auc": float(np.mean(diffs)), "ci95": [float(lo), float(hi)],
            "p_improvement": float(np.mean(np.asarray(diffs) > 0)), "n_boot": n_boot}


def aggregate(frame: pd.DataFrame, expected_clusters: int = 10) -> pd.DataFrame:
    pmf = load_pmf(PMF_FILE)
    frame = frame.copy()
    frame["pmf"] = frame.cluster.map(pmf)
    frame["BE"] = frame.vina_affinity + frame.pmf
    rows = []
    for mid, part in frame.groupby("molecule_id"):
        part = part.sort_values("cluster")
        best = part.loc[part.BE.idxmin()]
        rows.append({
            "molecule_id": mid,
            "n_clusters": int(part.cluster.nunique()),
            "cluster00_vina": float(part.loc[part.cluster.eq(0), "vina_affinity"].iloc[0])
                if part.cluster.eq(0).any() else np.nan,
            "BE_min": float(part.BE.min()),
            "BE_avg": float(part.BE.mean()),
            "best_cluster": int(best.cluster),
            "vina_mean": float(part.vina_affinity.mean()),
            "vina_sd": float(part.vina_affinity.std(ddof=0)),
            "vina_range": float(part.vina_affinity.max() - part.vina_affinity.min()),
            "pocket_coverage_mean": float(part.native_pocket_coverage.mean()),
            "pocket_coverage_min": float(part.native_pocket_coverage.min()),
            "complete_ensemble": bool(part.cluster.nunique() == expected_clusters),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=4000)
    args = ap.parse_args()
    source = ROOT / "validation_primary_docking_protocol_v127_ph7_rank1" / "framewise_scores.csv"
    frame = pd.read_csv(source)
    frame = frame[frame.error.isna() | frame.error.eq("")].dropna(subset=["vina_affinity"])
    agg = aggregate(frame)
    val = pd.read_csv(ROOT / "validation_set_primary.csv")
    out = agg.merge(val, left_on="molecule_id", right_on="canonical_molecule_id", validate="one_to_one")
    out = out[out.complete_ensemble].copy()
    out.to_csv(ROOT / "validation_ensemble_scores.csv", index=False)

    summary = {
        "n_complete": int(len(out)),
        "n_positive": int(out.target.sum()),
        "n_negative": int((1 - out.target).sum()),
        "score_definition": "BE(cluster) = Vina affinity + GaMD PMF; BE_min=min; BE_avg=arithmetic mean",
        "claim_boundary": "Measures binding/enrichment evidence, not PAM efficacy or cooperativity.",
        "metrics": {}, "paired_bootstrap_vs_cluster00": {},
    }
    y = out.target.to_numpy(int)
    for col in ["cluster00_vina", "BE_min", "BE_avg"]:
        summary["metrics"][col] = metrics(y, out[col].to_numpy(float))
    for col in ["BE_min", "BE_avg"]:
        summary["paired_bootstrap_vs_cluster00"][col] = stratified_paired_bootstrap(
            out, "cluster00_vina", col, args.n_boot)
    (ROOT / "validation_ensemble_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    report = [
        "# M4R GaMD ensemble docking validation", "",
        f"- Complete molecules: {len(out)} ({int(out.target.sum())} PAM / {int((1-out.target).sum())} inactive)",
        "- Reweighting: `BE_i = Vina_i + PMF_i`; `BEmin = min(BE_i)`; `BEavg = mean(BE_i)`.",
        "- Boundary: this validates structural enrichment only; it does not prove PAM efficacy, cooperativity, or probe dependence.", "",
        "| Method | ROC-AUC | PR-AUC | EF10% |", "|---|---:|---:|---:|",
    ]
    for col, label in [("cluster00_vina", "Single GaMD cluster 00"), ("BE_min", "GaMD BEmin"), ("BE_avg", "GaMD BEavg")]:
        m = summary["metrics"][col]
        report.append(f"| {label} | {m['roc_auc']:.3f} | {m['pr_auc']:.3f} | {m['ef_10pct']['ef']:.2f} |")
    report += ["", "## Paired bootstrap vs cluster 00", ""]
    for col in ["BE_min", "BE_avg"]:
        b = summary["paired_bootstrap_vs_cluster00"][col]
        report.append(f"- {col}: ΔAUC={b['delta_auc']:+.3f}, 95% CI [{b['ci95'][0]:+.3f}, {b['ci95'][1]:+.3f}], P(Δ>0)={b['p_improvement']:.3f}.")
    (ROOT / "VALIDATION_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
