"""Nested development test of the frozen adaptive anchor-budget policy."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import adapted, base_model, choose_diverse, features, safe_rho, train_delta


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = PROJECT / "results" / "pacer_adaptive_anchors_v01"
SEED = 20260830
MIN_GROUP = 12
SPAN_THRESHOLD = .75
MAX_ANCHORS = 7


def stable_value(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:16], 16)


def pool_query_split(indices: np.ndarray, ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ordered = sorted(indices, key=lambda i: stable_value(str(ids[i])))
    n_pool = max(7, int(np.floor(.60 * len(ordered))))
    n_pool = min(n_pool, len(ordered) - 4)
    return np.asarray(ordered[:n_pool], int), np.asarray(ordered[n_pool:], int)


def expand_extremes(support, pool, base, y, X, D, bv, delta_pack):
    remaining = np.asarray([i for i in pool if i not in set(support)], int)
    if len(remaining) <= 2:
        return np.r_[support, remaining]
    pred = adapted("DeltaSAR", base, support, remaining, y, X, D, bv, delta_pack)
    order = np.argsort(pred)
    chosen = np.asarray([remaining[order[0]], remaining[order[-1]]], int)
    return np.r_[support, chosen]


def evaluate(name, group, support, query, base, y, X, D, bv, delta_pack, refused=False):
    pred = adapted("DeltaSAR", base, support, query, y, X, D, bv, delta_pack)
    return {
        "group": group, "method": name, "n_group": len(support) + len(query),
        "n_anchor": len(support), "n_query": len(query),
        "anchor_span": float(np.ptp(y[support])), "refused": bool(refused),
        "Spearman": safe_rho(y[query], pred),
        "MAE": float(mean_absolute_error(y[query], pred)),
        "support_ids": "|".join(map(str, support)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles.tolist())
    y = d.pEC50.to_numpy(float)
    groups = d.source_component.astype(str).to_numpy()
    ids = d.canonical_molecule_id.astype(str).to_numpy()
    eligible = [g for g, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []

    for fold, group in enumerate(eligible):
        target = np.where(groups == group)[0]
        train = np.where(groups != group)[0]
        pool, query = pool_query_split(target, ids)
        base = np.full(len(d), np.nan)
        base[target] = base_model(X, D, y, train, target)
        delta_pack = train_delta(train, X, D, y, groups, bv, seed=SEED + fold, cliff_weight=0)

        fixed_supports = {}
        for shot in (3, 5, 7):
            fixed_supports[shot] = choose_diverse(pool, min(shot, len(pool)), bv)
            rows.append(evaluate(f"FixedDiverse-{shot}", group, fixed_supports[shot],
                                 query, base, y, X, D, bv, delta_pack))

        support = fixed_supports[3].copy()
        while np.ptp(y[support]) < SPAN_THRESHOLD and len(support) < min(MAX_ANCHORS, len(pool)):
            support = expand_extremes(support, pool, base, y, X, D, bv, delta_pack)
        refused = bool(np.ptp(y[support]) < SPAN_THRESHOLD)
        rows.append(evaluate("PACER-AdaptiveBracket", group, support, query,
                             base, y, X, D, bv, delta_pack, refused=refused))

    episodes = pd.DataFrame(rows)
    episodes.to_csv(OUT / "episodes.csv", index=False)
    summary = episodes.groupby("method").agg(
        n_series=("group", "size"), macro_Spearman=("Spearman", "mean"),
        median_Spearman=("Spearman", "median"), worst_Spearman=("Spearman", "min"),
        macro_MAE=("MAE", "mean"), mean_anchor_budget=("n_anchor", "mean"),
        refusal_rate=("refused", "mean"), mean_anchor_span=("anchor_span", "mean"),
    ).reset_index().sort_values("macro_Spearman", ascending=False)
    summary.to_csv(OUT / "summary.csv", index=False)
    adaptive = summary[summary.method == "PACER-AdaptiveBracket"].iloc[0]
    fixed3 = summary[summary.method == "FixedDiverse-3"].iloc[0]
    fixed7 = summary[summary.method == "FixedDiverse-7"].iloc[0]
    audit = {
        "status": "historical_nested_development",
        "query_labels_hidden_from_anchor_selection_and_stopping": True,
        "span_threshold_log10": SPAN_THRESHOLD, "maximum_anchors": MAX_ANCHORS,
        "adaptive_macro_delta_vs_fixed3": float(adaptive.macro_Spearman - fixed3.macro_Spearman),
        "adaptive_worst_delta_vs_fixed3": float(adaptive.worst_Spearman - fixed3.worst_Spearman),
        "adaptive_mean_budget": float(adaptive.mean_anchor_budget),
        "fixed7_macro_spearman": float(fixed7.macro_Spearman),
        "preregistered_development_pass": bool(
            adaptive.macro_Spearman > fixed3.macro_Spearman and
            adaptive.worst_Spearman >= fixed3.worst_Spearman and
            adaptive.mean_anchor_budget < 7),
        "claim_boundary": "Historical nested development only; requires untouched campaign validation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER adaptive anchor-budget development", "",
              "Every series used a fixed, never-measured query partition. Only anchor-pool labels could drive expansion.", "",
              summary.to_markdown(index=False, floatfmt=".3f"), "",
              f"Preregistered development pass: **{audit['preregistered_development_pass']}**.", "",
              "This is historical nested development, not external SOTA evidence."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
