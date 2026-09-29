#!/usr/bin/env python3
"""Frozen, label-free fusion of DrugCLIP and ECFP-to-pocket transfer scores.

Scale factors are estimated on the non-M4 strict panel.  The M4 labels are
used once, for final evaluation only; no M4 threshold or fusion weight is fit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


KEYS = ["canonical_smiles", "positive_subtype", "negative_subtype"]


def collapse_drugclip(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    return frame.groupby(KEYS, as_index=False).agg(
        drugclip_delta=("targeted_delta", "mean")
    )


def accuracy(x: np.ndarray) -> float:
    return float(np.mean(x > 0))


def clustered_delta_ci(frame: pd.DataFrame, left: str, right: str,
                       seed: int = 20260925, draws: int = 10000) -> list[float]:
    rng = np.random.default_rng(seed)
    molecules = frame.canonical_smiles.unique()
    grouped = {m: frame.index[frame.canonical_smiles == m].to_numpy() for m in molecules}
    values = np.empty(draws, dtype=float)
    for i in range(draws):
        sampled = rng.choice(molecules, size=len(molecules), replace=True)
        idx = np.concatenate([grouped[m] for m in sampled])
        values[i] = accuracy(frame.loc[idx, left].to_numpy()) - accuracy(frame.loc[idx, right].to_numpy())
    return [float(v) for v in np.quantile(values, [0.025, 0.5, 0.975])]


def paired_sign_p(frame: pd.DataFrame, left: str, right: str,
                  seed: int = 20260925, draws: int = 10000) -> float:
    """Molecule-level paired sign-flip permutation test."""
    rng = np.random.default_rng(seed)
    per_mol = frame.groupby("canonical_smiles").apply(
        lambda g: accuracy(g[left].to_numpy()) - accuracy(g[right].to_numpy()),
        include_groups=False,
    ).to_numpy()
    observed = float(per_mol.mean())
    null = np.asarray([(per_mol * rng.choice([-1.0, 1.0], len(per_mol))).mean()
                       for _ in range(draws)])
    return float((1 + np.sum(null >= observed)) / (draws + 1))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--drugclip-strict", type=Path, required=True)
    p.add_argument("--drugclip-m4", type=Path, required=True)
    p.add_argument("--ecfp-strict", type=Path, required=True)
    p.add_argument("--ecfp-m4", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    ds = collapse_drugclip(args.drugclip_strict)
    dm = pd.read_csv(args.drugclip_m4)[KEYS + ["gpcr_triplet_delta"]].rename(
        columns={"gpcr_triplet_delta": "drugclip_delta"})
    es = pd.read_csv(args.ecfp_strict)[KEYS + ["ecfp_pocket_delta"]]
    em = pd.read_csv(args.ecfp_m4)
    strict = ds.merge(es, on=KEYS, validate="one_to_one")
    m4 = dm.merge(em, on=KEYS, validate="one_to_one")

    # Positive homogeneous scaling preserves each model's zero decision boundary.
    scales = {
        "drugclip": float(np.sqrt(np.mean(strict.drugclip_delta.to_numpy() ** 2))),
        "ecfp_pocket": float(np.sqrt(np.mean(strict.ecfp_pocket_delta.to_numpy() ** 2))),
    }
    m4["drugclip_scaled"] = m4.drugclip_delta / scales["drugclip"]
    m4["ecfp_scaled"] = m4.ecfp_pocket_delta / scales["ecfp_pocket"]
    m4["fusion_delta"] = 0.5 * (m4.drugclip_scaled + m4.ecfp_scaled)

    metrics = {}
    for name, col in [("drugclip_triplet", "drugclip_delta"),
                      ("ecfp_to_pocket", "ecfp_pocket_delta"),
                      ("fixed_equal_fusion", "fusion_delta")]:
        metrics[name] = {
            "pair_accuracy": accuracy(m4[col].to_numpy()),
            "balanced_direction_accuracy": float(m4.assign(ok=m4[col] > 0).groupby("m4_direction").ok.mean().mean()),
            "delta_vs_pchembl_gap_spearman": float(spearmanr(m4[col], m4.delta_pchembl).statistic),
        }

    report = {
        "protocol": "M4 labels untouched until final evaluation; fixed 0.5/0.5 fusion; scales frozen on non-M4 strict panel",
        "n_pairs": int(len(m4)),
        "n_molecules": int(m4.canonical_smiles.nunique()),
        "strict_panel_rms_scales": scales,
        "metrics": metrics,
        "fusion_minus_drugclip_molecule_cluster_95ci": clustered_delta_ci(m4, "fusion_delta", "drugclip_delta"),
        "fusion_minus_ecfp_molecule_cluster_95ci": clustered_delta_ci(m4, "fusion_delta", "ecfp_pocket_delta"),
        "fusion_vs_drugclip_one_sided_molecule_signflip_p": paired_sign_p(m4, "fusion_delta", "drugclip_delta"),
        "fusion_vs_ecfp_one_sided_molecule_signflip_p": paired_sign_p(m4, "fusion_delta", "ecfp_pocket_delta"),
        "drugclip_ecfp_spearman": float(spearmanr(m4.drugclip_delta, m4.ecfp_pocket_delta).statistic),
        "claim_boundary": "Unseen-M4 subtype activity ranking; not binding affinity, allostery, PAM identity, or efficacy.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    m4.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
