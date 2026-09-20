# -*- coding: utf-8 -*-
"""Build a descriptor-matched experimental PAM/inactive validation set."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.optimize import linear_sum_assignment


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "pam_vs_inactive.csv"
OUT = P / "results" / "m4_gamd_ensemble" / "validation_set.csv"
RAW = P.parent / "CHRM4_PAM_Modeling_Handoff_v1.0" / "data" / "modeling_pam_vs_inactive.csv"


def desc(smiles):
    m = Chem.MolFromSmiles(smiles)
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m), Descriptors.TPSA(m), Descriptors.HeavyAtomCount(m)]


def build(d, out_path, scope):
    pos = d[d.target == 1].copy().reset_index(drop=True)
    neg = d[d.target == 0].copy().reset_index(drop=True)
    xp = np.asarray([desc(x) for x in pos.canonical_smiles], float)
    xn = np.asarray([desc(x) for x in neg.canonical_smiles], float)
    scale = np.std(np.vstack([xp, xn]), axis=0); scale[scale < 1e-8] = 1
    cost = np.linalg.norm((xn[:, None, :] - xp[None, :, :]) / scale, axis=2)
    ni, pi = linear_sum_assignment(cost)
    matched_pos = pos.iloc[pi].copy(); matched_neg = neg.iloc[ni].copy()
    matched_pos["match_pair"] = np.arange(len(pi)); matched_neg["match_pair"] = np.arange(len(ni))
    out = pd.concat([matched_pos, matched_neg], ignore_index=True)
    out = out[["canonical_molecule_id", "canonical_smiles", "target", "match_pair", "source_component", "murcko_scaffold"]]
    out.to_csv(out_path, index=False)
    audit = {
        "n_PAM": int((out.target == 1).sum()), "n_experimental_inactive": int((out.target == 0).sum()),
        "matching": "one-to-one Hungarian matching on standardized MW, logP, TPSA, heavy-atom count",
        "mean_standardized_pair_distance": float(cost[ni, pi].mean()),
        "scope": scope,
        "purpose": "experimental negative validation; no artificial decoys",
    }
    out_path.with_name(out_path.stem + "_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


def main():
    d = pd.read_csv(DATA)
    build(d, OUT, "all strict labels, including database-text negatives")
    raw = pd.read_csv(RAW, low_memory=False)
    primary_ids = set(raw.loc[
        raw.negative_evidence_tier.eq("A_primary_confirmed_inactive"), "canonical_molecule_id"
    ].dropna())
    primary = d[(d.target == 1) | d.canonical_molecule_id.isin(primary_ids)].copy()
    build(primary, OUT.with_name("validation_set_primary.csv"),
          "A-tier primary-source-confirmed inactive only")


if __name__ == "__main__":
    main()
