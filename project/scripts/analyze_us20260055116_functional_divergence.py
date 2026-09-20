"""Quantify species/readout divergence and local ordinal cliffs in US20260055116."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from scipy.stats import binomtest, spearmanr


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "pacer_external_us20260055116_v01" / "external_unique_active_moieties.csv"
OUT = PROJECT / "results" / "us20260055116_functional_divergence_v01"
ORDER = {"A": 3, "B": 2, "C": 1, "D": 0}
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def paired_summary(frame: pd.DataFrame, endpoint: str) -> dict:
    part = frame[frame[endpoint].isin(ORDER)].copy()
    human = part.human_m4_perk.map(ORDER).to_numpy(int)
    other = part[endpoint].map(ORDER).to_numpy(int)
    delta = human - other
    non_ties = delta[delta != 0]
    return {
        "n": int(len(part)),
        "spearman": float(spearmanr(human, other).statistic),
        "exact_bin_agreement": float(np.mean(delta == 0)),
        "human_stronger_bins": int(np.sum(delta > 0)),
        "other_stronger_bins": int(np.sum(delta < 0)),
        "ties": int(np.sum(delta == 0)),
        "median_human_minus_other_bin": float(np.median(delta)),
        "two_sided_sign_test_p_non_ties": float(
            binomtest(int(np.sum(non_ties > 0)), len(non_ties), .5).pvalue)
            if len(non_ties) else 1.0,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(SOURCE)
    fps = [FPGEN.GetFingerprint(Chem.MolFromSmiles(x)) for x in frame.canonical_smiles]
    mw = np.asarray([Descriptors.MolWt(Chem.MolFromSmiles(x)) for x in frame.canonical_smiles])
    pairs = []
    for i in range(len(frame)):
        for j in range(i + 1, len(frame)):
            similarity = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            if similarity < .55 or abs(mw[i] - mw[j]) > 90:
                continue
            item = {"compound_i": frame.iloc[i].compound_id,
                    "compound_j": frame.iloc[j].compound_id,
                    "tanimoto": similarity, "abs_MW_delta": abs(mw[i] - mw[j])}
            for endpoint in ["human_m4_perk", "rat_m4_perk", "human_m4_gtpgs"]:
                a, b = frame.iloc[i][endpoint], frame.iloc[j][endpoint]
                item[f"{endpoint}_i"] = a
                item[f"{endpoint}_j"] = b
                item[f"{endpoint}_bin_delta"] = (
                    ORDER[a] - ORDER[b] if a in ORDER and b in ORDER else np.nan)
            pairs.append(item)
    pair_frame = pd.DataFrame(pairs)
    if len(pair_frame):
        pair_frame["max_abs_endpoint_bin_delta"] = pair_frame[
            [c for c in pair_frame if c.endswith("_bin_delta")]].abs().max(axis=1)
        pair_frame = pair_frame.sort_values(
            ["max_abs_endpoint_bin_delta", "tanimoto"], ascending=[False, False])
    pair_frame.to_csv(OUT / "local_functional_cliffs.csv", index=False)

    summary = {
        "human_pERK_vs_rat_pERK": paired_summary(frame, "rat_m4_perk"),
        "human_pERK_vs_human_GTPgammaS": paired_summary(frame, "human_m4_gtpgs"),
        "local_pair_definition": "ECFP4 Tanimoto >=0.55 and absolute MW delta <=90 Da",
        "n_local_pairs": int(len(pair_frame)),
        "n_local_pairs_with_ge2_bin_endpoint_cliff": int(
            (pair_frame.max_abs_endpoint_bin_delta >= 2).sum()) if len(pair_frame) else 0,
        "claim_boundary": (
            "Ordinal bins establish context dependence, not pathway bias, probe dependence, "
            "or molecular mechanism without matched quantitative assays."),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = ["# US20260055116 functional divergence", "",
              "## Paired endpoint summary", "", "```json",
              json.dumps(summary, indent=2), "```", "",
              "## Highest local endpoint cliffs", "",
              (pair_frame.head(15).to_markdown(index=False, floatfmt=".3f")
               if len(pair_frame) else "No local pairs met the frozen definition."), "",
              summary["claim_boundary"]]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if len(pair_frame):
        print(pair_frame.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
