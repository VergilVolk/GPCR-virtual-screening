"""Matched-study specificity test of the frozen M4 PAM structural signature."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import test_consensus_pam_signature as sigmod
from validate_external_static_pam_signature import ca_coords, frozen_score


PROJECT = Path(__file__).resolve().parents[1]
EXT = PROJECT / "data" / "pdb" / "m4_external_structures"
OUT = PROJECT / "results" / "pacer_matched_static_specificity_v02"
STRUCTURES = {
    "7V68_LY2119620_PAM": ("7V68.pdb", "PAM"),
    "7V69_iperoxo_only": ("7V69.pdb", "no_PAM"),
    "7V6A_compound110": ("7V6A.pdb", "allosteric_agonist"),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    q, p, s, _ = sigmod.load_static()
    feats = {name: sigmod.pair_features(ca_coords(EXT / file, "R"))
             for name, (file, _) in STRUCTURES.items()}
    rows = []
    for threshold in (.05, .10, .20):
        for name, x in feats.items():
            score, n = frozen_score(q, p, s, x, threshold)
            rows.append({"threshold_A": threshold, "structure": name,
                         "role": STRUCTURES[name][1], "signature_score": score,
                         "features": n})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "matched_structure_scores.csv", index=False)
    checks = []
    for threshold in (.05, .10, .20):
        v = table[table.threshold_A == threshold].set_index("structure").signature_score
        pam, no, agonist = (float(v[k]) for k in STRUCTURES)
        checks.append({"threshold_A": threshold, "PAM_above_matched_no_PAM": pam > no,
                       "PAM_above_allosteric_agonist": pam > agonist,
                       "PAM_specificity_supported": pam > no and pam > agonist})
    checks = pd.DataFrame(checks)
    checks.to_csv(OUT / "specificity_checks.csv", index=False)
    audit = {"signature_refit": False, "same_study_same_construct": True,
             "probe_matched_PAM_vs_no_PAM": True,
             "allosteric_agonist_specificity_control": True,
             "matched_probe_support": bool(checks.PAM_above_matched_no_PAM.all()),
             "PAM_specificity_support": bool(checks.PAM_specificity_supported.all()),
             "claim_boundary": "Independent matched static specificity; not efficacy prediction."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = "# PACER Matched Static PAM-Specificity v0.2\n\n" + table.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n" + checks.to_markdown(index=False) + "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n"
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(table.to_string(index=False)); print(checks.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
