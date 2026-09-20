"""Blind external validation of the frozen cross-PAM structural signature."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import test_consensus_pam_signature as sigmod


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "pacer_external_static_pam_signature_v01"
EXT = PROJECT / "data" / "pdb" / "m4_external_structures"
STRUCTURES = {
    "7V68_LY2119620_PAM": (EXT / "7V68.pdb", "R", "PAM", "iperoxo"),
    "7TRK_no_PAM": (EXT / "7TRK.pdb", "R", "no_PAM", "iperoxo"),
    "5DSG_inactive": (EXT / "5DSG.pdb", "A", "no_PAM_inactive", "tiotropium"),
}
THRESHOLDS = [0.05, 0.10, 0.20]


def ca_coords(path: Path, chain: str) -> np.ndarray:
    atoms = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM  ") or line[21:22] != chain or line[12:16].strip() != "CA":
            continue
        atoms[int(line[22:26])] = (line[17:20].strip(), np.array([
            float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    missing = [r for r in sigmod.STATIC_RESIDS if r not in atoms]
    if missing:
        raise RuntimeError(f"{path.stem}/{chain}: missing {missing}")
    return np.stack([atoms[r][1] for r in sigmod.STATIC_RESIDS])


def frozen_score(q, p, s, x, threshold):
    dq, dp = q - s, p - s
    mask = (dq * dp > 0) & (np.abs(dq) >= threshold) & (np.abs(dp) >= threshold)
    d = (dq + dp) / 2
    agreement = np.minimum(np.abs(dq), np.abs(dp)) / np.maximum(np.abs(dq), np.abs(dp))
    v = d[mask] * np.sqrt(agreement[mask])
    return float(np.dot((x[mask] - s[mask]) * np.sqrt(agreement[mask]), v) / np.dot(v, v)), int(mask.sum())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    q, p, s, _ = sigmod.load_static()
    external = {name: sigmod.pair_features(ca_coords(path, chain))
                for name, (path, chain, _, _) in STRUCTURES.items()}
    rows = []
    for threshold in THRESHOLDS:
        for name, values in external.items():
            score, nfeat = frozen_score(q, p, s, values, threshold)
            _, _, label, probe = STRUCTURES[name]
            rows.append({"threshold_A": threshold, "structure": name, "label": label,
                         "orthosteric_probe": probe, "signature_score": score,
                         "signature_features": nfeat})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "external_structure_scores.csv", index=False)
    conditions = []
    for threshold in THRESHOLDS:
        sub = table[table.threshold_A == threshold].set_index("structure").signature_score
        pam = float(sub["7V68_LY2119620_PAM"])
        ipx = float(sub["7TRK_no_PAM"])
        inactive = float(sub["5DSG_inactive"])
        conditions.append({"threshold_A": threshold, "pam_positive": pam > 0,
                           "pam_above_probe_matched_no_pam": pam > ipx,
                           "pam_above_inactive": pam > inactive,
                           "all_conditions": pam > 0 and pam > ipx and pam > inactive})
    cond = pd.DataFrame(conditions)
    cond.to_csv(OUT / "validation_conditions.csv", index=False)
    audit = {
        "external_structures": 3,
        "signature_refit": False,
        "primary_threshold_A": .10,
        "primary_supported": bool(cond.set_index("threshold_A").loc[.10, "all_conditions"]),
        "sensitivity_supported": bool(cond.all_conditions.all()),
        "external_structural_direction_supported": bool(cond.all_conditions.all()),
        "claim_boundary": "Independent static structural direction; not PAM efficacy or dynamic-state validation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = "# PACER External Static PAM-Signature Validation v0.1\n\n"
    report += table.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n## Frozen conditions\n\n" + cond.to_markdown(index=False)
    report += "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n"
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(table.to_string(index=False)); print(cond.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
