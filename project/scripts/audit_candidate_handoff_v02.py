#!/usr/bin/env python3
"""Audit the actual PACER-200 handoff between old and family-augmented scores."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


def sha256(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return h


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--old-ranking", type=Path, required=True)
    p.add_argument("--famaug-scores", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    old, new = pd.read_csv(args.old_ranking), pd.read_csv(args.famaug_scores)
    required_old = {"candidate_id", "canonical_smiles", "drugclip_rank", "fusion_rank"}
    required_new = {"candidate_id", "candidate_smiles", "new_rank", "famaug_m4_score"}
    if required_old - set(old) or required_new - set(new):
        raise ValueError("candidate handoff columns are incomplete")
    merged = old.merge(new, on="candidate_id", how="outer", indicator=True, suffixes=("_old", "_new"))
    smiles_match = merged.eval("canonical_smiles == candidate_smiles").fillna(False)
    matched = merged.loc[merged._merge.eq("both") & smiles_match].copy()
    rho = float(spearmanr(matched.drugclip_rank, matched.new_rank).statistic)
    overlaps = {}
    for k in (20, 50, 100):
        a = set(matched.nsmallest(k, "drugclip_rank").candidate_id)
        b = set(matched.nsmallest(k, "new_rank").candidate_id)
        overlaps[str(k)] = {"intersection": len(a & b), "jaccard": len(a & b) / len(a | b)}
    dual = matched[(matched.drugclip_rank <= 50) & (matched.new_rank <= 25)].copy()
    dual = dual.sort_values(["new_rank", "drugclip_rank"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matched.sort_values("new_rank").to_csv(args.output_dir / "candidate_handoff_all200_v02.csv", index=False)
    dual.to_csv(args.output_dir / "candidate_dual_model_shortlist_v02.csv", index=False)
    report = {
        "status": "PASS" if len(matched) == 200 and smiles_match.sum() == 200 else "FAIL",
        "old_rows": int(len(old)), "new_rows": int(len(new)), "matched_ids_and_smiles": int(len(matched)),
        "spearman_old_drugclip_vs_famaug_rank": rho, "topk_overlap": overlaps,
        "dual_rule": "old DrugCLIP rank <= 50 AND family-augmented rank <= 25",
        "dual_count": int(len(dual)), "dual_candidate_ids": dual.candidate_id.tolist(),
        "provenance": {"old": str(args.old_ranking), "old_sha256": sha256(args.old_ranking),
                       "new": str(args.famaug_scores), "new_sha256": sha256(args.famaug_scores)},
        "correction": "v01 pipeline_continuity_check compared PACER-200 to an unrelated 28,519-molecule library and is superseded.",
    }
    (args.output_dir / "candidate_handoff_audit_v02.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
