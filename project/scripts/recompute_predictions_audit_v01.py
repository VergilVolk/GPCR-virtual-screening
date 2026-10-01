#!/usr/bin/env python
"""Independent metric recomputation from raw prediction CSVs.

The strongest integration check: do NOT trust committed metric JSONs; recompute
per-target ROC-AUC / PR-AUC from the 37k-row prediction files and compare.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from sklearn.metrics import roc_auc_score, average_precision_score

BIND = Path(r"D:\CLC\project\results\drugclip_science2026")
OUT = Path(r"D:\CLC_geom2vec_pilot\project\results\project_wide_integration_benchmark_v01")

results = []


def recompute(csv_path: Path, score_col: str):
    by_target = defaultdict(lambda: ([], []))
    with csv_path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            by_target[row["target"]][0].append(float(row[score_col]))
            by_target[row["target"]][1].append(int(float(row["label"])))
    per_target = {}
    for t, (scores, labels) in sorted(by_target.items()):
        per_target[t] = {
            "roc_auc": roc_auc_score(labels, scores),
            "pr_auc": average_precision_score(labels, scores),
            "n": len(labels),
        }
    macro_roc = sum(v["roc_auc"] for v in per_target.values()) / len(per_target)
    macro_pr = sum(v["pr_auc"] for v in per_target.values()) / len(per_target)
    return {"per_target": per_target, "macro_roc": macro_roc, "macro_pr": macro_pr}


def compare(name: str, recomputed: dict, committed_node: dict):
    committed = dict(committed_node.get("macro", committed_node))
    committed["per_target"] = committed_node.get("per_target", committed.get("per_target", {}))
    ok_roc = abs(recomputed["macro_roc"] - committed["roc_auc"]) < 2e-3
    ok_pr = abs(recomputed["macro_pr"] - committed["pr_auc"]) < 2e-3
    bad_targets = [
        t for t, v in recomputed["per_target"].items()
        if abs(v["roc_auc"] - committed["per_target"][t]["roc_auc"]) > 2e-3
    ] if committed.get("per_target") else []
    results.append({
        "check": name,
        "recomputed_macro_roc": round(recomputed["macro_roc"], 4),
        "recomputed_macro_pr": round(recomputed["macro_pr"], 4),
        "committed_macro_roc": round(committed["roc_auc"], 4),
        "committed_macro_pr": round(committed["pr_auc"], 4),
        "macro_roc_match": ok_roc, "macro_pr_match": ok_pr,
        "per_target_mismatches": bad_targets,
    })
    print(("PASS " if ok_roc and ok_pr else "FAIL ") + name,
          f"roc {recomputed['macro_roc']:.4f} vs {committed['roc_auc']:.4f}",
          f"pr {recomputed['macro_pr']:.4f} vs {committed['pr_auc']:.4f}",
          ("targets>2e-3: " + str(bad_targets)) if bad_targets else "")


# 13-target ep80 seed25: official + tuned
p = BIND / "migration_13target_v01" / "loso_preserve0.2_seed20260925.predictions.csv"
j = json.loads((BIND / "migration_13target_v01" / "loso_preserve0.2_seed20260925.json").read_text())
compare("13T ep80 seed25 / official", recompute(p, "official"), j["metrics"]["official"])
compare("13T ep80 seed25 / tuned(bce_retrieval)", recompute(p, "tuned"), j["metrics"]["bce_retrieval"])

# 13-target ep40 ensemble summary
p = BIND / "migration_13target_v01" / "summary_preserve0_3seed.predictions.csv"
j = json.loads((BIND / "migration_13target_v01" / "summary_preserve0_3seed.json").read_text())
compare("13T ep40 ensemble / tuned", recompute(p, "tuned"), j["metrics"]["tuned"])
compare("13T ep40 ensemble / ecfp4_logistic", recompute(p, "official"), j["metrics"]["ecfp4_logistic"] if "ecfp4_logistic" in j["metrics"] else None) if False else None
# ecfp predictions live in ecfp_loso.predictions.csv (separate file)

# 13-target ECFP has its own predictions file
p = BIND / "migration_13target_v01" / "ecfp_loso.predictions.csv"
j2 = json.loads((BIND / "migration_13target_v01" / "ecfp_loso.json").read_text())
with p.open(encoding="utf-8") as f:
    hdr = next(csv.reader(f))
compare("13T ECFP logistic", recompute(p, hdr[-1]), j2["metrics"])

# family augmentation seed25
p = BIND / "family_aug_v01" / "famaug_ep80_seed20260925.predictions.csv"
j = json.loads((BIND / "family_aug_v01" / "famaug_ep80_seed20260925.json").read_text())
compare("13T famaug seed25 / bce_retrieval", recompute(p, "tuned"), j["metrics"]["bce_retrieval"])

# 20-target ECFP baseline
p = BIND / "extended20_finetune_v01" / "ecfp20_loso.predictions.csv"
j = json.loads((BIND / "extended20_finetune_v01" / "ecfp20_loso.json").read_text())
first_col = None
with p.open(encoding="utf-8") as f:
    first_col = next(csv.reader(f))
score_col = "ecfp" if "ecfp" in first_col else first_col[-1]
compare(f"20T ECFP (col={score_col})", recompute(p, score_col), j["metrics"])

(OUT / "recompute_verification_v01.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
print("\nsaved recompute_verification_v01.json |",
      sum(1 for r in results if r["macro_roc_match"] and r["macro_pr_match"]),
      "/", len(results), "PASS")

