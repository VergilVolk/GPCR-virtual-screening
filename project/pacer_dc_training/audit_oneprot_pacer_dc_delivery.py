#!/usr/bin/env python
"""Audit the OneProt-MD to PACER-DC handoff without training or inventing evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "candidate_id", "chemotype", "replicate_id", "window_id", "context",
    "pam_label", "agonism_label", "split",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_head(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def audit_features(path: Path, contexts: set[str]) -> dict:
    result = {"path": str(path), "present": path.is_file(), "schema_valid": False}
    if not path.is_file():
        return result
    frame = pd.read_csv(path)
    missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
    feature_cols = sorted(c for c in frame.columns if c.startswith("f_"))
    result.update({"rows": int(len(frame)), "missing_columns": missing,
                   "embedding_dimension": len(feature_cols)})
    if missing or not feature_cols:
        return result
    values = frame[feature_cols].to_numpy(float)
    result["features_finite"] = bool(np.isfinite(values).all())
    result["features_nonconstant"] = bool(np.nanstd(values, axis=0).max() > 0)

    keys = ["candidate_id", "replicate_id", "window_id"]
    context_sets = frame.groupby(keys).context.agg(lambda x: set(map(str, x)))
    result["four_context_unit_count"] = int((context_sets == contexts).sum())
    result["incomplete_unit_count"] = int((context_sets != contexts).sum())

    candidate_leak = frame.groupby("candidate_id").split.nunique()
    chemotype_leak = frame.groupby("chemotype").split.nunique()
    result["candidate_split_leakage"] = sorted(map(str, candidate_leak[candidate_leak > 1].index))
    result["chemotype_split_leakage"] = sorted(map(str, chemotype_leak[chemotype_leak > 1].index))
    molecules = frame.drop_duplicates("candidate_id")
    train = molecules[molecules.split.eq("train")]
    result["candidate_counts"] = {str(k): int(v) for k, v in molecules.split.value_counts().items()}
    result["train_pam_class_counts"] = {
        str(int(k)): int(v) for k, v in train.pam_label.value_counts().items()
    }
    result["train_agonism_class_counts"] = {
        str(int(k)): int(v) for k, v in train.agonism_label.value_counts().items()
    }
    class_gate = all(
        set(train[label].dropna().astype(int)) == {0, 1}
        and int(train[label].value_counts().min()) >= 2
        for label in ["pam_label", "agonism_label"]
    )
    result["schema_valid"] = bool(
        result["features_finite"] and result["features_nonconstant"]
        and result["incomplete_unit_count"] == 0
        and not result["candidate_split_leakage"]
        and not result["chemotype_split_leakage"]
    )
    result["training_data_gate"] = bool(
        result["schema_valid"] and set(molecules.split) >= {"train", "val", "test"}
        and train.candidate_id.nunique() >= 6 and class_gate
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oneprot-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--config", type=Path,
                        default=Path("project/config/oneprot_pacer_dc_delivery_v01.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    head = git_head(args.oneprot_root)
    checkpoint_present = args.checkpoint.is_file()
    mdgen = args.oneprot_root / cfg["required_submodule"]
    feature_audit = audit_features(args.features, set(cfg["required_contexts"]))
    audit = {
        "oneprot_root": str(args.oneprot_root),
        "oneprot_commit": head,
        "oneprot_commit_matches": head == cfg["expected_oneprot_commit"],
        "mdgen_submodule_present": mdgen.is_dir() and any(mdgen.iterdir()),
        "checkpoint": str(args.checkpoint),
        "checkpoint_present": checkpoint_present,
        "checkpoint_sha256": sha256(args.checkpoint) if checkpoint_present else None,
        "features": feature_audit,
        "technical_delivery_gate": bool(
            head == cfg["expected_oneprot_commit"] and mdgen.is_dir()
            and any(mdgen.iterdir()) and checkpoint_present and feature_audit.get("schema_valid", False)
        ),
        "adapter_training_gate": bool(feature_audit.get("training_data_gate", False)),
        "claim_boundary": cfg["claim_boundary"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if not audit["technical_delivery_gate"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
