#!/usr/bin/env python
"""Molecule-level PACER-DC data inventory and TRAINING_GATE decision.

Counts molecules, chemotypes, complete replicas and complete four-context windows
from the *actual* trajectory tree, joins the frozen manifest roles as labels, and
decides the training gate. Windows and frames are never counted as independent
samples.

    python audit_pacer_dc_training_gate.py \
        --manifest project/config/pacer_dc_pilot_md_manifest_v2.csv \
        --runs project/results/pacer_dc_production_v01 \
        --window-frames 100 \
        --output project/results/pacer_dc_four_context_v01/training_gate.json
"""
from __future__ import annotations

import argparse
import csv
import json
import struct
from pathlib import Path

CONTEXTS = ("candidate_probe", "candidate_no_probe", "probe_only", "apo")

# Frozen manifest roles -> functional labels. Anything not listed is unlabelled.
ROLE_LABELS = {
    "known_PAM_positive_control": {"pam_label": 1, "agonism_label": 0},
    "same_source_PAM_positive_control": {"pam_label": 1, "agonism_label": 0},
    "allosteric_agonist_specificity_control": {"pam_label": 0, "agonism_label": 1},
    "strict_experimental_inactive_binding_to_be_tested": {"pam_label": 0, "agonism_label": 0},
    "prospective_computational_hypothesis": {"pam_label": None, "agonism_label": None},
}


def frames(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open("rb") as fh:
        head = fh.read(16)
    if len(head) < 12:
        return 0
    if head[4:8] == b"CORD":
        return int(struct.unpack("<i", head[8:12])[0])
    if head[:4] == b"CORD":
        return int(struct.unpack("<i", head[4:8])[0])
    return 0


def system_dir(candidate: str, context: str) -> str:
    return context if context in ("probe_only", "apo") else f"{candidate}__{context}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path,
                    default=Path("project/config/pacer_dc_pilot_md_manifest_v2.csv"))
    ap.add_argument("--runs", type=Path, default=Path("project/results/pacer_dc_production_v01"))
    ap.add_argument("--window-frames", type=int, default=100)
    ap.add_argument("--output", type=Path,
                    default=Path("project/results/pacer_dc_four_context_v01/training_gate.json"))
    args = ap.parse_args()

    with args.manifest.open(encoding="utf-8") as fh:
        manifest = [r for r in csv.DictReader(fh)]

    molecules: dict[str, dict] = {}
    for row in manifest:
        cand, ctx = row["candidate_id"], row["context"]
        rep = row.get("replicate_id") or row.get("replica_id")
        if cand == "SHARED_CONTROL":
            continue
        m = molecules.setdefault(cand, {
            "candidate_id": cand, "role": row["role"],
            # The frozen manifest has NO chemotype column; Murcko scaffolds are not
            # computed here, so chemotype stays unknown rather than invented.
            "chemotype": row.get("chemotype") or None,
            "replica_frames": {},
        })
        d = args.runs / system_dir(cand, ctx) / f"replica_{int(row['paired_seed_group']):02d}"
        n = frames(d / "trajectory.dcd")
        m["replica_frames"].setdefault(rep, {})[ctx] = n

    # shared controls, counted once
    shared = {}
    for row in manifest:
        if row["candidate_id"] != "SHARED_CONTROL":
            continue
        rep = row.get("replicate_id") or row.get("replica_id")
        d = args.runs / row["context"] / f"replica_{int(row['paired_seed_group']):02d}"
        shared.setdefault(rep, {})[row["context"]] = frames(d / "trajectory.dcd")

    need = args.window_frames
    inventory = []
    for cand, m in sorted(molecules.items()):
        labels = ROLE_LABELS.get(m["role"], {"pam_label": None, "agonism_label": None})
        complete_reps, windows_total = 0, 0
        for rep, got in sorted(m["replica_frames"].items()):
            sh = shared.get(rep, {})
            # probe_only and apo are SHARED controls: they are not rows of this
            # molecule in the manifest, so completeness must consult `shared`.
            all_ctx = all(got.get(c, 0) >= need for c in ("candidate_probe", "candidate_no_probe")) \
                and all(sh.get(c, 0) >= need for c in ("probe_only", "apo"))
            own = min([got.get("candidate_probe", 0), got.get("candidate_no_probe", 0)] or [0])
            sh_min = min([sh.get("probe_only", 0), sh.get("apo", 0)] or [0])
            avail = min(own, sh_min)
            if all_ctx:
                complete_reps += 1
            windows_total += avail // need
        inventory.append({
            "candidate_id": cand, "role": m["role"], "chemotype": m["chemotype"],
            "pam_label": labels["pam_label"], "agonism_label": labels["agonism_label"],
            "replicas_with_all_four_contexts": complete_reps,
            "complete_four_context_windows": windows_total,
            "frames_per_context": {r: m["replica_frames"][r] for r in sorted(m["replica_frames"])},
        })

    labelled = [m for m in inventory if m["pam_label"] is not None]
    chemotypes = {m["chemotype"] for m in inventory if m["chemotype"]}
    pam_pos = sum(1 for m in labelled if m["pam_label"] == 1)
    pam_neg = sum(1 for m in labelled if m["pam_label"] == 0)
    ago_pos = sum(1 for m in labelled if m["agonism_label"] == 1)
    with_windows = [m for m in inventory if m["complete_four_context_windows"] > 0]

    # Gate: >=6 independent training molecules, both classes >=2 in train, AND a
    # leak-free molecule/chemotype-aware train/val/test. Counted per molecule.
    reasons = []
    if len(inventory) < 6:
        reasons.append(f"only {len(inventory)} independent candidate molecules (need >=6)")
    if len(with_windows) < 6:
        reasons.append(f"only {len(with_windows)} molecules have any complete four-context window")
    if pam_pos < 2 or pam_neg < 2:
        reasons.append(f"train PAM classes need >=2 each: pam_label=1 ->{pam_pos}, =0 ->{pam_neg}")
    if ago_pos < 2:
        reasons.append(f"agonism_label=1 molecules: {ago_pos} (need >=2)")
    if len(inventory) + 0 < 7:
        reasons.append("cannot form train/val/test at molecule level with this few molecules")
    gate = "OPEN" if not reasons else "CLOSED"

    report = {
        "manifest": str(args.manifest),
        "runs_root": str(args.runs),
        "window_frames": need,
        "independent_candidates": len(inventory),
        "independent_chemotypes": len(chemotypes),
        "molecules_with_complete_four_context_windows": len(with_windows),
        "total_complete_four_context_windows": sum(m["complete_four_context_windows"] for m in inventory),
        "labelled_molecules": len(labelled),
        "pam_label_1": pam_pos, "pam_label_0": pam_neg, "agonism_label_1": ago_pos,
        "shared_control_frames": shared,
        "windows_are_not_independent_samples": True,
        "inventory": inventory,
        "gate_reasons": reasons,
        "TRAINING_GATE": gate,
        "claim_boundary": (
            "Inventory and gate only. No training, no AUC, no PAM claim. Windows and "
            "replica frames are augmentation, never independent molecules."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"{'candidate':<34}{'role':<46}{'reps4ctx':<10}{'windows':<9}pam/ago")
    print("-" * 110)
    for m in inventory:
        print(f"{m['candidate_id']:<34}{m['role']:<46}{m['replicas_with_all_four_contexts']:<10}"
              f"{m['complete_four_context_windows']:<9}{m['pam_label']}/{m['agonism_label']}")
    print("-" * 110)
    print(f"independent candidates          : {report['independent_candidates']}")
    print(f"independent chemotypes          : {report['independent_chemotypes']}")
    print(f"molecules with 4-context windows: {report['molecules_with_complete_four_context_windows']}")
    print(f"labelled molecules              : {report['labelled_molecules']} "
          f"(pam+= {pam_pos}, pam-= {pam_neg}, ago+= {ago_pos})")
    for r in reasons:
        print(f"  BLOCKER: {r}")
    print(f"\nTRAINING_GATE = {gate}")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
