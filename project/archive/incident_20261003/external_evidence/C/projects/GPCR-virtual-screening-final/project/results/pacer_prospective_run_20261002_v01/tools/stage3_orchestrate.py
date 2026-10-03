"""Stage-3 orchestration: validate Model A, run Model B, apply the frozen dual-model rule."""
from __future__ import annotations

import csv, hashlib, json, sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
STAGE3 = RUN / "stage3"
OUT = RUN / "STAGE3_OUTPUTS"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "project"))

from project.drugclip_freeze import dual_model, model_a, model_b  # noqa: E402


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    archive = STAGE3 / "embeddings.npz"

    # --- Model A validation (fail closed) ---
    validation = model_a.validate_archive(archive)

    import numpy as np
    npz = np.load(archive, allow_pickle=False)
    # The frozen extractor labels molecules with their SMILES (sample["smi_name"]).
    # Resolve them to the deterministic PACER candidate ids via the frozen portfolio,
    # requiring a complete bijection.
    raw_ids = [str(x) for x in npz["molecule_ids"]]
    portfolio_rows = list(csv.DictReader(
        (REPO / "project/results/pacer_candidates_v01/predock_portfolio.csv").read_text(encoding="utf-8").splitlines()))
    smiles_to_cid = {r["canonical_smiles"]: r["candidate_id"] for r in portfolio_rows}
    if set(raw_ids) != set(smiles_to_cid):
        unmapped = sorted(set(raw_ids) - set(smiles_to_cid))[:5]
        raise SystemExit(f"STAGE3_ID_MAPPING_FAILED: {len(unmapped)} unmapped ids e.g. {unmapped}")
    molecule_ids = [smiles_to_cid[x] for x in raw_ids]
    pocket_ids = [str(x) for x in npz["pocket_ids"]]
    scores = npz["scores"].astype(np.float64)          # (pockets, molecules)

    pocket_index = {p: i for i, p in enumerate(pocket_ids)}
    m4 = [p for p in model_a.M4_POCKET_IDS]
    missing = [p for p in m4 if p not in pocket_index]
    if missing:
        raise SystemExit(f"STAGE3_BLOCKED_FROZEN_ASSET: missing M4 pockets {missing}")

    per_pocket = {p: scores[pocket_index[p]] for p in m4}
    state_mean = np.mean(np.stack([per_pocket[p] for p in m4]), axis=0)
    state_max = np.max(np.stack([per_pocket[p] for p in m4]), axis=0)

    # Model A decision-line score = official max-pooled cosine over the declared M4 pockets
    order = sorted(range(len(molecule_ids)), key=lambda i: (-state_max[i], molecule_ids[i]))
    old_rank = [0] * len(molecule_ids)
    for r, i in enumerate(order, 1):
        old_rank[i] = r

    model_a_rows = []
    portfolio = {r["candidate_id"]: r for r in portfolio_rows}
    for i, mid in enumerate(molecule_ids):
        model_a_rows.append({
            "candidate_id": mid,
            "canonical_smiles": portfolio[mid]["canonical_smiles"],
            **{f"model_a_{p}": float(per_pocket[p][i]) for p in m4},
            "model_a_state_mean": float(state_mean[i]),
            "model_a_state_max": float(state_max[i]),
            "model_a_old_rank": old_rank[i],
        })
    model_a_rows.sort(key=lambda r: r["model_a_old_rank"])
    with (OUT / "model_a_scores.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_a_rows[0].keys()))
        w.writeheader(); w.writerows(model_a_rows)

    # --- Model B (frozen projection-only) ---
    projections = model_b.default_projection_paths(REPO)
    b_report = model_b.score(archive, projections)
    model_b_csv = model_b.write_csv(b_report, OUT / "model_b_scores.csv")
    (OUT / "model_b_scores.json").write_text(json.dumps(b_report, indent=2) + "\n", encoding="utf-8")

    new_rank_by_id = {smiles_to_cid[mid]: int(r) for mid, r in zip(b_report["molecule_ids"], b_report["new_rank"])}
    model_b_score_by_id = {smiles_to_cid[mid]: float(s)
                           for mid, s in zip(b_report["molecule_ids"], b_report["ensemble_score"])}

    records = [{"candidate_id": r["candidate_id"], "old_rank": r["model_a_old_rank"],
                "new_rank": new_rank_by_id[r["candidate_id"]]} for r in model_a_rows]
    decision = dual_model.decide(records)
    (OUT / "STAGE3_DUAL_MODEL_DECISION_v01.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")

    old_top50 = sorted([r["candidate_id"] for r in records if r["old_rank"] <= 50])
    new_top25 = sorted([r["candidate_id"] for r in records if r["new_rank"] <= 25])
    # decision["dual_top_candidates"] is already ordered by (old_rank, new_rank)
    intersection = [r["candidate_id"] for r in decision["dual_top_candidates"]]

    shortlist_rows = []
    for r in model_a_rows:
        cid = r["candidate_id"]
        if cid in set(intersection):
            shortlist_rows.append({
                "shortlist_rank": intersection.index(cid) + 1,
                "candidate_id": cid,
                "canonical_smiles": r["canonical_smiles"],
                "model_a_state_max": r["model_a_state_max"],
                "model_a_old_rank": r["model_a_old_rank"],
                "model_b_score": model_b_score_by_id[cid],
                "model_b_new_rank": new_rank_by_id[cid],
                "selection_rule": "old_rank <= 50 AND new_rank <= 25",
            })
    with (OUT / "STAGE3_DUAL_MODEL_SHORTLIST_v01.csv").open("w", newline="", encoding="utf-8") as fh:
        if shortlist_rows:
            w = csv.DictWriter(fh, fieldnames=list(shortlist_rows[0].keys()))
            w.writeheader(); w.writerows(shortlist_rows)

    receipt = {
        "schema": "pacer.prospective.stage3_drugclip_receipt.v1",
        "stage": "STAGE3_DRUGCLIP",
        "status": "STAGE3_DRUGCLIP_COMPLETE",
        "run_id": "PACER-PROSPECTIVE-20261002",
        "stage2_input_molecule_count": len(molecule_ids),
        "model_a": {
            "validation": validation,
            "decision_line_score": "official max-pooled cosine over the three declared M4 pockets",
            "pockets": m4,
            "per_molecule": model_a_rows,
            "old_top50_count": len(old_top50),
            "old_top50": old_top50,
        },
        "model_b": {
            "deployed_seeds": b_report["seeds"],
            "projection_files": b_report["projection_files"],
            "pooling": b_report["pooling"],
            "new_top25_count": len(new_top25),
            "new_top25": new_top25,
            "per_molecule": [{"candidate_id": m, "model_b_score": model_b_score_by_id[m],
                              "model_b_new_rank": new_rank_by_id[m]} for m in molecule_ids],
        },
        "dual_model_decision": {**decision, "intersection_count": len(intersection),
                                "intersection": intersection},
        "dual_model_shortlist": shortlist_rows,
        "deterministic_candidate_ids": sorted(r["candidate_id"] for r in model_a_rows) ==
                                       [f"PACER{i:04d}" for i in range(1, len(molecule_ids) + 1)],
        "candidate_id_count": len({r["candidate_id"] for r in model_a_rows}),
        "outputs": {p.name: {"sha256": sha(p), "bytes": p.stat().st_size}
                    for p in sorted(OUT.glob("*"))},
        "forbidden_operations_performed": {
            "pocket_realign_fitting": False, "adapter_training": False, "finetuning": False,
            "optimizer_created": False, "backward_pass": False, "model_weight_update": False,
            "threshold_tuning": False, "score_fusion_invented": False},
        "claim_boundary": ("Frozen DrugCLIP apply-only binding-compatibility evidence and the frozen "
                           "dual-model prioritization rule only. Not a PAM classifier, not potency, "
                           "not experimental confirmation."),
    }
    (RUN / "STAGE3_DRUGCLIP_RECEIPT_v01.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in receipt.items()
                      if k not in ("model_a", "model_b", "dual_model_shortlist")}, indent=2)[:2500])


if __name__ == "__main__":
    main()
