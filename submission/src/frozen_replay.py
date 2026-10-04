from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from module4_four_context import export_evidence

TRACK = "Track 3: small-molecule virtual screening for ion channels and GPCR drug targets"
RUN_VERSION = "pacer-m4-frozen-screen-2026-10-04"
MODEL_VERSION = "DrugCLIP-2023 GPCR-adapted 3-seed ensemble; Glide/Vina PDB-BEmin-BEavg fusion; four-context dynamics"
CHANNELS = ["Glide_PDB", "Glide_BEmin", "Glide_BEavg", "Vina_PDB", "Vina_BEmin", "Vina_BEavg"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def verify_manifest(data_dir: Path) -> dict:
    manifest = json.loads((data_dir / "manifest.json").read_text(encoding="utf-8"))
    for name, record in manifest["files"].items():
        path = data_dir / name
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"Frozen data hash mismatch: {name}")
    return manifest


def replay_frozen(data_dir: Path, work_dir: Path, output: Path) -> pd.DataFrame:
    verify_manifest(data_dir)
    generated = pd.read_csv(data_dir / "generated_candidates_2605.csv")
    audit = pd.read_csv(data_dir / "candidate_audit_2592.csv")
    binding = pd.read_csv(data_dir / "drugclip_rerun_top200.csv")
    docking = pd.read_csv(data_dir / "six_channel_top5_plus_anchors.csv")

    if len(generated) != 2605 or generated.smiles.map(canonical).nunique() != 2605:
        raise ValueError("Generated library must contain 2,605 unique valid molecular graphs")
    if len(audit) != 2592 or audit.canonical_smiles.map(canonical).nunique() != 2592:
        raise ValueError("Candidate audit must contain 2,592 unique valid molecular graphs")
    if len(binding) != 200 or binding.smiles.map(canonical).nunique() != 200:
        raise ValueError("DrugCLIP rerun table must contain 200 unique valid molecular graphs")
    if len(docking) != 8 or docking.canonical_smiles.map(canonical).nunique() != 8:
        raise ValueError("Six-channel table must contain eight unique molecular graphs")

    binding = binding.copy()
    binding["canonical_graph"] = binding.smiles.map(canonical)
    docking = docking.copy()
    docking["canonical_graph"] = docking.canonical_smiles.map(canonical)
    joined = docking.merge(
        binding[["gen_id", "canonical_graph", "score", "rankpct", "final_rank"]],
        left_on=["candidate_id", "canonical_graph"],
        right_on=["gen_id", "canonical_graph"],
        how="left",
        validate="one_to_one",
    )
    if joined.gen_id.isna().any():
        raise ValueError("A docked candidate is absent from the frozen DrugCLIP top-200 table")

    rank_columns = [f"{channel}_rankpct" for channel in CHANNELS]
    raw_columns = [f"{channel}_raw" for channel in CHANNELS]
    if joined[rank_columns + raw_columns].isna().any().any():
        raise ValueError("All eight candidates require six finite docking channels")
    joined["six_channel_score_recomputed"] = joined[rank_columns].mean(axis=1)
    if not np.allclose(joined["PACER_XR"], joined["six_channel_score_recomputed"], atol=1e-12):
        raise ValueError("Six-channel score does not match the frozen mean-rank definition")

    mapping = pd.read_csv(data_dir / "candidate_id_mapping.csv")
    joined = joined.drop(columns=["legacy_id"], errors="ignore").merge(
        mapping, on="candidate_id", how="left", validate="one_to_one")
    stage4 = export_evidence(data_dir / "four_context_results.json", work_dir / "module4_function.csv")
    joined = joined.merge(stage4, on="legacy_id", how="left", validate="many_to_one")
    joined["four_context_status"] = joined["four_context_status"].fillna("not_run")
    joined["functional_prediction"] = joined["functional_prediction"].fillna("not_evaluated")
    joined["experimental_pam_validation"] = joined["experimental_pam_validation"].fillna("not_performed")

    joined = joined.sort_values("PACER_XR_comparison_rank", kind="stable").reset_index(drop=True)
    result = pd.DataFrame({
        "candidate_rank": np.arange(1, len(joined) + 1),
        "candidate_id": joined["candidate_id"],
        "historical_id": joined["legacy_id"].fillna(""),
        "track": TRACK,
        "smiles": joined["canonical_smiles"],
        "candidate_group": joined["group"],
        "drugclip_score": joined["score"],
        "drugclip_rank_within_2605": joined["final_rank"].astype(int),
        "six_channel_score": joined["PACER_XR"],
        "six_channel_rank_within_8": joined["PACER_XR_comparison_rank"].astype(int),
        "glide_pdb_score": joined["Glide_PDB_raw"],
        "glide_bemin_score": joined["Glide_BEmin_raw"],
        "glide_beavg_score": joined["Glide_BEavg_raw"],
        "vina_pdb_kcal_mol": joined["Vina_PDB_raw"],
        "vina_bemin_kcal_mol": joined["Vina_BEmin_raw"],
        "vina_beavg_kcal_mol": joined["Vina_BEavg_raw"],
        "four_context_status": joined["four_context_status"],
        "delta_int_r1_r3_direction_cosine": joined.get("delta_int_r1_r3_direction_cosine"),
        "functional_prediction": joined["functional_prediction"],
        "experimental_pam_validation": joined["experimental_pam_validation"],
        "structure_file": "not_redistributed",
        "model_version": MODEL_VERSION,
        "run_version": RUN_VERSION,
    })
    result["remarks"] = np.where(
        result.four_context_status.eq("complete_unresolved"),
        "Four-context dynamics completed; functional identity remains unresolved and requires experiment.",
        "Binding and six-channel docking evidence only; four-context dynamics not run.",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    return result
