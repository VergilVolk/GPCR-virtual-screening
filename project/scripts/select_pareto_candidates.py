#!/usr/bin/env python3
"""Select a diverse, non-dominated PACER-M4 candidate portfolio.

The Pareto objectives are pocket coverage, QED, bounded novelty and a weak
inactive-risk safety flag. Vina affinity and unvalidated potency estimates are
reported for provenance but deliberately excluded from selection.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDConfig, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_BASE = PROJECT / "results" / "pacer_candidates_v01"
DEFAULT_KNOWN = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
PARETO_COLUMNS = ["pocket_residue_coverage", "QED", "novelty", "inactive_safety_ref"]
OUTPUT_COLUMNS = [
    "final_rank", "candidate_id", "canonical_smiles", "murcko_scaffold",
    "portfolio_role", "generator", "claim_level", "pareto_front",
    "diversity_cluster", "nearest_known_id", "nearest_known_pEC50",
    "max_tanimoto", "n_neighbors_045", "local_knn_pEC50", "local_neighbor_sd",
    "pocket_residue_coverage", "unique_pocket_residue_contacts",
    "contacted_residues", "vina_affinity", "QED", "SA_score", "MW", "logP",
    "TPSA", "PAINS", "druglike", "strict_inactive_risk_ref",
]

sys.path.append(str(Path(RDConfig.RDContribDir) / "SA_Score"))
import sascorer  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def fingerprint(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return AllChem.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=2048)


def pareto_mask(values: np.ndarray) -> np.ndarray:
    """Return the non-dominated rows when every column is maximized."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("Pareto input must be a two-dimensional array")
    keep = np.ones(len(values), dtype=bool)
    for index, row in enumerate(values):
        weakly_better = np.all(values >= row, axis=1)
        strictly_better = np.any(values > row, axis=1)
        if np.any(weakly_better & strictly_better):
            keep[index] = False
    return keep


def cluster_fingerprints(fingerprints: Sequence, distance_cutoff: float = 0.58):
    """Cluster fingerprints with the historical Butina distance threshold."""
    if not fingerprints:
        return tuple()
    distances: list[float] = []
    for index in range(1, len(fingerprints)):
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[index], fingerprints[:index]
        )
        distances.extend(1.0 - value for value in similarities)
    return Butina.ClusterData(
        distances, len(fingerprints), distance_cutoff, isDistData=True
    )


def assign_pareto_fronts(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"Missing Pareto columns: {missing}")
    output = frame.copy()
    output["pareto_front"] = -1
    front = 0
    while (output.pareto_front < 0).any():
        remaining = output.index[output.pareto_front < 0]
        mask = pareto_mask(output.loc[remaining, columns].to_numpy(float))
        output.loc[remaining[mask], "pareto_front"] = front
        front += 1
    return output


def choose_diverse(frame: pd.DataFrame, count: int) -> pd.DataFrame:
    """Choose one candidate per ECFP cluster before filling remaining slots."""
    if frame.empty or count <= 0:
        return frame.iloc[0:0].copy()
    output = frame.reset_index(drop=True).copy()
    fingerprints = [fingerprint(smiles) for smiles in output.canonical_smiles]
    clusters = cluster_fingerprints(fingerprints)
    cluster_by_row = {
        row_index: cluster_index
        for cluster_index, cluster in enumerate(clusters)
        for row_index in cluster
    }
    output["diversity_cluster"] = [cluster_by_row[index] for index in range(len(output))]

    primary_order = output.sort_values(
        ["pareto_front", "pocket_residue_coverage", "QED", "SA_score"],
        ascending=[True, False, False, True],
        kind="stable",
    )
    selected_indices: list[int] = []
    used_clusters: set[int] = set()
    for index, row in primary_order.iterrows():
        cluster = int(row.diversity_cluster)
        if cluster not in used_clusters:
            selected_indices.append(index)
            used_clusters.add(cluster)
        if len(selected_indices) >= count:
            break

    if len(selected_indices) < count:
        secondary_order = output.sort_values(
            ["pareto_front", "pocket_residue_coverage", "QED"],
            ascending=[True, False, False],
            kind="stable",
        )
        selected = set(selected_indices)
        for index in secondary_order.index:
            if index not in selected:
                selected_indices.append(index)
                selected.add(index)
            if len(selected_indices) >= count:
                break
    return output.loc[selected_indices].reset_index(drop=True)


def _murcko_smiles(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(molecule))


def select_portfolio(
    predock_path: Path,
    docking_features_path: Path,
    known_potency_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    candidates = pd.read_csv(predock_path)
    structures = pd.read_csv(docking_features_path)
    if "error" in structures:
        structures = structures[structures.error.fillna("") == ""]
    if "seed" not in structures:
        raise ValueError("Docking feature table must contain a seed column")
    structures = structures.sort_values("seed", kind="stable").drop_duplicates("candidate_id")
    candidates = candidates.merge(
        structures, on="candidate_id", validate="one_to_one", suffixes=("", "_docking")
    )
    if candidates.empty:
        raise ValueError("No candidates remain after joining docking features")

    candidates["SA_score"] = [
        float(sascorer.calculateScore(Chem.MolFromSmiles(smiles)))
        for smiles in candidates.canonical_smiles
    ]
    candidates["novelty"] = 1.0 - candidates.max_tanimoto
    candidates["inactive_safety_ref"] = 1.0 - candidates.strict_inactive_risk_ref
    candidates = assign_pareto_fronts(candidates, PARETO_COLUMNS)

    allocations = [
        ("local_exploitation", candidates[(candidates.domain_tier == "local") & (candidates.max_tanimoto >= 0.65)], 12),
        ("local_diversification", candidates[(candidates.domain_tier == "local") & (candidates.max_tanimoto < 0.65)], 6),
        ("exploratory_hypothesis", candidates[candidates.domain_tier == "exploratory"], 6),
    ]
    selected_parts: list[pd.DataFrame] = []
    for role, subset, count in allocations:
        selected = choose_diverse(subset, count)
        selected["portfolio_role"] = role
        selected_parts.append(selected)
    final = pd.concat(selected_parts, ignore_index=True).drop_duplicates("candidate_id")
    if final.empty:
        raise ValueError("No candidates satisfied the portfolio allocation rules")
    final.insert(0, "final_rank", np.arange(1, len(final) + 1))
    final["claim_level"] = "computational PAM hypothesis; requires functional ternary-complex assay"
    final["murcko_scaffold"] = [_murcko_smiles(smiles) for smiles in final.canonical_smiles]

    known = pd.read_csv(known_potency_path).set_index("canonical_molecule_id")
    final["nearest_known_pEC50"] = [
        known.at[identifier, "pEC50"] if identifier in known.index else np.nan
        for identifier in final.nearest_known_id
    ]
    missing_output = [column for column in OUTPUT_COLUMNS if column not in final]
    if missing_output:
        raise ValueError(f"Missing output columns after selection: {missing_output}")

    output_dir.mkdir(parents=True, exist_ok=True)
    final[OUTPUT_COLUMNS].to_csv(output_dir / "final_candidate_hypotheses.csv", index=False)
    candidates.to_csv(output_dir / "all_predock_with_structure.csv", index=False)
    audit: dict[str, object] = {
        "n_docked": int(len(candidates)),
        "n_final": int(len(final)),
        "roles": {str(key): int(value) for key, value in final.portfolio_role.value_counts().items()},
        "unique_ecfp_clusters": int(final.diversity_cluster.nunique()),
        "unique_murcko_scaffolds": int(final.murcko_scaffold.nunique()),
        "similarity_range": [float(final.max_tanimoto.min()), float(final.max_tanimoto.max())],
        "pocket_coverage_range": [
            float(final.pocket_residue_coverage.min()),
            float(final.pocket_residue_coverage.max()),
        ],
        "policy": (
            "Pareto objectives exclude Vina affinity and potency predictions. "
            "Local kNN and inactive classifier values are references only. "
            "Candidates require a functional PAM assay."
        ),
    }
    (output_dir / "selection_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predock", type=Path, default=DEFAULT_BASE / "predock_portfolio.csv")
    parser.add_argument(
        "--docking-features", type=Path, default=DEFAULT_BASE / "7trs_docking" / "features.csv"
    )
    parser.add_argument("--known-potency", type=Path, default=DEFAULT_KNOWN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_BASE / "final")
    args = parser.parse_args()
    audit = select_portfolio(
        args.predock, args.docking_features, args.known_potency, args.output_dir
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
