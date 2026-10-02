#!/usr/bin/env python3
"""Build the endpoint-stratified PACER-M4 external functional benchmark.

The statistical unit remains a unique molecule within an assay endpoint.  The
script deliberately does not pool raw values from different probes/readouts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


ORDINAL = {"A": 4.0, "B": 3.0, "C": 2.0, "D": 1.0}


def canonicalize(value: object) -> str | None:
    mol = Chem.MolFromSmiles(str(value)) if pd.notna(value) else None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) if mol else None


def molecule_id(smiles: str) -> str:
    return "EXT_" + hashlib.sha1(smiles.encode()).hexdigest()[:12]


def add(rows: list[dict], *, dataset: str, record_id: object, smiles: object,
        endpoint: str, task: str, value: object, role: str = "primary_function",
        exact_overlap: bool = False, source: str, probe: str, readout: str) -> None:
    smi = canonicalize(smiles)
    val = pd.to_numeric(value, errors="coerce")
    if smi is None or pd.isna(val):
        return
    rows.append({
        "dataset": dataset, "record_id": str(record_id), "canonical_smiles": smi,
        "external_molecule_id": molecule_id(smi), "endpoint": endpoint,
        "task_type": task, "value": float(val), "endpoint_role": role,
        "exact_training_overlap_reported": bool(exact_overlap), "source": source,
        "orthosteric_probe": probe, "readout": readout,
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-root", type=Path, required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.benchmark_root
    rows: list[dict] = []

    monash = pd.read_csv(root / "external_monash_ly2033298/external_allostery.csv")
    for _, r in monash.iterrows():
        add(rows, dataset="monash_2015", record_id=r.get("compound_id", r.name),
            smiles=r.canonical_smiles, endpoint="reported_PAM_vs_inactive",
            task="binary", value=r.pam_label, source="Monash LY2033298 series",
            probe="ACh", readout="pERK functional allostery")

    acadia = pd.read_csv(root / "external_acadia_2025/external_functional_benchmark.csv")
    for _, r in acadia.iterrows():
        common = dict(dataset="acadia_2025", record_id=r.compound_id,
                      smiles=r.canonical_smiles, exact_overlap=r.exact_training_overlap,
                      source="WO2025122811A1", probe="ACh EC20",
                      readout="CHO-K1/M4/Galpha15 calcium flux")
        add(rows, endpoint="PAM_RE_ge_50pct", task="binary", value=r.pam_re_ge_50, **common)
        add(rows, endpoint="PAM_relative_efficacy_pct", task="continuous",
            value=r.hM4_PAM_RE_pct, **common)
        if str(r.hM4_PAM_EC50_qualifier).strip() == "=":
            add(rows, endpoint="PAM_pEC50", task="continuous", value=r.pEC50, **common)
        add(rows, endpoint="intrinsic_agonist_RE_pct", task="continuous",
            value=r.hM4_agonist_RE_pct, role="liability", **common)
        add(rows, endpoint="intrinsic_agonist_RE_ge_50pct", task="binary",
            value=r.intrinsic_agonist_re_ge_50, role="liability", **common)

    suven = pd.read_csv(root / "external_suven_2025/external_patent_potency.csv")
    for _, r in suven.iterrows():
        common = dict(dataset="suven_2025", record_id=r.compound_id, smiles=r.smiles,
                      exact_overlap=r.exact_training_overlap, source="WO2025099660A1",
                      probe="ACh EC20", readout="human M4 CRE-luciferase")
        add(rows, endpoint="CRE_luc_PAM_pEC50", task="continuous",
            value=r.cre_luc_pEC50, **common)
        add(rows, endpoint="GloSensor_PAM_pEC50", task="continuous",
            value=r.glosensor_pEC50, readout="human M4 GloSensor cAMP",
            **{k: v for k, v in common.items() if k != "readout"})
        if pd.notna(r.m4_m2_selectivity_fold_lower_bound) and float(r.m4_m2_selectivity_fold_lower_bound) > 0:
            add(rows, endpoint="log10_M4_over_M2_selectivity", task="continuous",
                value=np.log10(float(r.m4_m2_selectivity_fold_lower_bound)),
                role="selectivity", readout="M4/M2 CRE-luciferase",
                **{k: v for k, v in common.items() if k != "readout"})

    vu = pd.read_csv(root / "external_2026_vu6025733/external_potency.csv")
    for _, r in vu.iterrows():
        add(rows, dataset="vu6025733_2026", record_id=r.compound_id,
            smiles=r.canonical_smiles, endpoint="calcium_PAM_pEC50", task="continuous",
            value=r.pEC50, exact_overlap=r.exact_training_overlap,
            source=str(r.source_doi), probe="ACh EC20",
            readout="hM4/Gqi5-CHO calcium mobilization")

    patent = pd.read_csv(root / "external_us20260055116/named_structure_functional_subset.csv")
    contexts = [
        ("human_pERK_PAM_bin", "human_m4_perk", "human M4 pERK"),
        ("rat_pERK_PAM_bin", "rat_m4_perk", "rat M4 pERK"),
        ("human_GTPgammaS_PAM_bin", "human_m4_gtpgs", "human M4 GTPgammaS"),
    ]
    for _, r in patent.iterrows():
        for endpoint, column, readout in contexts:
            label = ORDINAL.get(str(r[column]).strip().upper())
            add(rows, dataset="us20260055116_2026", record_id=r.compound_id,
                smiles=r.canonical_smiles, endpoint=endpoint, task="ordinal", value=label,
                exact_overlap=r.exact_train_overlap, source="US20260055116A1",
                probe="ACh", readout=readout)

    table = pd.DataFrame(rows)
    training = pd.read_csv(args.training)
    train_smiles = [canonicalize(v) for v in training.canonical_smiles]
    train_smiles = [v for v in train_smiles if v]
    train_set = set(train_smiles)
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    train_fps = [fpgen.GetFingerprint(Chem.MolFromSmiles(v)) for v in train_smiles]
    unique = table[["external_molecule_id", "canonical_smiles"]].drop_duplicates().copy()
    unique["exact_training_overlap_computed"] = unique.canonical_smiles.isin(train_set)
    unique["murcko_scaffold"] = unique.canonical_smiles.map(
        lambda s: MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s)))
    unique["max_train_tanimoto_ecfp4"] = unique.canonical_smiles.map(
        lambda s: max(DataStructs.BulkTanimotoSimilarity(
            fpgen.GetFingerprint(Chem.MolFromSmiles(s)), train_fps)))
    table = table.merge(unique, on=["external_molecule_id", "canonical_smiles"], how="left")
    table["eligible_zero_shot"] = ~(
        table.exact_training_overlap_reported | table.exact_training_overlap_computed)

    # Duplicate salts/rows within an endpoint are one statistical unit.
    keys = ["dataset", "endpoint", "external_molecule_id"]
    dedup = table.sort_values(keys).drop_duplicates(keys, keep="first")
    summary = []
    for (dataset, endpoint), g in dedup.groupby(["dataset", "endpoint"], sort=False):
        eligible = g[g.eligible_zero_shot]
        task = str(g.task_type.iloc[0])
        counts = eligible.value.value_counts().to_dict() if task in {"binary", "ordinal"} else {}
        primary = str(g.endpoint_role.iloc[0]) == "primary_function"
        confirmatory = primary and len(eligible) >= 15 and (
            task != "binary" or (eligible.value.eq(0).sum() >= 10 and eligible.value.eq(1).sum() >= 10))
        summary.append({
            "dataset": dataset, "endpoint": endpoint, "task_type": task,
            "endpoint_role": g.endpoint_role.iloc[0], "n_rows": int(len(g)),
            "n_zero_shot": int(len(eligible)), "n_scaffolds": int(eligible.murcko_scaffold.nunique()),
            "class_counts": json.dumps({str(k): int(v) for k, v in counts.items()}),
            "confirmatory_macro_eligible": bool(confirmatory),
        })
    summary_df = pd.DataFrame(summary)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "endpoint_records.csv", index=False)
    unique.to_csv(args.output_dir / "molecules.csv", index=False)
    summary_df.to_csv(args.output_dir / "endpoint_summary.csv", index=False)
    audit = {
        "version": "PACER-M4 external functional benchmark v02",
        "endpoint_records": int(len(table)),
        "unique_molecules": int(unique.external_molecule_id.nunique()),
        "unique_scaffolds": int(unique.murcko_scaffold.nunique()),
        "zero_shot_endpoint_units": int(dedup.eligible_zero_shot.sum()),
        "endpoints": int(summary_df.shape[0]),
        "confirmatory_macro_endpoints": summary_df.loc[
            summary_df.confirmatory_macro_eligible, "endpoint"].tolist(),
        "pooling_rule": "Raw labels are never pooled across assay endpoints.",
        "independence_rule": "A molecule is counted once per endpoint; repeated endpoints on one molecule are not independent molecules.",
        "claim_boundary": "Retrospective public functional benchmark; not prospective or wet-lab confirmation of a new PAM.",
    }
    (args.output_dir / "benchmark.audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
