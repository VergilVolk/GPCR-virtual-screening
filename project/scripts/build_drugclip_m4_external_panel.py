#!/usr/bin/env python3
"""Build one deduplicated DrugCLIP input panel for frozen M4 external tests.

The script preserves each source endpoint separately.  It never pools assay
labels and independently recomputes exact overlap with the M4 training table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem


def canonicalize(smiles: object) -> str | None:
    mol = Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) if mol else None


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def add_rows(rows: list[dict], path: Path, dataset: str) -> None:
    table = pd.read_csv(path)
    for _, row in table.iterrows():
        smiles_column = "canonical_smiles" if "canonical_smiles" in table else "smiles"
        smiles = canonicalize(row.get(smiles_column))
        if not smiles:
            continue
        common = {
            "dataset": dataset,
            "record_id": str(row.get("compound_id", row.get("example_id", len(rows)))),
            "canonical_smiles": smiles,
            "label_binary": None,
            "potency_value": None,
            "ordinal_class": None,
            "eligible": True,
        }
        if dataset == "monash_allostery":
            common["endpoint"] = "reported_PAM_vs_inactive"
            common["label_binary"] = int(row["pam_label"])
        elif dataset == "acadia_functional":
            common["endpoint"] = "hM4_PAM_RE_ge_50pct"
            common["label_binary"] = int(row["pam_re_ge_50"])
            common["eligible"] = str(row.get("structure_status", "")) == "high_confidence"
        elif dataset == "suven_potency":
            common["endpoint"] = "ACh_EC20_CRE_luc_pEC50"
            common["potency_value"] = pd.to_numeric(row.get("cre_luc_pEC50"), errors="coerce")
            common["eligible"] = truthy(row.get("eligible_external", True))
        elif dataset == "vu6025733_potency":
            common["endpoint"] = "ACh_EC20_calcium_pEC50"
            common["potency_value"] = pd.to_numeric(row.get("pEC50"), errors="coerce")
            common["eligible"] = truthy(row.get("eligible_zero_shot", True))
        elif dataset == "us20260055116_ordinal":
            common["endpoint"] = "human_M4_pERK_patent_bin"
            category = str(row.get("human_m4_perk", "")).strip().upper()
            common["ordinal_class"] = {"A": 4, "B": 3, "C": 2, "D": 1}.get(category)
            common["eligible"] = truthy(row.get("high_confidence_named_structure", True))
        rows.append(common)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-root", type=Path, required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    specs = [
        ("monash_allostery", "external_monash_ly2033298/external_allostery.csv"),
        ("acadia_functional", "external_acadia_2025/external_functional_benchmark.csv"),
        ("suven_potency", "external_suven_2025/external_patent_potency.csv"),
        ("vu6025733_potency", "external_2026_vu6025733/external_potency.csv"),
        ("us20260055116_ordinal", "external_us20260055116/named_structure_functional_subset.csv"),
    ]
    rows: list[dict] = []
    for dataset, relative in specs:
        add_rows(rows, args.benchmark_root / relative, dataset)
    manifest = pd.DataFrame(rows)
    training = pd.read_csv(args.training)
    training_smiles = {canonicalize(v) for v in training["canonical_smiles"]}
    training_smiles.discard(None)
    manifest["exact_training_overlap_computed"] = manifest.canonical_smiles.isin(training_smiles)
    unique = list(dict.fromkeys(manifest.canonical_smiles))
    id_map = {s: "EXT_" + hashlib.sha1(s.encode()).hexdigest()[:12] for s in unique}
    manifest["external_molecule_id"] = manifest.canonical_smiles.map(id_map)
    molecules = pd.DataFrame({
        "external_molecule_id": [id_map[s] for s in unique],
        "canonical_smiles": unique,
    })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output_dir / "endpoint_manifest.csv", index=False)
    molecules.to_csv(args.output_dir / "molecules.csv", index=False)
    counts = {}
    for dataset, group in manifest.groupby("dataset"):
        counts[dataset] = {
            "rows": int(len(group)),
            "eligible": int(group.eligible.sum()),
            "exact_training_overlap": int(group.exact_training_overlap_computed.sum()),
            "unique_molecules": int(group.canonical_smiles.nunique()),
        }
    audit = {
        "training_table": str(args.training),
        "endpoint_rows": int(len(manifest)),
        "unique_external_molecules": int(len(molecules)),
        "datasets": counts,
        "pooling_rule": "No cross-assay label pooling; every endpoint is evaluated separately.",
        "external_rule": "Exact training overlaps are encoded but excluded from zero-shot metrics.",
    }
    (args.output_dir / "panel.audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
