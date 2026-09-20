# -*- coding: utf-8 -*-
"""Build and audit M4-PAM-Benchmark v1.0.

This script does not train a model. It freezes molecule-level benchmark tables,
group labels and leakage audits before PACER-M4 development starts.

Outputs (project/data/benchmarks/m4_pam_v1):
  potency_molecules.csv       assay-harmonized pEC50 regression view
  pam_vs_inactive.csv         strict functional classification view
  threshold_10uM.csv          potency-threshold classification view
  source_metadata.csv         source/year metadata used by the benchmark
  audit_report.json           coverage, overlap and split diagnostics

Optional PubMed enrichment fills missing publication years in the benchmark
copy only; the frozen handoff dataset is never modified.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import time
import urllib.parse
import urllib.request
import urllib.error
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import StratifiedGroupKFold

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(__file__).resolve().parents[1]
REPO = PROJECT.parent
FROZEN = REPO / "CHRM4_PAM_Modeling_Handoff_v1.0" / "data"
DEFAULT_OUT = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
SEED = 42


class UnionFind:
    def __init__(self):
        self.parent: dict[str, str] = {}

    def find(self, x: str) -> str:
        self.parent.setdefault(x, x)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return "INVALID"
    value = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    return value or "ACYCLIC"


def balanced_group_folds(groups: list[str], n_splits: int = 5) -> np.ndarray:
    """Deterministic greedy assignment balancing molecule counts, not group count."""
    members: dict[str, list[int]] = defaultdict(list)
    for idx, group in enumerate(groups):
        members[str(group)].append(idx)
    rng = np.random.RandomState(SEED)
    tie = {g: float(rng.random()) for g in members}
    ordered = sorted(members, key=lambda g: (-len(members[g]), tie[g], g))
    fold_sizes = [0] * n_splits
    assignment: dict[str, int] = {}
    for group in ordered:
        fold = min(range(n_splits), key=lambda f: (fold_sizes[f], f))
        assignment[group] = fold
        fold_sizes[fold] += len(members[group])
    return np.asarray([assignment[str(g)] for g in groups], dtype=int)


def stratified_group_folds(groups: list[str], labels: list[int],
                           requested_splits: int = 5) -> tuple[np.ndarray, int]:
    """Stratified group folds, reducing k when too few groups carry a class."""
    frame = pd.DataFrame({"group": list(map(str, groups)), "label": labels})
    support = frame.groupby(["group", "label"]).size().unstack(fill_value=0)
    class_group_counts = [(support.get(c, pd.Series(dtype=int)) > 0).sum()
                          for c in sorted(frame["label"].unique())]
    n_splits = min([requested_splits] + class_group_counts)
    if n_splits < 2:
        raise ValueError("Fewer than two independent groups carry every class")
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    fold_ids = np.full(len(frame), -1, dtype=int)
    dummy = np.zeros((len(frame), 1), dtype=np.float32)
    for fold, (_, test_idx) in enumerate(splitter.split(dummy, frame["label"], frame["group"])):
        fold_ids[test_idx] = fold
    return fold_ids, n_splits


def pubmed_years(pmids: list[str]) -> dict[str, int]:
    """Fetch publication years through NCBI ESummary in reproducible batches."""
    result: dict[str, int] = {}
    for start in range(0, len(pmids), 100):
        batch = pmids[start:start + 100]
        query = urllib.parse.urlencode({
            "db": "pubmed", "id": ",".join(batch), "retmode": "json",
            "tool": "PACER-M4-benchmark",
        })
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + query
        with urllib.request.urlopen(url, timeout=30) as response:
            payload = json.load(response)
        for pmid in batch:
            item = payload.get("result", {}).get(pmid, {})
            match = re.search(r"(?:19|20)\d{2}", str(item.get("pubdate", "")))
            if match:
                result[pmid] = int(match.group(0))
        time.sleep(0.34)
    return result


def source_table(fetch_pubmed: bool) -> tuple[pd.DataFrame, list[str]]:
    src = pd.read_csv(FROZEN / "sources.csv")
    src["year_original"] = pd.to_numeric(src.get("year"), errors="coerce")
    src["pmid_clean"] = src.get("pmid", pd.Series(index=src.index, dtype=object)).apply(
        lambda x: str(int(float(x))) if pd.notna(x) else ""
    )
    src["year_benchmark"] = src["year_original"]
    src["year_provenance"] = np.where(src["year_original"].notna(), "frozen_sources.csv", "missing")
    warnings: list[str] = []
    if fetch_pubmed:
        missing = sorted({p for p, y in zip(src["pmid_clean"], src["year_benchmark"])
                          if p and pd.isna(y)})
        try:
            fetched = pubmed_years(missing)
        except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
            fetched = {}
            warnings.append(
                "PubMed enrichment failed; temporal split remains blocked: "
                f"{type(exc).__name__}: {exc}"
            )
        for idx, row in src.iterrows():
            year = fetched.get(row["pmid_clean"])
            if pd.isna(row["year_benchmark"]) and year is not None:
                src.at[idx, "year_benchmark"] = year
                src.at[idx, "year_provenance"] = "NCBI_PubMed_ESummary"
    return src, warnings


def molecule_source_components(activity_rows: pd.DataFrame) -> tuple[dict[str, str], dict[str, int]]:
    """Connect molecules sharing sources so source-held-out folds cannot leak."""
    uf = UnionFind()
    for row in activity_rows[["canonical_molecule_id", "source_id"]].dropna().itertuples(index=False):
        uf.union("M:" + str(row.canonical_molecule_id), "S:" + str(row.source_id))
    mol_component = {
        str(m): uf.find("M:" + str(m))
        for m in activity_rows["canonical_molecule_id"].dropna().unique()
    }
    counts = Counter(mol_component.values())
    ordered = {root: idx for idx, root in enumerate(sorted(counts, key=lambda r: (-counts[r], r)))}
    return ({m: f"SCOMP{ordered[root]:04d}" for m, root in mol_component.items()}, dict(counts))


def add_common_columns(mol: pd.DataFrame, source_rows: pd.DataFrame,
                       sources: pd.DataFrame, target_col: str | None = None) -> pd.DataFrame:
    source_map = source_rows.groupby("canonical_molecule_id")["source_id"].agg(
        lambda x: "|".join(sorted(set(map(str, x.dropna()))))
    )
    primary = source_rows.groupby(["canonical_molecule_id", "source_id"]).size().reset_index(name="n")
    primary = primary.sort_values(["canonical_molecule_id", "n", "source_id"],
                                  ascending=[True, False, True]).drop_duplicates("canonical_molecule_id")
    primary_map = primary.set_index("canonical_molecule_id")["source_id"]
    year_map = sources.set_index("source_id")["year_benchmark"].to_dict()

    mol = mol.copy()
    mol["source_ids"] = mol["canonical_molecule_id"].map(source_map).fillna("")
    mol["primary_source_id"] = mol["canonical_molecule_id"].map(primary_map).fillna("")
    mol["publication_year"] = mol["primary_source_id"].map(year_map)
    mol["murcko_scaffold"] = mol["canonical_smiles"].map(scaffold)
    if target_col:
        mol["scaffold_fold"], mol["scaffold_n_splits"] = stratified_group_folds(
            mol["murcko_scaffold"].tolist(), mol[target_col].astype(int).tolist(), 5
        )
    else:
        mol["scaffold_fold"] = balanced_group_folds(mol["murcko_scaffold"].tolist())
        mol["scaffold_n_splits"] = 5

    component_map, _ = molecule_source_components(source_rows)
    mol["source_component"] = mol["canonical_molecule_id"].map(component_map).fillna("NO_SOURCE")
    if target_col:
        mol["source_fold"], mol["source_n_splits"] = stratified_group_folds(
            mol["source_component"].tolist(), mol[target_col].astype(int).tolist(), 5
        )
    else:
        mol["source_fold"] = balanced_group_folds(mol["source_component"].tolist())
        mol["source_n_splits"] = 5
    return mol


def potency_view(compounds: pd.DataFrame, sources: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(FROZEN / "modeling_potency_exact_calcium.csv")
    raw = raw.drop_duplicates("measurement_group_id")
    rows = raw.merge(compounds, on="canonical_molecule_id", how="inner")
    # The frozen modeling view already applies final_identity_valid=yes. Do not
    # re-filter on compounds.identity_status: statuses such as "confirmed" or
    # "conflated_member_excluded" can coexist with a valid canonical molecule
    # after an invalid historical member was excluded.
    rows = rows[rows["canonical_smiles"].notna()].copy()
    agg = rows.groupby("canonical_molecule_id").agg(
        canonical_smiles=("canonical_smiles", "first"),
        pEC50=("aggregated_pEC50", "mean"),
        pEC50_sd=("aggregated_pEC50", "std"),
        n_measurements=("aggregated_pEC50", "count"),
        assay_families=("assay_family", lambda x: "|".join(sorted(set(map(str, x.dropna()))))),
        orthosteric_agonists=("orthosteric_agonist", lambda x: "|".join(sorted(set(map(str, x.dropna()))))),
        cell_systems=("cell_system", lambda x: "|".join(sorted(set(map(str, x.dropna()))))),
    ).reset_index()
    agg["pEC50_sd"] = agg["pEC50_sd"].fillna(0.0)
    return add_common_columns(agg, rows, sources), rows


def classification_view(filename: str, positive_label: str, label_column: str,
                        compounds: pd.DataFrame, sources: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(FROZEN / filename)
    rows = raw.merge(compounds, on="canonical_molecule_id", how="inner")
    rows = rows[rows["canonical_smiles"].notna()].copy()
    # Conflicted molecules are excluded when the frozen view exposes the flag.
    if "molecule_label_conflict" in rows.columns:
        rows = rows[rows["molecule_label_conflict"].fillna(False) != True].copy()  # noqa: E712
    rows["target"] = (rows[label_column].astype(str) == positive_label).astype(int)
    conflicts = rows.groupby("canonical_molecule_id")["target"].nunique()
    valid = conflicts[conflicts == 1].index
    rows = rows[rows["canonical_molecule_id"].isin(valid)]
    mol = rows.groupby("canonical_molecule_id").agg(
        canonical_smiles=("canonical_smiles", "first"),
        target=("target", "first"),
        n_records=("target", "size"),
    ).reset_index()
    return add_common_columns(mol, rows, sources, "target"), rows


def fold_audit(df: pd.DataFrame, fold_col: str, target_col: str | None = None) -> dict:
    split_count_col = fold_col.removesuffix("_fold") + "_n_splits"
    report = {
        "n_splits": int(df[split_count_col].iloc[0]),
        "fold_sizes": df[fold_col].value_counts().sort_index().astype(int).to_dict(),
    }
    group_col = "murcko_scaffold" if fold_col == "scaffold_fold" else "source_component"
    leakage = df.groupby(group_col)[fold_col].nunique()
    report["groups_crossing_folds"] = int((leakage > 1).sum())
    if target_col:
        report["target_by_fold"] = pd.crosstab(df[fold_col], df[target_col]).to_dict()
    return report


def view_audit(df: pd.DataFrame, target_col: str | None = None) -> dict:
    result = {
        "n_rows": int(len(df)),
        "n_unique_molecules": int(df["canonical_molecule_id"].nunique()),
        "n_scaffolds": int(df["murcko_scaffold"].nunique()),
        "n_primary_sources": int(df["primary_source_id"].nunique()),
        "publication_year_coverage": float(df["publication_year"].notna().mean()),
        "scaffold_split": fold_audit(df, "scaffold_fold", target_col),
        "source_split": fold_audit(df, "source_fold", target_col),
    }
    if target_col:
        result["target_counts"] = {str(k): int(v) for k, v in df[target_col].value_counts().items()}
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--fetch-pubmed", action="store_true",
                        help="Fill missing source years from NCBI in benchmark metadata only")
    parser.add_argument("--temporal-cutoff", type=int, default=2020)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    compounds = pd.read_csv(FROZEN / "compounds.csv")[[
        "canonical_molecule_id", "canonical_smiles", "identity_status"
    ]]
    sources, metadata_warnings = source_table(args.fetch_pubmed)

    potency, potency_rows = potency_view(compounds, sources)
    strict, strict_rows = classification_view(
        "modeling_pam_vs_inactive.csv", "positive", "molecule_strict_label", compounds, sources
    )
    threshold, threshold_rows = classification_view(
        "modeling_threshold_10uM.csv", "threshold_positive", "class_label", compounds, sources
    )

    # A scientifically interpretable strict stress test: independently hold out
    # each source component containing at least five experimental inactives.
    # Positive-only components remain training-only in every outer split.
    negative_by_component = strict.groupby("source_component")["target"].apply(
        lambda x: int((x == 0).sum())
    )
    eligible_components = sorted(negative_by_component[negative_by_component >= 5].index)
    series_fold_map = {component: fold for fold, component in enumerate(eligible_components)}
    strict["series_holdout_fold"] = strict["source_component"].map(series_fold_map).fillna(-1).astype(int)
    strict["series_holdout_n_splits"] = len(eligible_components)

    for frame in (potency, strict, threshold):
        frame["temporal_split"] = np.where(
            frame["publication_year"].isna(), "unassigned_missing_year",
            np.where(frame["publication_year"] <= args.temporal_cutoff, "train", "test")
        )

    potency.to_csv(args.out / "potency_molecules.csv", index=False)
    strict.to_csv(args.out / "pam_vs_inactive.csv", index=False)
    threshold.to_csv(args.out / "threshold_10uM.csv", index=False)
    sources.to_csv(args.out / "source_metadata.csv", index=False)

    used_sources = set(potency_rows["source_id"].dropna()) | set(strict_rows["source_id"].dropna()) \
        | set(threshold_rows["source_id"].dropna())
    used_meta = sources[sources["source_id"].isin(used_sources)]
    missing_year = used_meta[used_meta["year_benchmark"].isna()]["source_id"].astype(str).tolist()
    temporal_ready = len(missing_year) == 0

    report = {
        "benchmark": "M4-PAM-Benchmark-v1.0-development",
        "seed": SEED,
        "temporal_cutoff": args.temporal_cutoff,
        "temporal_split_ready": temporal_ready,
        "missing_year_source_ids": sorted(missing_year),
        "potency": view_audit(potency),
        "pam_vs_inactive": view_audit(strict, "target"),
        "threshold_10uM": view_audit(threshold, "target"),
        "hard_warnings": list(metadata_warnings),
    }
    series_eval = strict[strict["series_holdout_fold"] >= 0]
    report["pam_vs_inactive"]["negative_bearing_series_holdout"] = {
        "n_splits": len(eligible_components),
        "training_only_molecules": int((strict["series_holdout_fold"] < 0).sum()),
        "test_molecules": int(len(series_eval)),
        "test_components": eligible_components,
        "target_by_fold": pd.crosstab(
            series_eval["series_holdout_fold"], series_eval["target"]
        ).to_dict(),
    }
    if not temporal_ready:
        report["hard_warnings"].append(
            "Temporal evaluation is NOT valid until every used source has a verified publication year."
        )
    if strict["target"].value_counts().min() < 100:
        report["hard_warnings"].append(
            "Strict PAM classification has fewer than 100 molecules in the minority class; use MCC, "
            "balanced accuracy, calibration and confidence intervals, not accuracy or PR-AUC alone."
        )
    strict_neg_by_source_fold = strict.groupby("source_fold")["target"].apply(lambda x: int((x == 0).sum()))
    if int(strict_neg_by_source_fold.min()) < 5:
        report["hard_warnings"].append(
            "At least one strict source-held-out fold has fewer than 5 negatives. Source-split "
            "classification is a high-variance stress test; report fold composition and bootstrap CI."
        )
    with open(args.out / "audit_report.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)

    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not temporal_ready:
        print("\nTemporal split blocked. Re-run with --fetch-pubmed, then manually verify source years.")


if __name__ == "__main__":
    main()
