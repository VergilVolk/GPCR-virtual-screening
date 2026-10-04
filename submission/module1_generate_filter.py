from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


ROOT = Path(__file__).resolve().parent


def build_filter_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


def evaluate_library(input_path: Path, output_path: Path) -> pd.DataFrame:
    table = pd.read_csv(input_path)
    catalog = build_filter_catalog()
    records = []
    for row in table.itertuples(index=False):
        selection_rank = getattr(row, "selection_rank", None)
        molecule = Chem.MolFromSmiles(row.smiles)
        if molecule is None:
            records.append({"selection_rank": selection_rank, "candidate_id": row.candidate_id, "smiles": row.smiles,
                            "valid": False, "chemistry_pass": False})
            continue
        canonical = Chem.MolToSmiles(molecule)
        mw = Descriptors.MolWt(molecule)
        logp = Crippen.MolLogP(molecule)
        hbd = Lipinski.NumHDonors(molecule)
        hba = Lipinski.NumHAcceptors(molecule)
        tpsa = Descriptors.TPSA(molecule)
        pains = catalog.GetFirstMatch(molecule) is not None
        chemistry_pass = (
            180 <= mw <= 650 and -1.0 <= logp <= 6.5 and hbd <= 5
            and hba <= 12 and tpsa <= 160 and not pains
        )
        records.append({
            "selection_rank": selection_rank, "candidate_id": row.candidate_id,
            "smiles": canonical, "valid": True,
            "mw": round(mw, 3), "logp": round(logp, 3), "hbd": hbd, "hba": hba,
            "tpsa": round(tpsa, 3), "qed": round(QED.qed(molecule), 4),
            "pains_alert": pains, "chemistry_pass": chemistry_pass,
        })
    result = pd.DataFrame(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Medicinal chemistry checks for a candidate library")
    parser.add_argument("--input", type=Path, default=ROOT / "data/demo/candidate_library.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/module1_chemistry.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = evaluate_library(args.input, args.output)
    print(f"Chemistry checks: {int(result.chemistry_pass.sum())}/{len(result)} passed")
    print(args.output)


if __name__ == "__main__":
    main()
