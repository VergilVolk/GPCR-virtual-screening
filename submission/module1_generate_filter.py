from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

ROOT = Path(__file__).resolve().parent

def evaluate_library(input_path: Path, output_path: Path) -> pd.DataFrame:
    table = pd.read_csv(input_path).rename(columns={"molecule_id": "candidate_id"})
    if "candidate_id" not in table.columns:
        table.insert(0, "candidate_id", [f"PACERGEN{i:05d}" for i in range(1, len(table) + 1)])
    rows = []
    for row in table.itertuples(index=False):
        mol = Chem.MolFromSmiles(row.smiles)
        if mol is None:
            rows.append({"candidate_id": row.candidate_id, "smiles": row.smiles,
                         "valid": False, "chemistry_pass": False})
            continue
        mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
        hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
        passed = (250 <= mw <= 550 and logp <= 5.5 and hbd <= 5 and hba <= 10
                  and rotb <= 10 and 40 <= tpsa <= 140 and aromatic >= 2
                  and hba >= 2 and hbd >= 1)
        rows.append({"candidate_id": row.candidate_id, "smiles": Chem.MolToSmiles(mol),
                     "valid": True, "chemistry_pass": passed, "mw": round(mw, 3),
                     "logp": round(logp, 3), "hbd": hbd, "hba": hba,
                     "rotatable_bonds": rotb, "tpsa": round(tpsa, 3),
                     "aromatic_rings": aromatic, "qed": round(QED.qed(mol), 4)})
    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    return result

def main() -> None:
    parser = argparse.ArgumentParser(description="Recheck the frozen generated library")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "results/module1_chemistry.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = evaluate_library(args.input, args.output)
    print(f"Chemistry checks: {int(result.chemistry_pass.sum())}/{len(result)} passed")

if __name__ == "__main__":
    main()
