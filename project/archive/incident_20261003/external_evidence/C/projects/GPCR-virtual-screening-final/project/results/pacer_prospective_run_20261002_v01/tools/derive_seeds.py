"""PACER-SEED-DERIVATION-v01 - deterministic seed set from frozen tracked data.

Rule (declared before execution, see STAGE1_RUN_SPEC_v01.json):
  per canonical_molecule_id take max(pEC50) over measurement groups,
  sort by (-pEC50, canonical_molecule_id), take the top N,
  map to canonical_smiles via compounds.csv, canonicalise with RDKit,
  write one SMILES per line.

No random component. No fallback seed list.
"""
import csv, hashlib, json, sys
from pathlib import Path
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
REPO = Path(__file__).resolve().parents[4]  # <repo>/project/results/<run>/tools/derive_seeds.py
HANDOFF = REPO / "CHRM4_PAM_Modeling_Handoff_v1.0" / "data"
OUT_SMI = REPO / "project" / "data" / "generated" / "actives.smi"
TOPN = 12


def lf_sha256(p):
    return hashlib.sha256(Path(p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    pot_path = HANDOFF / "modeling_potency_exact_calcium.csv"
    cmp_path = HANDOFF / "compounds.csv"
    with pot_path.open(newline="", encoding="utf-8-sig") as fh:
        pot_rows = list(csv.DictReader(fh))
    with cmp_path.open(newline="", encoding="utf-8-sig") as fh:
        cmp_rows = list(csv.DictReader(fh))
    smiles_by_id = {r["canonical_molecule_id"]: (r.get("canonical_smiles") or "").strip() for r in cmp_rows}

    best = {}
    for r in pot_rows:
        mid = r["canonical_molecule_id"]
        try:
            v = float(r["pEC50"])
        except (TypeError, ValueError):
            continue
        if mid not in best or v > best[mid]:
            best[mid] = v

    ranked = sorted(best.items(), key=lambda kv: (-kv[1], kv[0]))
    selected, skipped = [], []
    for mid, pec50 in ranked:
        if len(selected) >= TOPN:
            break
        smi = smiles_by_id.get(mid, "")
        if not smi or Chem.MolFromSmiles(smi) is None:
            skipped.append({"canonical_molecule_id": mid, "pEC50": pec50, "reason": "no usable canonical SMILES"})
            continue
        selected.append({"canonical_molecule_id": mid, "pEC50": pec50,
                         "smiles_input": smi, "smiles_canonical": Chem.MolToSmiles(Chem.MolFromSmiles(smi))})

    OUT_SMI.parent.mkdir(parents=True, exist_ok=True)
    OUT_SMI.write_text("\n".join(r["smiles_canonical"] for r in selected) + "\n", encoding="utf-8")
    receipt = {
        "schema": "pacer.prospective.stage1_seed_derivation.v1",
        "protocol_id": "PACER-SEED-DERIVATION-v01",
        "inputs": {
            "potency": {"path": str(pot_path.relative_to(REPO)), "sha256_lf": lf_sha256(pot_path),
                        "rows": len(pot_rows), "distinct_molecules": len(best)},
            "compounds": {"path": str(cmp_path.relative_to(REPO)), "sha256_lf": lf_sha256(cmp_path),
                          "rows": len(cmp_rows)},
        },
        "top_n_requested": TOPN,
        "n_selected": len(selected),
        "skipped": skipped,
        "fallback_seeds_used": False,
        "seed_file": {"path": str(OUT_SMI.relative_to(REPO)),
                      "sha256": hashlib.sha256(OUT_SMI.read_bytes()).hexdigest(),
                      "bytes": OUT_SMI.stat().st_size},
        "selected": selected,
    }
    out = Path(__file__).resolve().parents[1] / "STAGE1_SEED_DERIVATION_RECEIPT_v01.json"
    out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"n_selected": len(selected), "seed_file_sha256": receipt["seed_file"]["sha256"],
                      "selected_ids": [r["canonical_molecule_id"] for r in selected]}, indent=2))


if __name__ == "__main__":
    main()
