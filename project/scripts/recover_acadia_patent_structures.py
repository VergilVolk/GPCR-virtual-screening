"""Recover WO2025122811 examples from patent names and validate by LC-MS.

Names are transcribed from Table 1 before model evaluation. OPSIN supplies a
machine-readable structure. A structure is high-confidence only when RDKit can
parse it and calculated [M+H]+ agrees with the patent value within 0.6 Da.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


P = Path(__file__).resolve().parents[1]
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_acadia_2025"
PROTON = 1.007276466621

CORE_R = "(R)-(4-(azetidin-1-yl)-2-methyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)"
CORE_S = "(S)-(4-(azetidin-1-yl)-2-methyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)"
CORE_RR = "((R)-4-(azetidin-1-yl)-2,5-dimethyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)"


def ketone(core: str, substituent: str) -> str:
    return f"{core}(({substituent})pyrrolidin-3-yl)methanone"


def ketone_rr(aryl: str) -> str:
    return f"{CORE_RR}((R)-1-({aryl})pyrrolidin-3-yl)methanone"


ENTRIES = [
    (1, ketone(CORE_S, "1-(pyrimidin-4-yl)"), 366.2),
    (2, ketone(CORE_R, "1-(pyrimidin-4-yl)"), 366.2),
    (3, ketone(CORE_R, "1-(5-chloro-2-(trifluoromethyl)pyridin-4-yl)"), 467.1),
    (4, ketone(CORE_R, "1-(5-chloro-2-methoxypyridin-4-yl)"), 429.1),
    (5, ketone(CORE_R, "1-(2-(difluoromethoxy)pyridin-4-yl)"), 431.3),
    (6, ketone(CORE_R, "1-(2-(trifluoromethyl)pyridin-4-yl)"), 433.1),
    (7, ketone(CORE_R, "1-(5-(difluoromethoxy)pyridin-3-yl)"), 431.1),
    (8, ketone(CORE_R, "1-(2-(difluoromethyl)pyridin-4-yl)"), 415.3),
    (9, ketone(CORE_R, "1-(5-(trifluoromethyl)pyridin-3-yl)"), 433.1),
    (10, ketone(CORE_R, "1-(6-(trifluoromethyl)pyridin-3-yl)"), 433.1),
    (11, ketone(CORE_R, "1-(2-methoxypyridin-4-yl)"), 395.3),
    (12, ketone(CORE_R, "1-(2-chloropyridin-4-yl)"), 399.3),
    (13, ketone(CORE_R, "1-(2-chloropyridin-4-yl)"), 399.1),
    (14, ketone(CORE_R, "1-(pyrazolo[1,5-a]pyridin-3-yl)"), 404.2),
    (15, "2-(3-((4-(azetidin-1-yl)-2-methyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)carbonyl)pyrrolidin-1-yl)-6-(trifluoromethyl)benzonitrile", 457.2),
    (16, ketone(CORE_R, "1-(2,5-dichloropyridin-4-yl)"), 433.1),
    (17, ketone(CORE_R, "1-(5-chloro-3-fluoropyridin-2-yl)"), 417.1),
    (18, "6-(3-((4-(azetidin-1-yl)-2-methyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)carbonyl)pyrrolidin-1-yl)-5-chloronicotinonitrile", 424.1),
    (19, ketone(CORE_R, "1-(2-(trifluoromethyl)pyrimidin-4-yl)"), 434.1),
    (20, ketone(CORE_R, "1-(6-(trifluoromethyl)pyrimidin-4-yl)"), 434.1),
    (21, "6-((R)-3-((4-(azetidin-1-yl)-2-methyl-6,7-dihydro-5H-pyrrolo[3,4-d]pyrimidine-6-carbonyl)pyrrolidin-1-yl)pyrimidine-2-carbonitrile", 391.1),
    (22, ketone(CORE_R, "1-(6-(trifluoromethyl)pyridin-2-yl)"), 433.1),
    (23, ketone(CORE_R, "1-(4-(trifluoromethyl)pyridin-2-yl)"), 433.1),
    (24, ketone(CORE_R, "1-(5-methoxy-2-(trifluoromethyl)pyrimidin-4-yl)"), 464.2),
    (25, ketone(CORE_R, "1-(4-(difluoromethyl)-3-fluoropyridin-2-yl)"), 433.1),
    (26, "6-(3-((4-(azetidin-1-yl)-2-methyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)carbonyl)pyrrolidin-1-yl)nicotinonitrile", 390.1),
    (27, ketone(CORE_R, "1-(2-methyl-6-(trifluoromethyl)pyridin-4-yl)"), 447.4),
    (28, ketone_rr("pyrimidin-4-yl"), 380.3),
    (29, ketone(CORE_R, "1-(2-(difluoromethyl)-6-methylpyridin-4-yl)"), 429.3),
    (30, ketone_rr("2-chloropyridin-4-yl"), 413.2),
    (31, ketone_rr("2-(trifluoromethyl)pyridin-4-yl"), 447.3),
    (32, ketone_rr("2-(difluoromethoxy)pyridin-4-yl"), 445.3),
    (33, ketone_rr("2-(difluoromethoxy)pyridin-4-yl"), 445.2),
    (34, "4-((R)-3-(((R)-4-(azetidin-1-yl)-2,5-dimethyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)carbonyl)pyrrolidin-1-yl)pyrimidine-2-carbonitrile", 405.3),
    (35, ketone_rr("6-(trifluoromethyl)pyrimidin-4-yl"), 448.3),
    (36, "4-((R)-3-((4-(azetidin-1-yl)-2-methyl-6,7-dihydro-5H-pyrrolo[3,4-d]pyrimidine-6-carbonyl)pyrrolidin-1-yl)pyrimidine-4-carbonitrile", 391.3),
    (37, ketone(CORE_R, "1-(5-chloropyridin-3-yl)"), 399.3),
    (38, ketone_rr("5-chloropyridin-3-yl"), 413.3),
    (39, ketone_rr("2-methoxypyridin-4-yl"), 409.2),
    (40, ketone_rr("2-(difluoromethyl)pyridin-4-yl"), 429.3),
    (41, "(R)-(4-(azetidin-1-yl)-2-ethyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)(1-(2-chloropyridin-4-yl)pyrrolidin-3-yl)methanone", 413.1),
    (42, "(R)-(4-(azetidin-1-yl)-2-ethyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)(1-(5-chloropyridin-3-yl)pyrrolidin-3-yl)methanone", 413.1),
    (43, ketone_rr("2,5-dichloropyridin-4-yl"), 447.1),
    (44, "6-((R)-3-(((R)-4-(azetidin-1-yl)-2,5-dimethyl-5,7-dihydro-6H-pyrrolo[3,4-d]pyrimidin-6-yl)carbonyl)pyrrolidin-1-yl)-5-chloronicotinonitrile", 438.2),
    (45, ketone_rr("5-chloro-3-fluoropyridin-2-yl"), 431.2),
    (46, ketone_rr("5-(trifluoromethyl)pyridin-3-yl"), 447.2),
    (47, ketone_rr("6-(trifluoromethyl)pyridin-3-yl"), 447.2),
    (48, ketone_rr("5-chloro-6-methylpyridin-3-yl"), 427.2),
    (49, ketone_rr("2-methyl-6-(trifluoromethyl)pyridin-4-yl"), 461.2),
]


def opsin(name: str) -> dict:
    url = "https://opsin.ch.cam.ac.uk/opsin/" + urllib.parse.quote(name, safe="") + ".json"
    request = urllib.request.Request(url, headers={"User-Agent": "PACER-M4/2.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> None:
    cache_path = OUT / "structures_opsin_audit.csv"
    cache = {}
    if cache_path.exists():
        previous = pd.read_csv(cache_path).fillna("")
        cache = {str(row.patent_name): row._asdict() for row in previous.itertuples(index=False)
                 if row.structure_status == "high_confidence"}
    rows = []
    for example_id, name, reported_mh in ENTRIES:
        row = {"example_id": example_id, "patent_name": name, "reported_mh": reported_mh}
        if name in cache and abs(float(cache[name]["reported_mh"]) - reported_mh) < 1e-6:
            preserved = cache[name]
            preserved.update(row)
            rows.append(preserved)
            print(example_id, "high_confidence", preserved.get("mass_delta_da"), "cached", flush=True)
            continue
        try:
            result = opsin(name)
            smiles = result.get("smiles")
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            if molecule is None:
                raise ValueError("OPSIN returned no RDKit-parseable SMILES")
            canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
            calculated_mh = float(Descriptors.ExactMolWt(molecule) + PROTON)
            delta = calculated_mh - reported_mh
            row.update({
                "canonical_smiles": canonical,
                "molecular_formula": rdMolDescriptors.CalcMolFormula(molecule),
                "calculated_mh": calculated_mh,
                "mass_delta_da": delta,
                "structure_status": "high_confidence" if abs(delta) <= 0.6 else "mass_mismatch",
                "error": "",
            })
        except Exception as exc:  # preserve failures for manual repair
            row.update({"canonical_smiles": "", "molecular_formula": "", "calculated_mh": None,
                        "mass_delta_da": None, "structure_status": "opsin_failed", "error": str(exc)})
        rows.append(row)
        print(example_id, row["structure_status"], row.get("mass_delta_da"), row.get("error", ""), flush=True)
        time.sleep(0.05)

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "structures_opsin_audit.csv", index=False)
    summary = {
        "n_examples": len(frame),
        "status_counts": frame.structure_status.value_counts().to_dict(),
        "entry_policy": "Only high_confidence rows may enter the primary external benchmark.",
        "claim_boundary": "Name-to-structure recovery validated by patent LC-MS; stereochemical image spot-check remains required.",
    }
    (OUT / "structure_recovery_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
