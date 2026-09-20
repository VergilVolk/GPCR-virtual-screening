"""Build the preregistration-ready Suven WO2025099660A1 external M4 PAM set.

The script reconstructs structures from patent IUPAC names through OPSIN and
checks each result with RDKit.  It deliberately does not load or run PACER.
"""

from __future__ import annotations

import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "project" / "data" / "papers" / "WO2025099660A1"
OUT_DIR = ROOT / "project" / "data" / "benchmarks" / "m4_pam_v1" / "external_suven_2025"
CACHE_DIR = SOURCE_DIR / "opsin_cache"


CRE_LUC = {
    1: 24, 3: 102, 4: 11.3, 5: 100, 6: 10, 7: 19.3, 9: 64, 11: 469,
    15: 475, 16: 238, 17: 154.3, 42: 205, 44: 16, 46: 67, 62: 32,
    65: 63, 73: 17, 79: 10, 80: 24, 86: 48.3, 87: 44, 88: 48,
    18: 136, 19: 48, 20: 30, 21: 49, 22: 355, 23: 85, 24: 12.7,
    26: 21, 28: 49, 31: 16, 33: 56, 41: 165, 91: 12, 95: 65,
    99: 17, 100: 43, 101: 112, 106: 85, 108: 32, 117: 124,
    118: 66, 119: 98, 123: 53,
}

GLOSENSOR = {
    1: 58, 4: 3.4, 7: 4, 17: 22, 19: 17, 20: 20, 21: 95,
    22: 92, 23: 31, 26: 57, 42: 36, 86: 28, 88: 19, 99: 4,
}

M2_CRE_LUC = {
    1: 3001, 7: 10000, 17: 10000, 18: 10000, 21: 3366,
    23: 10000, 26: 2282, 86: 10000,
}
M2_CENSORED = {7, 17, 18, 23, 86}


# Exact names were transcribed from the compound tables.  The common combinatorial
# series makes scaffold/headgroup swaps auditable against the source drawings.
NAMES = {
    1: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    3: "6-[4-(6-Methoxy-pyridin-3-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    4: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    5: "6-[4-(4-Methoxy-phenoxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    6: "6-[4-(2,3-Dihydro-benzofuran-5-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    7: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    9: "3-[1-(7,8-Dimethyl-[1,2,4]triazolo[4,3-b]pyridazin-6-yl)-piperidin-4-yloxy]-benzonitrile",
    11: "6-[4-(6-Fluoro-pyridin-3-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    15: "6-[4-(6-Fluoro-pyridin-3-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    16: "7-Methyl-6-[4-(6-methyl-pyridin-3-yloxy)-piperidin-1-yl]-[1,2,4]triazolo[4,3-b]pyridazine",
    17: "6-[4-(4-Methoxy-phenoxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    18: "3-[1-(7-Methyl-[1,2,4]triazolo[4,3-b]pyridazin-6-yl)-piperidin-4-yloxy]-benzonitrile",
    19: "6-[4-(2,3-Dihydro-benzofuran-5-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    20: "6-[4-(2,3-Dihydro-benzofuran-6-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    21: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    22: "6-[4-(6-Methoxy-pyridin-3-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[4,3-b]pyridazine",
    23: "6-[4-(4-Methoxy-phenoxy)-piperidin-1-yl]-3,7,8-trimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    24: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-3,7,8-trimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    26: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-3,7,8-trimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    28: "6-[4-(2,3-Dihydro-benzofuran-6-yloxy)-piperidin-1-yl]-3,7,8-trimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    31: "6-[4-(2,3-Dihydro-benzofuran-5-yloxy)-piperidin-1-yl]-3,7,8-trimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    33: "6-[4-(6-Methoxy-pyridin-3-yloxy)-piperidin-1-yl]-3,7,8-trimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    41: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-3,7-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    42: "6-[4-(6-Methoxy-pyridin-3-yloxy)-piperidin-1-yl]-3,7-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    44: "6-[4-(2,3-Dihydro-benzofuran-5-yloxy)-piperidin-1-yl]-3,7-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    46: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-3,7-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    62: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-3-methoxymethyl-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    65: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-3-methoxymethyl-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    73: "3-Difluoromethyl-6-[4-(2,3-dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[4,3-b]pyridazine",
    79: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    80: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-7,8-dimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    86: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[1,5-b]pyridazine",
    87: "6-[4-(2,3-Dihydro-benzofuran-5-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[1,5-b]pyridazine",
    88: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-7-methyl-[1,2,4]triazolo[1,5-b]pyridazine",
    91: "6-[4-(2,3-Dihydro-benzofuran-6-yloxy)-piperidin-1-yl]-2,7,8-trimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    95: "6-[4-(6-Methoxy-pyridin-3-yloxy)-piperidin-1-yl]-2,7,8-trimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    99: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-2,7,8-trimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    100: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-2,7-dimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    101: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-2,7-dimethyl-[1,2,4]triazolo[1,5-b]pyridazine",
    106: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-2,7-dimethyl-imidazo[1,2-b]pyridazine",
    108: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-2,7,8-trimethyl-imidazo[1,2-b]pyridazine",
    117: "6-[4-(1,3-Dihydro-isobenzofuran-5-yloxy)-piperidin-1-yl]-7-methyl-imidazo[1,2-b]pyridazine",
    118: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7-methyl-imidazo[1,2-b]pyridazine",
    119: "6-[4-(2,3-Dihydro-benzofuran-5-yloxy)-piperidin-1-yl]-7-methyl-imidazo[1,2-b]pyridazine",
    123: "6-[4-(2,3-Dihydro-benzo[1,4]dioxin-6-yloxy)-piperidin-1-yl]-7,8-dimethyl-imidazo[1,2-b]pyridazine",
}


def opsin(name: str, compound_id: int) -> dict[str, str]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / f"example_{compound_id:03d}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    url = "https://opsin.ch.cam.ac.uk/opsin/" + urllib.parse.quote(name, safe="") + ".json"
    request = urllib.request.Request(url, headers={"User-Agent": "PACER-M4 reproducibility audit"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        payload = {"status": "FAILURE", "message": str(error), "smiles": ""}
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    time.sleep(0.15)
    return payload


def main() -> None:
    if set(CRE_LUC) != set(NAMES):
        raise RuntimeError(f"Name/label mismatch: labels={len(CRE_LUC)} names={len(NAMES)}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = []
    for compound_id in sorted(CRE_LUC):
        payload = opsin(NAMES[compound_id], compound_id)
        smiles = payload.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True) if mol else ""
        mw = Descriptors.MolWt(mol) if mol else math.nan
        records.append(
            {
                "compound_id": compound_id,
                "name": NAMES[compound_id],
                "smiles": canonical,
                "cre_luc_ec50_nM": CRE_LUC[compound_id],
                "cre_luc_pEC50": 9 - math.log10(CRE_LUC[compound_id]),
                "glosensor_ec50_nM": GLOSENSOR.get(compound_id, math.nan),
                "glosensor_pEC50": 9 - math.log10(GLOSENSOR[compound_id]) if compound_id in GLOSENSOR else math.nan,
                "m2_cre_luc_ec50_nM": M2_CRE_LUC.get(compound_id, math.nan),
                "m2_censoring": ">" if compound_id in M2_CENSORED else ("exact" if compound_id in M2_CRE_LUC else ""),
                "m4_m2_selectivity_fold_lower_bound": (
                    M2_CRE_LUC[compound_id] / CRE_LUC[compound_id] if compound_id in M2_CRE_LUC else math.nan
                ),
                "molecular_formula": rdMolDescriptors.CalcMolFormula(mol) if mol else "",
                "mw": mw,
                "opsin_status": payload.get("status", ""),
                "structure_qc": "opsin_rdkit_pass" if mol else "failed",
                "source": "WO2025099660A1 Table 1A/1B/1C",
                "assay_probe": "acetylcholine EC20",
            }
        )
    frame = pd.DataFrame(records)
    frame.to_csv(OUT_DIR / "external_patent_potency.csv", index=False)
    audit = {
        "patent": "WO2025099660A1",
        "dataset_role": "independent_external_holdout",
        "n_cre_luc": int(frame.cre_luc_ec50_nM.notna().sum()),
        "n_glosensor": int(frame.glosensor_ec50_nM.notna().sum()),
        "n_m2": int(frame.m2_cre_luc_ec50_nM.notna().sum()),
        "n_opsin_rdkit_pass": int((frame.structure_qc == "opsin_rdkit_pass").sum()),
        "n_failed": int((frame.structure_qc == "failed").sum()),
        "labels_used_for_model_tuning": False,
        "predictions_run": False,
        "primary_endpoint": "human M4 CRE-Luc PAM EC50 at ACh EC20",
        "secondary_endpoints": ["human M4 GloSensor PAM EC50 at ACh EC20", "human M2 CRE-Luc EC50"],
        "claim_boundary": "Patent functional data; structures require image/name/mass QC before unblinding models.",
    }
    (OUT_DIR / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
