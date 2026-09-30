from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import dock_candidate_portfolio as docking  # noqa: E402

CID = "CM00734"
SMILES = "COc1ccc(C)c2sc(NC(=O)NC3CCC3)nc12"

SEED = 42
EXHAUSTIVENESS = 4

OUT = ROOT / "results" / "pacer_dc_reference_complexes_v01"
WORK = OUT / "cm00734_stage_b_pose_v01"

FINAL_SDF = OUT / "CM00734_7TRS_redocked.sdf"
AUDIT = OUT / "CM00734_STAGE_B_POSE_AUDIT_v01.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    required = [
        docking.RECEPTOR,
        docking.RECEPTOR_PDB,
        docking.VINA,
        SCRIPTS / "dock_candidate_portfolio.py",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing required frozen docking asset(s):\n" + "\n".join(missing))

    # Reuse the existing frozen pose-proposal implementation exactly.
    docking.POSES = WORK

    result = docking.one((CID, SMILES, SEED, EXHAUSTIVENESS))
    if result.get("error"):
        raise RuntimeError(f"CM00734 docking failed: {result['error']}")

    key = hashlib.sha1(f"{CID}|{SEED}".encode()).hexdigest()[:14]
    pose_pdbqt = WORK / f"{key}.pose.pdbqt"
    ligand_pdbqt = WORK / f"{key}.lig.pdbqt"

    if not pose_pdbqt.exists():
        raise RuntimeError(f"Expected Vina pose missing: {pose_pdbqt}")

    # Meeko reconstructs chemistry using the SMILES/index remarks retained
    # through Vina, avoiding bond-order guessing from PDBQT coordinates.
    pdbqt_mol = PDBQTMolecule.from_file(str(pose_pdbqt), skip_typing=True)
    molecules = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    molecules = [m for m in molecules if m is not None]

    if len(molecules) != 1:
        raise RuntimeError(
            f"Expected exactly one reconstructed ligand; got {len(molecules)}"
        )

    pose = molecules[0]
    pose.SetProp("_Name", CID)

    writer = Chem.SDWriter(str(FINAL_SDF))
    writer.write(pose)
    writer.close()

    check = Chem.SDMolSupplier(str(FINAL_SDF), removeHs=False)[0]
    if check is None or check.GetNumConformers() != 1:
        raise RuntimeError("Written SDF failed RDKit replay")

    audit = {
        "schema": "pacer_dc.cm00734.stage_b_pose.v1",
        "status": "CM00734_STAGE_B_POSE_PROPOSAL_COMPLETE",
        "candidate_id": CID,
        "canonical_smiles": SMILES,
        "role": "strict_experimental_inactive_binding_to_be_tested",
        "interpretation": (
            "Frozen 7TRS allosteric-pocket pose proposal only; "
            "Vina affinity and geometry are not PAM efficacy evidence."
        ),
        "protocol": {
            "implementation": str(SCRIPTS / "dock_candidate_portfolio.py"),
            "implementation_sha256": sha256(
                SCRIPTS / "dock_candidate_portfolio.py"
            ),
            "receptor_pdbqt": str(docking.RECEPTOR),
            "receptor_pdbqt_sha256": sha256(docking.RECEPTOR),
            "receptor_pdb": str(docking.RECEPTOR_PDB),
            "receptor_pdb_sha256": sha256(docking.RECEPTOR_PDB),
            "vina": str(docking.VINA),
            "vina_sha256": sha256(docking.VINA),
            "center_A": list(docking.CENTER),
            "size_A": list(docking.SIZE),
            "seed": SEED,
            "exhaustiveness": EXHAUSTIVENESS,
            "num_modes": 1,
            "pocket_residues": list(docking.POCKET),
            "contact_cutoff_A": docking.CONTACT,
        },
        "pose_result": result,
        "artifacts": {
            "input_ligand_pdbqt": str(ligand_pdbqt),
            "input_ligand_pdbqt_sha256": sha256(ligand_pdbqt),
            "vina_pose_pdbqt": str(pose_pdbqt),
            "vina_pose_pdbqt_sha256": sha256(pose_pdbqt),
            "pose_sdf": str(FINAL_SDF),
            "pose_sdf_sha256": sha256(FINAL_SDF),
        },
        "selection_policy": {
            "pose_selection": "single Vina mode from pre-existing frozen protocol",
            "outcome_driven_pose_selection": False,
            "LY2119620_result_used_for_pose_selection": False,
            "new_geometry_threshold": False,
        },
    }

    AUDIT.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(json.dumps({
        "status": audit["status"],
        "vina_affinity": result.get("vina_affinity"),
        "pocket_residue_coverage": result.get("pocket_residue_coverage"),
        "contacted_residues": result.get("contacted_residues"),
        "pose_sdf": str(FINAL_SDF),
        "pose_sdf_sha256": audit["artifacts"]["pose_sdf_sha256"],
        "audit": str(AUDIT),
    }, indent=2))


if __name__ == "__main__":
    main()