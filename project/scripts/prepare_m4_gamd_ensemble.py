# -*- coding: utf-8 -*-
"""Split and prepare the published M4R GaMD structural ensemble.

The source multi-model PDB contains the M4R receptor and the co-simulated PAM
MK-97 as residue XXX.  The heavy-atom centroid of MK-97 defines the allosteric
docking box independently in every cluster; the ligand is removed before
receptor PDBQT preparation.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np


P = Path(__file__).resolve().parents[1]
SOURCE = P / "data" / "miao2026_m4r" / "M4R_ensemble_clusters.pdb"
PMF_FILE = P / "data" / "miao2026_m4r" / "M4R_ensemble_PMF.xvg"
CLUSTER_FILE = P / "data" / "miao2026_m4r" / "M4R_ensemble_clustering_summary.dat"
OUT = P / "results" / "m4_gamd_ensemble" / "receptors"


def models(path: Path):
    current = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines(True):
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            yield current
            current = []
        elif current is not None and line.startswith(("ATOM", "HETATM", "TER")):
            current.append(line)


def xyz(lines):
    return np.asarray([[float(x[30:38]), float(x[38:46]), float(x[46:54])] for x in lines])


def heavy(line):
    element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
    return element != "H"


def pmf_values():
    out = {}
    for line in PMF_FILE.read_text().splitlines():
        if not line.strip() or line[0] in "#@":
            continue
        cluster, value = map(float, line.split()[:2])
        if cluster >= 0:
            out[int(cluster)] = value
    return out


def populations():
    out = {}
    for line in CLUSTER_FILE.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0].isdigit():
            out[int(fields[0])] = float(fields[2])
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pmf, pop = pmf_values(), populations()
    report = []
    for i, block in enumerate(models(SOURCE)):
        ligand = [x for x in block if x.startswith(("ATOM", "HETATM")) and x[17:20].strip() == "XXX"]
        receptor = [x for x in block if not (x.startswith(("ATOM", "HETATM")) and x[17:20].strip() == "XXX")]
        ligand_heavy = [x for x in ligand if heavy(x)]
        if not ligand_heavy:
            raise RuntimeError(f"cluster {i}: bound MK-97/XXX ligand is missing")
        center = xyz(ligand_heavy).mean(0)
        stem = OUT / f"cluster_{i:02d}"
        receptor_pdb = stem.with_name(stem.name + "_receptor.pdb")
        ligand_pdb = stem.with_name(stem.name + "_MK97.pdb")
        receptor_pdb.write_text("".join(receptor) + "END\n", encoding="utf-8")
        ligand_pdb.write_text("".join(ligand) + "END\n", encoding="utf-8")
        cmd = [
            "mk_prepare_receptor.exe", "--read_pdb", str(receptor_pdb),
            "-o", str(stem), "-p", "--charge_model", "gasteiger", "-a",
        ]
        run = subprocess.run(cmd, capture_output=True, text=True)
        pdbqt = stem.with_suffix(".pdbqt")
        report.append({
            "cluster": i, "population": pop.get(i), "pmf_kcal_mol": pmf.get(i),
            "box_center": center.tolist(), "box_size": [30.0, 30.0, 30.0],
            "bound_ligand": "MK-97 (source residue XXX)",
            "bound_ligand_heavy_atoms": len(ligand_heavy),
            "receptor_atoms": sum(x.startswith(("ATOM", "HETATM")) for x in receptor),
            "pdbqt_ok": run.returncode == 0 and pdbqt.exists(),
            "pdbqt": str(pdbqt), "preparation_message": (run.stderr + run.stdout)[-1000:],
        })
        print(f"cluster {i:02d}: center={center.round(3)} pdbqt_ok={report[-1]['pdbqt_ok']}", flush=True)
    payload = {
        "source": str(SOURCE), "n_clusters": len(report),
        "box_definition": "30 A cube at each representative's co-simulated MK-97 heavy-atom centroid",
        "interpretation": "Published holo GaMD ensemble; receptor flexibility approximation, not PAM efficacy labels.",
        "clusters": report,
    }
    (OUT.parent / "ensemble_preparation_audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if not all(x["pdbqt_ok"] for x in report):
        raise RuntimeError("one or more receptor clusters failed preparation; inspect audit JSON")


if __name__ == "__main__":
    main()
