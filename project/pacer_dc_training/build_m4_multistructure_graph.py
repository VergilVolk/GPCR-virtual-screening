#!/usr/bin/env python
"""Build a frozen M4 residue graph and ligand-contact regions from PDB structures.

The graph is deliberately label-free.  It uses only experimental structures,
geometric contacts and a small, literature-declared residue panel.  Functional
labels and PACER candidate scores must not be used to alter this artifact.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


PAM_STRUCTURES = {"7TRP": "IUE", "7TRQ": "IUI", "7V68": "2CU"}
AGONIST_CONTROL = {"7V6A": "5XI"}
ORTHOSTERIC_STRUCTURES = {"7TRK": "IXO", "7TRS": "ACH"}

# Declared before candidate scoring.  Residues are M4 sequence numbers.
LITERATURE_REGIONS = {
    "cooperativity_mutagenesis": [89, 92, 93, 95, 96, 184, 186, 423, 432, 433, 435, 439, 443],
    "orthosteric_activation_core": [112, 113, 116, 117, 164, 203, 204, 413, 416, 417, 439, 442, 443],
    "intracellular_microswitches": [129, 130, 131, 409, 413, 415, 445, 450, 453],
    # Predeclared spatially distal control used by the previous G2 audit.
    "distal_control": [221, 222, 223, 224, 225, 226, 392, 393, 394, 395, 396, 397],
}


def parse_mapping(path: Path):
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    by_resid = {int(r["structure_resid"]): int(r["embedding_index"]) for r in rows}
    site = {int(r["structure_resid"]): f'{r["aa"]}{r["structure_resid"]}' for r in rows}
    return rows, by_resid, site


def parse_pdb(path: Path, ligand: str, chain: str = "R"):
    residues = defaultdict(list)
    ca = {}
    ligand_xyz = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        record = line[:6].strip()
        if line[21:22] != chain:
            continue
        resname = line[17:20].strip()
        try:
            resid = int(line[22:26])
            xyz = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        if record == "ATOM":
            residues[resid].append(xyz)
            if line[12:16].strip() == "CA":
                ca[resid] = xyz
        elif record == "HETATM" and resname == ligand:
            element = line[76:78].strip().upper()
            if element != "H":
                ligand_xyz.append(xyz)
    if not ligand_xyz:
        raise ValueError(f"{path.name}: ligand {ligand} not found in chain {chain}")
    return residues, ca, np.asarray(ligand_xyz)


def ligand_contacts(residues, ligand_xyz, cutoff):
    hits = []
    for resid, atoms in residues.items():
        d = np.linalg.norm(np.asarray(atoms)[:, None, :] - ligand_xyz[None, :, :], axis=-1).min()
        if d <= cutoff:
            hits.append((resid, float(d)))
    return dict(hits)


def entries(resids, by_resid, site):
    return [
        {"site": site[r], "structure_resid": r, "embedding_index": by_resid[r]}
        for r in sorted(set(resids)) if r in by_resid
    ]


def kabsch_align(mobile, reference):
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    return (mobile - mobile_center) @ rotation + reference_center


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb-dir", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--contact-cutoff", type=float, default=4.5)
    p.add_argument("--graph-cutoff", type=float, default=8.0)
    args = p.parse_args()

    rows, by_resid, site = parse_mapping(args.mapping)
    structure_specs = PAM_STRUCTURES | AGONIST_CONTROL | ORTHOSTERIC_STRUCTURES
    contacts = {}
    ca_by_structure = {}
    ligand_xyz_by_structure = {}
    for pdb_id, ligand in structure_specs.items():
        residues, ca, ligand_xyz = parse_pdb(args.pdb_dir / f"{pdb_id}.pdb", ligand)
        contacts[pdb_id] = ligand_contacts(residues, ligand_xyz, args.contact_cutoff)
        ca_by_structure[pdb_id] = ca
        ligand_xyz_by_structure[pdb_id] = ligand_xyz

    pam_counts = defaultdict(int)
    for pdb_id in PAM_STRUCTURES:
        for resid in contacts[pdb_id]:
            pam_counts[resid] += 1
    pam_union = set(pam_counts)
    pam_consensus = {r for r, n in pam_counts.items() if n >= 2}
    agonist_only = set(contacts["7V6A"]) - pam_consensus
    ortho_union = set().union(*(set(contacts[p]) for p in ORTHOSTERIC_STRUCTURES))

    regions = {
        "pam_contact_consensus": entries(pam_consensus, by_resid, site),
        "pam_contact_union": entries(pam_union, by_resid, site),
        "compound110_extension": entries(agonist_only, by_resid, site),
        "orthosteric_contact_union": entries(ortho_union, by_resid, site),
    }
    for name, resids in LITERATURE_REGIONS.items():
        regions[name] = entries(resids, by_resid, site)

    # Consensus residue graph: backbone edges plus non-local CA contacts observed
    # in at least two structures.  This avoids making any one ligand structure
    # the sole definition of the allosteric pathway.
    edge_counts = defaultdict(int)
    for ca in ca_by_structure.values():
        present = sorted(r for r in ca if r in by_resid)
        for i, a in enumerate(present):
            for b in present[i + 1:]:
                if b == a + 1:
                    edge_counts[(a, b)] += 1
                elif np.linalg.norm(ca[a] - ca[b]) <= args.graph_cutoff:
                    edge_counts[(a, b)] += 1
    edges = []
    for (a, b), count in sorted(edge_counts.items()):
        is_backbone = b == a + 1
        if is_backbone or count >= 2:
            edges.append({
                "source": by_resid[a], "target": by_resid[b],
                "source_site": site[a], "target_site": site[b],
                "support_n": count, "support_fraction": count / len(structure_specs),
                "edge_type": "backbone" if is_backbone else "spatial_consensus",
            })

    # A label-free, structurally stable control panel.  Structures are aligned
    # on common C-alpha atoms; residues near any ligand or declared mechanism
    # panel are excluded.  This is safer than using the flexible ICL3 truncation
    # boundary as the sole negative control.
    common = sorted(set.intersection(*(set(ca) for ca in ca_by_structure.values())) & set(by_resid))
    reference = np.stack([ca_by_structure["7TRS"][r] for r in common])
    aligned = {}
    for pdb_id, ca in ca_by_structure.items():
        mobile = np.stack([ca[r] for r in common])
        transformed = kabsch_align(mobile, reference)
        aligned[pdb_id] = {r: transformed[i] for i, r in enumerate(common)}
    declared = set().union(*LITERATURE_REGIONS.values(), pam_union, agonist_only, ortho_union)
    candidates = []
    for resid in common:
        if resid in declared:
            continue
        # Stay away from every crystallographic ligand.
        if any(np.linalg.norm(ca_by_structure[p][resid][None, :] - ligand_xyz_by_structure[p], axis=1).min() < 12.0
               for p in structure_specs):
            continue
        xyz = np.stack([aligned[p][resid] for p in structure_specs])
        variance = float(np.mean(np.sum((xyz - xyz.mean(axis=0)) ** 2, axis=1)))
        candidates.append((variance, resid))
    stable_control = [resid for _, resid in sorted(candidates)[:12]]
    regions["stable_core_control"] = entries(stable_control, by_resid, site)

    result = {
        "version": "M4_MULTISTRUCTURE_GRAPH_v01",
        "evidence_level": "label_free_structural_prior",
        "structures": {
            "pam": PAM_STRUCTURES,
            "allosteric_agonist_control": AGONIST_CONTROL,
            "orthosteric": ORTHOSTERIC_STRUCTURES,
        },
        "contact_cutoff_A": args.contact_cutoff,
        "graph_cutoff_A": args.graph_cutoff,
        "contact_residues": {
            pdb_id: [site[r] for r in sorted(hit) if r in site] for pdb_id, hit in contacts.items()
        },
        "regions": regions,
        "stable_control_definition": "12 lowest cross-structure aligned CA-variance residues, >=12 A from every ligand and outside declared mechanism panels",
        "nodes": [{"embedding_index": int(r["embedding_index"]), "site": f'{r["aa"]}{r["structure_resid"]}'} for r in rows],
        "edges": edges,
        "claim_boundary": "Structural prior only; contact membership does not establish PAM function or causal transmission.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "pam_consensus": [site[r] for r in sorted(pam_consensus)],
        "pam_union_n": len(pam_union),
        "compound110_extension": [site[r] for r in sorted(agonist_only) if r in site],
        "edges": len(edges),
    }, indent=2))


if __name__ == "__main__":
    main()
