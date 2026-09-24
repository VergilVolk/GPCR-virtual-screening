#!/usr/bin/env python3
"""PACER-DC reproducible trajectory QC, initial version."""

import argparse
import csv
import json
import tempfile
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import distance_array, minimize_vectors

CORE = [58, 61, 62, 77, 153, 154, 155, 231, 236, 239, 240, 243]


def statistics(values):
    a = np.asarray(values, dtype=float)
    return {
        "mean": float(a.mean()),
        "std": float(a.std()),
        "min": float(a.min()),
        "max": float(a.max()),
        "last": float(a[-1]),
    }


def windows(values):
    assert len(values) == 500
    return [statistics(values[i:i+100]) for i in range(0, 500, 100)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("context")
    parser.add_argument("--replica", type=int, required=True)
    parser.add_argument("--ach-chain", default=None)
    parser.add_argument("--candidate-chain", default=None)
    args = parser.parse_args()

    root = Path("project/results")
    top = root / "pacer_dc_membrane_reference_v01" / args.context / "minimized.pdb"
    prod = root / "pacer_dc_production_v01" / args.context / f"replica_{args.replica:02d}"
    traj = prod / "trajectory.dcd"
    state = prod / "state.csv"

    for p in [top, traj, state]:
        if not p.is_file():
            raise FileNotFoundError(p)

    # Temporary analysis topology; original PDB is unchanged.
    with tempfile.TemporaryDirectory() as tmp:
        clean = Path(tmp) / "topology.pdb"
        with top.open(errors="replace") as src, clean.open("w") as dst:
            for line in src:
                if not line.startswith("CONECT"):
                    dst.write(line)

        u = mda.Universe(str(clean), str(traj))

        ca = u.select_atoms("chainID E and name CA")
        receptor = u.select_atoms("chainID E and not name H*")

        assert len(ca) == 270
        assert len(u.trajectory) == 500

        ligands = {}
        if args.ach_chain:
            ligands["ACh"] = u.select_atoms(
                f"chainID {args.ach_chain} and not name H*"
            )
            assert len(ligands["ACh"]) == 10

        if args.candidate_chain:
            ligands["compound110"] = u.select_atoms(
                f"chainID {args.candidate_chain} and not name H*"
            )
            assert len(ligands["compound110"]) == 27

        u.trajectory[0]
        ref_ca = ca.positions.astype(float).copy()
        ref_lig = {
            name: atoms.positions.astype(float).copy()
            for name, atoms in ligands.items()
        }

        pockets = {}
        for name, atoms in ligands.items():
            d = distance_array(
                receptor.positions, atoms.positions,
                box=u.dimensions
            )
            pockets[name] = receptor[np.any(d < 5.0, axis=1)]

        core = receptor[
            np.isin(
                receptor.resindices,
                ca.resindices[np.array(CORE) - 1]
            )
        ]
        assert len(np.unique(core.resindices)) == 12

        receptor_rmsd = rms.RMSD(
            ca, ca, ref_frame=0, superposition=True
        ).run().results.rmsd[:, 2]

        results = {
            "context": args.context,
            "replica": args.replica,
            "frames": len(u.trajectory),
            "atoms": len(u.atoms),
            "receptor": {
                "ca_atoms": len(ca),
                "rmsd_A": statistics(receptor_rmsd),
                "windows": windows(receptor_rmsd),
            },
            "ligands": {},
        }

        for name, atoms in ligands.items():
            results["ligands"][name] = {
                "heavy_atoms": len(atoms),
                "initial_pocket_atoms": len(pockets[name]),
            }

        ligand_rmsd = {name: [] for name in ligands}
        contacts = {name: [] for name in ligands}
        min_dist = {name: [] for name in ligands}
        occupancy = {name: [] for name in ligands}
        raw_centroid = []
        corrected_centroid = []
        times = []
        finite = True

        for ts in u.trajectory:
            times.append(float(ts.time))
            finite &= bool(np.isfinite(ts.positions).all())

            mobile = ca.positions.astype(float)
            mc = mobile.mean(axis=0)
            rc = ref_ca.mean(axis=0)
            rot, _ = align.rotation_matrix(mobile - mc, ref_ca - rc)

            for name, atoms in ligands.items():
                xyz = atoms.positions.astype(float)
                fitted = (xyz - mc) @ rot.T + rc
                diff = fitted - ref_lig[name]
                ligand_rmsd[name].append(
                    float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))
                )

                d = distance_array(
                    pockets[name].positions,
                    atoms.positions,
                    box=ts.dimensions
                )
                contacts[name].append(
                    int(np.any(d < 4.5, axis=1).sum())
                )
                min_dist[name].append(float(d.min()))

                d_all = distance_array(
                    receptor.positions,
                    atoms.positions,
                    box=ts.dimensions
                )
                hit = set(receptor.resindices[
                    np.any(d_all < 4.5, axis=1)
                ])
                occupancy[name].append([
                    int(idx in hit) for idx in ca.resindices
                ])

            if "compound110" in ligands:
                cxyz = core.positions.astype(float)
                lxyz = ligands["compound110"].positions.astype(float)

                delta = minimize_vectors(
                    (lxyz.mean(axis=0) - cxyz.mean(axis=0))[None, :],
                    ts.dimensions
                )[0]
                raw_centroid.append(float(np.linalg.norm(delta)))

                anchor = cxyz[0]
                cu = anchor + minimize_vectors(
                    cxyz - anchor, ts.dimensions
                )
                lu = anchor + minimize_vectors(
                    lxyz - anchor, ts.dimensions
                )
                corrected_centroid.append(float(
                    np.linalg.norm(lu.mean(axis=0) - cu.mean(axis=0))
                ))

        times = np.asarray(times)
        results["trajectory_qc"] = {
            "coordinates_finite": finite,
            "first_time_ps": float(times[0]),
            "last_time_ps": float(times[-1]),
            "intervals_10ps": bool(
                np.allclose(np.diff(times), 10, atol=0.001)
            ),
        }

        for name in ligands:
            occ = np.asarray(occupancy[name])
            early = occ[:100].mean(axis=0)
            late = occ[-100:].mean(axis=0)

            records = [
                {
                    "sim_residue_index": i + 1,
                    "resname": ca[i].resname,
                    "early": float(early[i]),
                    "late": float(late[i]),
                }
                for i in range(270)
                if early[i] >= 0.2 or late[i] >= 0.2
            ]

            results["ligands"][name].update({
                "rmsd_A": statistics(ligand_rmsd[name]),
                "rmsd_windows": windows(ligand_rmsd[name]),
                "initial_contacts": contacts[name][0],
                "last_contacts": contacts[name][-1],
                "contact_windows": windows(contacts[name]),
                "last_min_distance_A": min_dist[name][-1],
                "residue_occupancy": records,
            })

        if raw_centroid:
            delta = np.abs(
                np.asarray(raw_centroid)
                - np.asarray(corrected_centroid)
            )
            results["core_centroid"] = {
                "initial_A": raw_centroid[0],
                "last_A": raw_centroid[-1],
                "windows": windows(raw_centroid),
                "pbc_max_difference_A": float(delta.max()),
                "pbc_mean_difference_A": float(delta.mean()),
            }

    with state.open(newline="") as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 500
    steps = np.array([int(r['#"Step"']) for r in rows])
    assert np.array_equal(steps, np.arange(5000, 2500001, 5000))

    results["thermodynamics"] = {}
    for key, column in {
        "potential_energy": "Potential Energy (kJ/mole)",
        "kinetic_energy": "Kinetic Energy (kJ/mole)",
        "temperature": "Temperature (K)",
        "box_volume": "Box Volume (nm^3)",
    }.items():
        results["thermodynamics"][key] = statistics(
            [float(row[column]) for row in rows]
        )

    output = root / "pacer_dc_qc_v01" / (
        f"{args.context}_replica_{args.replica:02d}_analysis_v01.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    if output.exists():
        raise FileExistsError(output)

    with output.open("x", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, allow_nan=False)

    print("OUTPUT:", output)
    print("FRAMES:", results["frames"])
    print("RECEPTOR_RMSD:", results["receptor"]["rmsd_A"]["mean"])
    for name, data in results["ligands"].items():
        print(name, "RMSD:", data["rmsd_A"]["mean"])
        print(name, "LAST_CONTACTS:", data["last_contacts"])
    if "core_centroid" in results:
        print("CORE_CENTROID:", results["core_centroid"]["last_A"])
        print("PBC_MAX_DIFFERENCE:",
              results["core_centroid"]["pbc_max_difference_A"])
    print("STATUS: ANALYSIS_COMPLETE")


if __name__ == "__main__":
    main()
