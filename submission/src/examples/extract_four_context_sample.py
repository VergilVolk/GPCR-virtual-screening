"""Create or verify a small, traceable excerpt of the actual Stage4 MD inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import MDAnalysis as mda
import numpy as np

SYSTEMS = (
    "cluster0__apo", "cluster0__probe_only",
    "PACER0073__candidate_no_probe", "PACER0073__candidate_probe",
)
FRAME_IDS = list(range(0, 100, 5))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(output: Path) -> dict:
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    if len(manifest["systems"]) != 4:
        raise ValueError("Expected four contexts")
    for record in manifest["systems"]:
        for artifact in record["sample_files"].values():
            path = output / artifact["path"]
            if path.stat().st_size != artifact["bytes"] or sha256(path) != artifact["sha256"]:
                raise ValueError(f"Sample identity mismatch: {path}")
        top = output / record["sample_files"]["topology"]["path"]
        traj = output / record["sample_files"]["trajectory"]["path"]
        universe = mda.Universe(str(top), str(traj))
        try:
            if len(universe.trajectory) != 20 or universe.atoms.n_atoms != record["atoms"]:
                raise ValueError("Sample shape mismatch")
            observed_times = [float(ts.time) for ts in universe.trajectory]
            if not np.allclose(observed_times, record["source_times_ps"], rtol=0, atol=1e-3):
                raise ValueError("Sample timestamps mismatch")
            for ts in universe.trajectory:
                if not np.isfinite(ts.positions).all():
                    raise ValueError("Nonfinite sample coordinates")
        finally:
            universe.trajectory.close()
    return {"status": "PASS", "contexts": 4, "frames_per_context": 20}


def extract(source: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    output.mkdir(parents=True)
    records = []
    for system in SYSTEMS:
        top = source / "systems" / system / "minimized.pdb"
        traj = source / "production" / system / "replica_01" / "trajectory.dcd"
        before = {"topology": sha256(top), "trajectory": sha256(traj)}
        folder = output / system
        folder.mkdir()
        sample_top = folder / "minimized.pdb"
        sample_traj = folder / "trajectory_sample.dcd"
        shutil.copyfile(top, sample_top)
        universe = mda.Universe(str(top), str(traj))
        coordinates, cells, times = [], [], []
        try:
            if len(universe.trajectory) != 1000 or not np.isclose(universe.trajectory.dt, 10, atol=1e-3):
                raise ValueError("Unexpected Stage4 temporal layout")
            for frame in FRAME_IDS:
                ts = universe.trajectory[frame]
                coordinates.append(ts.positions.copy())
                cells.append(ts.dimensions.copy())
                times.append(float(ts.time))
            # nsavc=5 and istart=1 preserve the source's 10 ps first time
            # and the frozen analysis stride of 50 ps without interpolation.
            with mda.Writer(str(sample_traj), n_atoms=universe.atoms.n_atoms,
                            dt=universe.trajectory.dt * 5, nsavc=5, istart=1) as writer:
                for frame in FRAME_IDS:
                    universe.trajectory[frame]
                    writer.write(universe.atoms)
            atoms = universe.atoms.n_atoms
        finally:
            universe.trajectory.close()
        sample = mda.Universe(str(sample_top), str(sample_traj))
        try:
            for index, ts in enumerate(sample.trajectory):
                if not np.array_equal(ts.positions, coordinates[index]):
                    raise ValueError("Sample differs from original coordinates")
                if not np.allclose(ts.dimensions, cells[index], atol=1e-4, rtol=0):
                    raise ValueError("Sample unit cell differs from original")
        finally:
            sample.trajectory.close()
        if before != {"topology": sha256(top), "trajectory": sha256(traj)}:
            raise ValueError("Source changed during extraction")
        files = {key: {"path": path.relative_to(output).as_posix(),
                       "bytes": path.stat().st_size, "sha256": sha256(path)}
                 for key, path in (("topology", sample_top), ("trajectory", sample_traj))}
        records.append({"system": system, "replica": 1, "seed": 27101, "atoms": atoms,
                        "source_topology": top.relative_to(source).as_posix(),
                        "source_trajectory": traj.relative_to(source).as_posix(),
                        "source_sha256": before, "raw_frame_ids": FRAME_IDS,
                        "source_times_ps": times, "sample_files": files,
                        "coordinates_bitwise_equal": True, "source_unchanged": True})
        print(f"Extracted and verified: {system}", flush=True)
    manifest = {"candidate": "PACER0073", "receptor_cluster": 0,
                "source_root_at_extraction": str(source.resolve()),
                "sample_kind": "first frozen-stride block of production R1",
                "claim_boundary": "Input demonstration only; not full Stage4 authentication, three-replica evaluation or PAM evidence.",
                "coordinate_unit": "angstrom", "time_unit": "ps", "interpolation": False,
                "MDAnalysis_version": mda.__version__, "systems": records}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return verify(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if not args.verify_only and args.source_root is None:
        parser.error("--source-root is required for extraction")
    result = verify(args.output_dir) if args.verify_only else extract(args.source_root, args.output_dir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
