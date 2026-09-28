from pathlib import Path
import json, tempfile
import numpy as np
import MDAnalysis as mda

ROOT = Path(r"C:\projects\GPCR-virtual-screening")
MD = Path(r"C:\projects\PACER_DC_MD_backup")
TOP = ROOT / "project/results/pacer_dc_long_md_v02/topology"
GRAPH = ROOT / "project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json"

SYSTEMS = (
    "apo",
    "probe_only",
    "compound110__candidate_no_probe",
    "compound110__candidate_probe",
)

AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","HSE":"H","HSD":"H","HSP":"H","ILE":"I","LEU":"L",
    "LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
    "TYR":"Y","VAL":"V",
}

graph = json.loads(GRAPH.read_text(encoding="utf-8"))
graph_seq = "".join(n["site"][0] for n in graph["nodes"])
signatures = {}

print("system                              R  frames  atoms    dt_ps   box  rec_res rec_atoms CA graph finite maxCAjump status")

with tempfile.TemporaryDirectory() as td:
    td = Path(td)

    for system in SYSTEMS:
        src = TOP / system / "minimized.pdb"
        clean = td / f"{system}.pdb"

        with src.open("r", encoding="utf-8", errors="replace") as fi, clean.open("w") as fo:
            for line in fi:
                if not line.startswith("CONECT"):
                    fo.write(line)

        for rep in (1, 2, 3):
            dcd = MD / system / f"replica_{rep:02d}" / "trajectory.dcd"

            try:
                u = mda.Universe(str(clean), str(dcd))
                n = len(u.trajectory)
                dt = float(u.trajectory.dt)

                rec = u.select_atoms("chainID E and protein")
                ca = rec.select_atoms("name CA")

                seq = "".join(AA3.get(r.resname, "?") for r in rec.residues)
                sig = tuple((r.resid, r.resname) for r in rec.residues)
                signatures.setdefault(system, sig)

                finite = True
                boxes = []
                jumps = []

                for fi in sorted(set((0, n // 2, n - 1))):
                    ts = u.trajectory[fi]
                    finite &= bool(np.isfinite(ts.positions).all())
                    boxes.append(
                        ts.dimensions is not None
                        and np.isfinite(ts.dimensions).all()
                        and np.all(ts.dimensions[:3] > 0)
                    )
                    xyz = ca.positions
                    jumps.append(
                        float(np.linalg.norm(np.diff(xyz, axis=0), axis=1).max())
                    )

                status = (
                    len(rec.residues) == 270
                    and len(ca) == 270
                    and len(rec.atoms) == 4349
                    and seq == graph_seq
                    and finite
                    and all(boxes)
                )

                print(
                    f"{system:35s} R{rep} "
                    f"{n:6d} {len(u.atoms):7d} {dt:8.3f} "
                    f"{'YES' if all(boxes) else 'NO ':3s} "
                    f"{len(rec.residues):7d} {len(rec.atoms):9d} {len(ca):3d} "
                    f"{'Y' if seq==graph_seq else 'N'}     "
                    f"{'Y' if finite else 'N'} "
                    f"{max(jumps):9.3f} "
                    f"{'PASS' if status else 'FAIL'}"
                )

            except Exception as e:
                print(f"{system:35s} R{rep} ERROR: {type(e).__name__}: {e}")

print()
base = signatures.get("apo")
for system in SYSTEMS:
    print(f"ORDER {system:35s}: {'MATCH' if signatures.get(system) == base else 'DIFF'}")
