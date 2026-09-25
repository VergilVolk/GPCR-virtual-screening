#!/usr/bin/env python
"""Unified PACER-DC four-context trajectory + embedding extractor.

One entry point, one code path, four contexts. Contract:
``project/pacer_dc_training/FOUR_CONTEXT_WINDOW_CONTRACT.md``.

    candidate × replica × window  ->  {candidate_probe, candidate_no_probe,
                                       probe_only, apo}
    d_PAM = z_CA - z_A        d_AGO = z_C - z_0

Why the extractor has three stages instead of one function: the two resources it
needs live on different machines in this project, and each stage uses only the
dependencies available where its data is.

  ``traj``     (WSL)       MDAnalysis only. Selects the receptor, extracts the
                           matched 100-frame window from each context's
                           production DCD and dumps a compact per-context
                           ``raw/<context>_w<NNN>.npz`` (coordinates + atom and
                           residue names; ~5 MB). The WSL conda env has no
                           dm-tree, so mdgen cannot be imported there.
  ``atom14``   (Windows)   mdgen only. Maps the dumped window onto atom14
                           (CHARMM quirks applied) -> ``atom14/<context>_w<NNN>.npy``
                           + the seqres CSV.
  ``embed``    (Windows)   torch + the nested oneprot-embeddings checkout. Runs
                           the frozen encoder for every context, writes the
                           embeddings, the long-form metadata table, the
                           differential table and (optionally) PCA/tICA/VAMP on
                           the identical window.

``--stage all`` chains them and is only usable on a host that has both.

No stage ever invents a differential: if any of the four contexts is missing or
too short, the unit is written as ``incomplete`` and d_PAM / d_AGO are not
computed.

Run from the repository root.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

CONTEXTS = ("candidate_probe", "candidate_no_probe", "probe_only", "apo")
PREPROCESSING_VERSION = "pacer_dc_atom14_v0.1"
RECEPTOR_SELECTION_DEFAULT = "chainID E and protein"
WINDOW_FRAMES_DEFAULT = 100
FRAME_SPACING_PS_DEFAULT = 10.0

# CHARMM -> atom14 quirks, identical to the G0 converter.
RESNAME_MAP = {"HSE": "HIS", "HSD": "HIS", "HSP": "HIS"}
CHARMM_ATOM_RENAMES = {"ILE": {"CD": "CD1"}}

ROW_FIELDS = [
    "candidate_id", "replica_id", "window_id", "start_frame", "end_frame",
    "frame_count", "context", "system", "trajectory_path", "topology_path",
    "topology_source_path", "topology_note",
    "trajectory_frames_total", "frame_spacing_ps", "receptor_selection",
    "preprocessing_version", "checkpoint_sha256", "raw_path", "atom14_path",
    "embedding_path", "embedding_dim", "finite", "norm", "status",
    "checkpoint_missing_keys", "checkpoint_unexpected_keys",
    "receptor_residues", "receptor_atoms", "ca_neighbour_min_A",
    "ca_neighbour_max_A",
]


def sanitized_topology(path: Path, out_dir: Path, name: str) -> tuple[Path, str]:
    """Return a CONECT-free copy of ``path`` for MDAnalysis.

    The PACER-DC membrane topologies carry hexadecimal CONECT serials
    (e.g. ``CONECT BF2C9 ...``) that MDAnalysis' PDB parser cannot read. We only
    need atom names and coordinates, so every CONECT record is dropped. The
    original file is left untouched; the copy lives beside the unit output.
    """
    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    n_conect = sum(1 for ln in lines if ln.startswith("CONECT"))
    dest_dir = out_dir / "topology"
    dest_dir.mkdir(parents=True, exist_ok=True)
    # Every context needs its own immutable topology artefact.  Reusing
    # ``minimized.pdb`` silently overwrote the previous context and made the
    # provenance rows point to the last processed system.
    dest = dest_dir / f"{name}_{path.name}"
    dest.write_text("\n".join(ln for ln in lines if not ln.startswith("CONECT")) + "\n",
                    encoding="utf-8")
    return dest, f"dropped {n_conect} CONECT records (hex serials break the PDB parser)"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def system_name(candidate: str, context: str) -> str:
    """Production directory name for a context.

    ``probe_only`` and ``apo`` are SHARED controls, so their directory is the
    context name itself (``apo``, ``probe_only``) - not prefixed by the candidate.
    Candidate contexts are ``<candidate>__<context>``.
    """
    if context in ("probe_only", "apo"):
        return context
    return f"{candidate}__{context}"


def dcd_frame_count(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open("rb") as fh:
        head = fh.read(16)
    if len(head) < 12:
        return None
    if head[4:8] == b"CORD":
        return int(np.frombuffer(head[8:12], dtype="<i4")[0])
    if head[:4] == b"CORD":
        return int(np.frombuffer(head[4:8], dtype="<i4")[0])
    return None


def sha256_file(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def load_g0_preprocessing(oneprot_pkg: Path):
    """Import the G0 ``run_g0_embedding`` module so build_batch is shared verbatim."""
    src = oneprot_pkg / "scripts" / "run_g0_embedding.py"
    spec = importlib.util.spec_from_file_location("g0_preproc", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_residue_constants(args):
    if args.residue_constants:
        spec = importlib.util.spec_from_file_location("rc_external", args.residue_constants)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    sys.path.insert(0, str(args.oneprot_root / "external" / "mdgen"))
    from mdgen import residue_constants as rc
    return rc


def base_row(args, context: str, start: int, end: int) -> dict:
    system = system_name(args.candidate, context)
    root = args.work_root
    ref = root / "project/results/pacer_dc_membrane_reference_v01" / system
    prod = root / "project/results/pacer_dc_production_v01" / system / f"replica_{args.replica:02d}"
    topology, trajectory = ref / "minimized.pdb", prod / "trajectory.dcd"
    if args.trajectory_override is not None:
        if args.contexts != ["apo"]:
            raise ValueError("--trajectory-override requires --contexts apo")
        trajectory = args.trajectory_override
    return {f: "" for f in ROW_FIELDS} | {
        "candidate_id": args.candidate, "replica_id": args.replica,
        "window_id": args.window, "start_frame": start, "end_frame": end,
        "context": context, "system": system,
        "trajectory_path": str(trajectory), "topology_path": str(topology),
        "frame_spacing_ps": FRAME_SPACING_PS_DEFAULT,
        "receptor_selection": args.receptor_selection,
        "preprocessing_version": PREPROCESSING_VERSION,
    }


# --------------------------------------------------------------------------- #
# stage 1 (WSL): production DCD -> matched raw window dump
# --------------------------------------------------------------------------- #
def stage_traj(args) -> dict:
    import MDAnalysis as mda

    start, end = args.start_frame, args.start_frame + args.window_frames - 1
    out_raw = args.out_dir / "raw"
    out_raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for context in args.contexts:
        row = base_row(args, context, start, end)
        total = dcd_frame_count(Path(row["trajectory_path"]))
        row["trajectory_frames_total"] = total if total is not None else ""
        if not Path(row["trajectory_path"]).is_file() or not Path(row["topology_path"]).is_file():
            row["status"] = "missing_trajectory"
            rows.append(row)
            print(f"  [{context}] missing trajectory/topology -> incomplete")
            continue
        if total is None or total <= end:
            row["status"] = "insufficient_frames"
            rows.append(row)
            print(f"  [{context}] only {total} frames, need {end + 1} -> incomplete")
            continue

        topo_source = Path(row["topology_path"])
        topo_used, topo_note = sanitized_topology(
            topo_source, args.out_dir, f"{context}_w{args.window:03d}")
        row["topology_source_path"] = str(topo_source)
        row["topology_path"] = str(topo_used)
        row["topology_note"] = topo_note

        u = mda.Universe(str(topo_used), row["trajectory_path"])
        sel = u.select_atoms(args.receptor_selection)
        residues = list(sel.residues)
        atom_names = np.array([a.name for a in sel.atoms])
        resnames = np.array([r.resname for r in residues])
        resids = np.array([r.resid for r in residues])
        res_index = np.concatenate([[k] * len(r.atoms) for k, r in enumerate(residues)])
        coords = np.zeros((args.window_frames, len(sel.atoms), 3), dtype=np.float32)
        for fi, _ts in enumerate(u.trajectory[start:end + 1]):
            coords[fi] = sel.positions      # selection only, not the whole 227k-atom box

        ca = atom_names == "CA"
        d = np.linalg.norm(np.diff(coords[:, ca, :], axis=1), axis=-1)
        raw = out_raw / f"{context}_w{args.window:03d}.npz"
        np.savez(raw, coords=coords, atom_names=atom_names, resnames=resnames,
                 resids=resids, res_index=res_index)
        row.update({
            "status": "ok", "frame_count": args.window_frames,
            "raw_path": f"raw/{raw.name}", "receptor_residues": len(residues),
            "receptor_atoms": int(len(sel.atoms)),
            "ca_neighbour_min_A": round(float(d.min()), 3),
            "ca_neighbour_max_A": round(float(d.max()), 3),
        })
        rows.append(row)
        print(f"  [{context}] {args.window_frames} frames x {len(residues)} residues "
              f"({len(sel.atoms)} atoms) -> {raw.name}; CA neighbour "
              f"{d.min():.2f}-{d.max():.2f} A")
    out = {"rows": rows}
    (args.out_dir / "trajectory_windows.json").write_text(json.dumps(out, indent=2),
                                                          encoding="utf-8")
    return out


# --------------------------------------------------------------------------- #
# stage 2 (Windows): raw window -> atom14 (mdgen mapping)
# --------------------------------------------------------------------------- #
def atom14_index_map(residues, atom_names, rc):
    """(resname, atom_index_array) list -> atom14 slot index map + one-letter seqres."""
    L = len(residues)
    idx = np.full((L, 14), -1, dtype=np.int64)
    letters = []
    for j, (resname, atom_idx) in enumerate(residues):
        std = RESNAME_MAP.get(resname, resname)
        names = rc.restype_name_to_atom14_names.get(std)
        if names is None:
            raise SystemExit(f"non-standard residue {resname} in receptor selection")
        letters.append(rc.restype_3to1[std])
        rename = CHARMM_ATOM_RENAMES.get(std, {})
        local = {rename.get(atom_names[k], atom_names[k]): k for k in atom_idx}
        for slot, nm in enumerate(names):
            if nm and nm in local:
                idx[j, slot] = local[nm]
    return idx, "".join(letters)


def stage_atom14(args) -> dict:
    rc = load_residue_constants(args)
    traj = json.loads((args.out_dir / "trajectory_windows.json").read_text(encoding="utf-8"))
    out_atom14 = args.out_dir / "atom14"
    out_atom14.mkdir(parents=True, exist_ok=True)
    for row in traj["rows"]:
        if row["status"] != "ok":
            continue
        z = np.load(args.out_dir / row["raw_path"])
        coords = z["coords"]
        atom_names = [str(x) for x in z["atom_names"]]
        resnames = [str(x) for x in z["resnames"]]
        res_index = z["res_index"]
        # resnames is per RESIDUE and res_index maps each atom to its residue index,
        # so group atoms by residue index rather than indexing resnames by atom.
        residues = [(resnames[k], np.where(res_index == k)[0])
                    for k in range(int(res_index.max()) + 1)]
        idx, seqres = atom14_index_map(residues, atom_names, rc)
        L = len(residues)
        valid = idx >= 0
        flat = coords.reshape(coords.shape[0], -1, 3)
        arr = np.where(valid[None, :, :, None], flat[:, np.where(valid, idx, 0), :], 0.0)
        arr = arr.astype(np.float32)
        npy = out_atom14 / f"{row['context']}_w{args.window:03d}.npy"
        np.save(npy, arr)
        (out_atom14 / f"{row['context']}_w{args.window:03d}.csv").write_text(
            f"name,seqres\n{args.candidate}_{row['context']},{seqres}\n", encoding="utf-8")
        row["atom14_path"] = f"atom14/{npy.name}"
        row["receptor_residues"] = L
        print(f"  [{row['context']}] atom14 {arr.shape}; seqres {len(seqres)}; "
              f"missing slots/残基 {int((~valid).sum(axis=1).min())}-{int((~valid).sum(axis=1).max())}")
    (args.out_dir / "trajectory_windows.json").write_text(json.dumps(traj, indent=2),
                                                          encoding="utf-8")
    return traj


# --------------------------------------------------------------------------- #
# baselines on the identical window
# --------------------------------------------------------------------------- #
def baseline_representations(ca: np.ndarray, k: int = 5, lag: int = 1) -> dict:
    """PCA / tICA / VAMP on the SAME frames. ca: (T, 3L) -> mean-centred here.

    The representation of a window is the **per-mode amplitude** (standard
    deviation of the projection over the window's frames), not the frame mean of
    the projection: the projection is linear and X is mean-centred, so
    ``mean(X @ V)`` is identically zero and would carry no information about the
    context. Amplitudes are the standard trajectory read-out and are what makes
    the four contexts comparable; the same function is used for every context,
    context differentials are formed from these vectors.
    """
    X = ca - ca.mean(axis=0, keepdims=True)
    T = X.shape[0]
    out = {}
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    out["pca"] = (X @ vt[:k].T).std(axis=0)
    X0, X1 = X[:-lag], X[lag:]
    C00 = X0.T @ X0 / (T - lag) + 1e-6 * np.eye(X.shape[1])
    C11 = X1.T @ X1 / (T - lag) + 1e-6 * np.eye(X.shape[1])
    C01 = X0.T @ X1 / (T - lag)
    w00 = np.linalg.inv(np.linalg.cholesky(C00)).T
    w11 = np.linalg.inv(np.linalg.cholesky(C11)).T
    U, S, Vt = np.linalg.svd(w00 @ C01 @ w11.T, full_matrices=False)
    out["vamp2_score"] = float((S ** 2).sum())
    out["vamp"] = (X @ (w11.T @ Vt[:k].T)).std(axis=0)
    out["tica"] = (X @ (w00.T @ U[:, :k])).std(axis=0)
    return out


# --------------------------------------------------------------------------- #
# stage 3 (Windows): atom14 window -> OneProt-MD embedding
# --------------------------------------------------------------------------- #
def stage_embed(args) -> dict:
    torch = __import__("torch")
    g0 = load_g0_preprocessing(args.oneprot_pkg)
    ckpt_sha = args.checkpoint_sha256 or sha256_file(args.checkpoint)
    out_emb = args.out_dir / "embeddings"
    out_base = args.out_dir / "baselines"
    out_emb.mkdir(parents=True, exist_ok=True)
    out_base.mkdir(parents=True, exist_ok=True)

    traj = json.loads((args.out_dir / "trajectory_windows.json").read_text(encoding="utf-8"))
    obj = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    md = {k: v for k, v in obj["state_dict"].items() if k.startswith("network.md.")}

    # Instantiate and audit one frozen encoder for the whole matched unit.  A
    # separate strict=False load per context could hide missing random weights
    # and contaminate the differential with initialization noise.
    enc = g0.TrajectoryEncoder(
        output_dim=1024, hidden_size=21, num_layers=4, num_heads=8,
        pretrained=False, frozen=True, proj_type="mlp",
        num_frames=args.window_frames, suffix="_i1",
    )
    incompat = enc.load_state_dict(
        {k[len("network.md."):]: v for k, v in md.items()}, strict=False
    )
    missing_keys = list(incompat.missing_keys)
    unexpected_keys = list(incompat.unexpected_keys)
    allowed_unexpected = {"norm.1.log_logit_scale"}
    if missing_keys or set(unexpected_keys) - allowed_unexpected:
        raise RuntimeError(
            "OneProt-MD checkpoint coverage failed: "
            f"missing={missing_keys}, unexpected={unexpected_keys}"
        )
    enc.eval()

    embeddings, seqres_seen, long_rows = {}, set(), []
    for row in traj["rows"]:
        row = {f: row.get(f, "") for f in ROW_FIELDS} | row
        row["checkpoint_sha256"] = ckpt_sha
        if row["status"] != "ok" or not row["atom14_path"]:
            long_rows.append(row)
            continue
        arr = np.load(args.out_dir / row["atom14_path"])
        csv_path = (args.out_dir / row["atom14_path"]).with_suffix(".csv")
        seqres = csv_path.read_text().strip().splitlines()[1].split(",")[1]
        seqres_seen.add(seqres)
        latents, kwargs = g0.build_batch(arr, seqres)
        with torch.no_grad():
            first = enc(latents, 0, **kwargs)
            second = enc(latents, 0, **kwargs)
        vec = first[0].cpu().numpy()
        path = out_emb / f"{row['context']}_w{row['window_id']:03d}.npy"
        np.save(path, vec)
        embeddings[row["context"]] = vec
        row.update({
            "embedding_path": f"embeddings/{path.name}",
            "embedding_dim": int(vec.size),
            "finite": bool(np.isfinite(vec).all()),
            "norm": round(float(np.linalg.norm(vec)), 6),
            "repeat_identical": bool(torch.equal(first, second)),
            "nonzero": int((vec != 0).sum()),
            "latents_shape": "x".join(str(x) for x in latents.shape),
            "checkpoint_missing_keys": json.dumps(missing_keys),
            "checkpoint_unexpected_keys": json.dumps(unexpected_keys),
        })
        long_rows.append(row)
        print(f"  [{row['context']}] z={tuple(vec.shape)} norm={row['norm']} "
              f"finite={row['finite']} repeat_identical={row['repeat_identical']}")
        if args.baselines:
            ca = arr[:, :, 1, :].reshape(arr.shape[0], -1).astype(np.float64)
            reps = baseline_representations(ca)
            np.savez(out_base / f"{row['context']}_w{row['window_id']:03d}.npz", **reps)
            print(f"      baselines: pca{reps['pca'].shape} tica{reps['tica'].shape} "
                  f"vamp{reps['vamp'].shape} vamp2={reps['vamp2_score']:.3f}")

    if len(seqres_seen) > 1:
        raise SystemExit("contexts disagree on seqres - this unit is not matched")

    missing = [c for c in CONTEXTS if c not in embeddings]
    diff = {
        "candidate_id": args.candidate, "replica_id": args.replica,
        "window_id": args.window, "start_frame": args.start_frame,
        "end_frame": args.start_frame + args.window_frames - 1,
        "n_contexts_present": len(embeddings), "missing_contexts": ";".join(missing),
        "completeness": "complete" if not missing else "incomplete",
        "d_PAM_path": "", "d_AGO_path": "", "d_PAM_norm": "", "d_AGO_norm": "",
        "cosine_dPAM_dAGO": "",
        "seqres_length": len(next(iter(seqres_seen))) if seqres_seen else "",
    }
    if not missing:
        d_pam = embeddings["candidate_probe"] - embeddings["probe_only"]
        d_ago = embeddings["candidate_no_probe"] - embeddings["apo"]
        p_pam = out_emb / f"dPAM_w{args.window:03d}.npy"
        p_ago = out_emb / f"dAGO_w{args.window:03d}.npy"
        np.save(p_pam, d_pam)
        np.save(p_ago, d_ago)
        diff.update({
            "d_PAM_path": f"embeddings/{p_pam.name}", "d_AGO_path": f"embeddings/{p_ago.name}",
            "d_PAM_norm": round(float(np.linalg.norm(d_pam)), 6),
            "d_AGO_norm": round(float(np.linalg.norm(d_ago)), 6),
            "cosine_dPAM_dAGO": round(float(np.dot(d_pam, d_ago) /
                                            (np.linalg.norm(d_pam) * np.linalg.norm(d_ago))), 6),
        })
    else:
        print(f"\n  unit INCOMPLETE: missing {missing} -> no d_PAM / d_AGO computed")
    return {
        "long_rows": long_rows,
        "differential": diff,
        "checkpoint_sha256": ckpt_sha,
        "checkpoint_missing_keys": missing_keys,
        "checkpoint_unexpected_keys": unexpected_keys,
    }


def write_tables(out_dir: Path, long_rows, diff) -> None:
    with (out_dir / "four_context_long.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(long_rows)
    with (out_dir / "four_context_differential.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(diff.keys()))
        w.writeheader()
        w.writerow(diff)
    print(f"\nwrote {out_dir / 'four_context_long.csv'}")
    print(f"wrote {out_dir / 'four_context_differential.csv'}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("traj", "atom14", "embed", "all"), default="all")
    p.add_argument("--candidate", default="LY2119620")
    p.add_argument("--replica", type=int, default=1)
    p.add_argument("--window", type=int, default=0, help="window index (start = index*frames)")
    p.add_argument("--window-frames", type=int, default=WINDOW_FRAMES_DEFAULT)
    p.add_argument("--start-frame", type=int, default=None,
                   help="override the derived start frame (default window*window_frames)")
    p.add_argument("--contexts", default=",".join(CONTEXTS),
                   help="comma-separated subset, e.g. apo for a single-context smoke")
    p.add_argument("--work-root", type=Path, default=Path("."),
                   help="repository root holding project/results (the WSL checkout)")
    p.add_argument(
        "--trajectory-override",
        type=Path, default=None,
        help="Explicit DCD path for single-context apo extraction"
    )
    p.add_argument("--oneprot-root", type=Path, default=Path("project/tools/oneprot-embeddings"))
    p.add_argument("--oneprot-pkg", type=Path, default=Path("project/pacer_dc_training/oneprot_g0_audit"))
    p.add_argument("--residue-constants", type=Path, default=None,
                   help="explicit path to mdgen residue_constants.py (hosts without mdgen)")
    p.add_argument("--checkpoint", type=Path,
                   default=Path("project/tools/oneprot-embeddings/artifacts/Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt"))
    p.add_argument("--checkpoint-sha256", default=None)
    p.add_argument("--receptor-selection", default=RECEPTOR_SELECTION_DEFAULT)
    p.add_argument("--baselines", action="store_true",
                   help="also compute PCA/tICA/VAMP on the identical window")
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    if args.baselines:
        raise SystemExit(
            "--baselines is temporarily disabled: the current implementation fits "
            "PCA/tICA/VAMP independently per context, so its axes are not comparable. "
            "Use a frozen train-only shared basis before reporting a baseline."
        )

    args.contexts = [c.strip() for c in args.contexts.split(",") if c.strip()]
    if set(args.contexts) - set(CONTEXTS):
        raise SystemExit(f"unknown contexts: {sorted(set(args.contexts) - set(CONTEXTS))}")
    if args.start_frame is None:
        args.start_frame = args.window * args.window_frames
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"unit: candidate={args.candidate} replica={args.replica} window={args.window} "
          f"frames=[{args.start_frame},{args.start_frame + args.window_frames - 1}]")
    print(f"contexts: {args.contexts}")

    if args.stage in ("traj", "all"):
        print("\n=== stage traj (MDAnalysis) ===")
        stage_traj(args)
    if args.stage in ("atom14", "all"):
        print("\n=== stage atom14 (mdgen mapping) ===")
        stage_atom14(args)
    if args.stage in ("embed", "all"):
        print("\n=== stage embed (OneProt-MD forward) ===")
        res = stage_embed(args)
        write_tables(args.out_dir, res["long_rows"], res["differential"])
        summary = {
            "candidate_id": args.candidate, "replica_id": args.replica,
            "window_id": args.window, "window_frames": args.window_frames,
            "start_frame": args.start_frame,
            "end_frame": args.start_frame + args.window_frames - 1,
            "contexts_requested": args.contexts,
            "contexts_embedded": [r["context"] for r in res["long_rows"] if r["status"] == "ok"],
            "completeness": res["differential"]["completeness"],
            "missing_contexts": res["differential"]["missing_contexts"],
            "checkpoint_sha256": res["checkpoint_sha256"],
            "checkpoint_missing_keys": res["checkpoint_missing_keys"],
            "checkpoint_unexpected_keys": res["checkpoint_unexpected_keys"],
            "receptor_selection": args.receptor_selection,
            "preprocessing_version": PREPROCESSING_VERSION,
            "differential": res["differential"],
        }
        (args.out_dir / "unit_summary.json").write_text(json.dumps(summary, indent=2),
                                                        encoding="utf-8")
        print(f"\nwrote {args.out_dir / 'unit_summary.json'}")
        print(f"COMPLETENESS = {summary['completeness']}")
    print("\nEXTRACTOR_DONE")


if __name__ == "__main__":
    main()
