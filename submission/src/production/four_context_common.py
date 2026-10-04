"""Stage4-only contracts and read-only authentication; no fitting path.

Imports do not require a downloaded trajectory or the historical cache tree.
Historical numeric functions are loaded into a private module instance, so
explicit asset-root/selection bindings never mutate the historical imports.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
VERSION = "PACER_STAGE4_PROSPECTIVE_FKG_V02_v01"
SOURCE_ROOT = Path(os.environ.get("PACER_STAGE4_SOURCE_ROOT", "C:/projects/PACER_STAGE4_MD_backup"))
FROZEN_ROOT = Path(os.environ.get("PACER_STAGE4_FROZEN_ROOT", "C:/projects/GPCR-virtual-screening-fkg-v02"))
RESULT_ROOT = REPO / "project/results/pacer_stage4_prospective_fkg_v02_v01"
CACHE_ROOT = REPO / "project/cache/pacer_stage4_prospective_fkg_v02_v01"
REQUIRED_FREEZE_SHA256 = "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd"
# Phase1/2/3 hashes: final integration frozen-module ledger. Encoder/readout/
# spec hashes: Git blobs of frozen encoder candidate 0a4b4f238c73add75a09a0875f698250a39bdf33.
FROZEN_CODE_SHA256 = {
    "project/pacer_fkg_v02/run_phase1_bs256.py": "d3ed9764ff813b9f164e139668cab2f55d7112d7b04283812d4a5fe3beee2271",
    "project/pacer_fkg_v02/run_phase2_calibration.py": "ccf8af5247ff2c7325f1c451d36943f2e293dee0dff48625031d8bb10a2d6e11",
    "project/pacer_fkg_v02/run_phase3_fkg.py": "147cb31363fe733c72719494dd282c80c86a61e1b733ad2568b94eca3b1d224c",
    "project/encoder_intermediate_v01/core.py": "2fd99945e641682df7de6d9623b344108ee20994723ee032f37312962e479237",
    "project/encoder_atom_readout_v01/experiment_d.py": "4d41168854c28c64f66c9905b17d9090cf1c489acad12ec65760a64457a9e311",
    "project/encoder_candidate_v01/CANDIDATE_SPEC.json": "f1e8337e3b6163a3ec3c9dde255ff4cf214c3305c8115a57fda2fbeba874adc3",
}
FREEZE_REL = "project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json"
MAPPING_REL = "project/results/pacer_dc_four_context_v01/compound110/G2_RESIDUE_MAPPING_v01.csv"
CANDIDATE_CLUSTERS = {"PACER0073": 0, "PACER0027": 4, "PACER0010": 9}
REPLICAS = (1, 2, 3)
SEEDS = {1: 27101, 2: 38201, 3: 49301}
ALIASES = ("A", "P", "C", "CP")
RAW_FRAMES = 1000
RAW_FRAME_SPACING_PS = 10
TEMPORAL_STRIDE = 5
ANALYSIS_FRAMES = 200
ANALYSIS_FRAME_SPACING_PS = 50
BLOCK_FRAMES = 20
BLOCK_DURATION_NS = 1.0  # nominal block interval; 19 internal lags span 950 ps
BLOCKS_PER_TRAJECTORY = 10
FINAL_STEP = 5_000_000
N_RESIDUES = 270
BS_WIDTH = 256
BRANCHES = ("STATE_MOTION", "SIGNED_DRIFT")
TEMPORAL_CONTRACT = {
    "raw_frames": RAW_FRAMES, "raw_frame_spacing_ps": RAW_FRAME_SPACING_PS,
    "stride": TEMPORAL_STRIDE, "stride_definition": "frames[::5]",
    "raw_frame_ids": list(range(0, RAW_FRAMES, TEMPORAL_STRIDE)),
    "analysis_frames": ANALYSIS_FRAMES, "analysis_frame_spacing_ps": ANALYSIS_FRAME_SPACING_PS,
    "block_frames": BLOCK_FRAMES, "blocks_per_trajectory": BLOCKS_PER_TRAJECTORY,
    "block_duration_ns": BLOCK_DURATION_NS, "internal_lag_ps": 50,
    "first_raw_frame_index": 0, "last_raw_frame_index": 995,
    "block_internal_span_ps": 950, "interpolation": False,
}
FORBIDDEN = {name: False for name in (
    "training", "calibration", "normalization_refit", "bandwidth_refit", "rff_refit",
    "graph_refit", "region_refit", "threshold_selection", "outcome_driven_tuning",
    "md_rerun", "source_modification",
)}
CLAIM = ("Prospective frozen four-context dynamic differential evaluation only. "
         "No PAM/ago-PAM labels, probability, efficacy or independent-block inference.")
AA = dict(zip(("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL").split(),
              "ARNDCQEGHILKMFPSTWYV"))
AA.update({"HSE": "H", "HSD": "H", "HSP": "H"})
RESNAME_ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CYM": "CYS"}
AA.update({name: AA[standard] for name, standard in RESNAME_ALIASES.items()})


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def artifact(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size}


def verify_artifact(record):
    p = Path(record["path"])
    if not p.is_file() or p.stat().st_size != record["bytes"] or sha256(p) != record["sha256"]:
        raise RuntimeError(f"artifact hash/size mismatch: {p}")
    return p


def guard_output(path):
    p = Path(path).resolve()
    # Check fixed logical roots and their resolved counterparts; junctions may
    # not redirect Stage4 output into a historical result/cache tree.
    for rel in ("project/results/pacer_stage4_prospective_fkg_v02_v01",
                "project/cache/pacer_stage4_prospective_fkg_v02_v01"):
        logical = REPO.resolve() / rel
        if logical.resolve() != logical:
            raise ValueError(f"Stage4 output root redirects through a junction: {logical}")
        if p.is_relative_to(logical):
            return p
    raise ValueError(f"outside dedicated Stage4 output roots: {p}")


def write_bytes(path, data):
    """Immutable publication: valid identical output can resume; never replace."""
    p = guard_output(path)
    if p.exists():
        if p.read_bytes() != data:
            raise FileExistsError(f"existing output differs; refusing overwrite: {p}")
        return
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = guard_output(p.with_name(f".{p.name}.partial.{os.getpid()}"))
    if tmp.exists():
        raise FileExistsError(tmp)
    try:
        with tmp.open("xb") as f:
            f.write(data)
        # Hard-link publication is atomic and fails if another writer won.
        os.link(tmp, p)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_json(path, value):
    write_bytes(path, (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode())


def write_array(path, array, compressed=False):
    import io
    import numpy as np
    b = io.BytesIO()
    if compressed:
        np.savez_compressed(b, **array)
    else:
        np.save(b, array, allow_pickle=False)
    write_bytes(path, b.getvalue())


def code_identity():
    return [artifact(Path(__file__).with_name(name)) for name in (
        "stage4_prospective_common_v01.py", "run_stage4_prospective_phase1_bs256_v01.py",
        "run_stage4_prospective_phase2a_frozen_apply_v01.py",
        "run_stage4_prospective_phase2b_graph_region_v01.py")]


def private_engine():
    """Read-only historical functions with an explicit asset-root binding."""
    path = REPO / "project/pacer_fkg_v02/run_phase3_fkg.py"
    spec = importlib.util.spec_from_file_location("_stage4_private_frozen_phase3", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.REPO = FROZEN_ROOT.resolve()
    return module


def verify_frozen_anchor():
    """Authenticate serialized historical state, without calibration replay.

    The old verify_frozen_anchor invokes historical Phase2.verify, which
    reconstructs fitted arrays and requires every historical MD cache. Stage4
    instead authenticates all consumed bytes against the SAME immutable SHA.
    Numeric loading/mapping/diffusion remain the unmodified historical functions.
    """
    import numpy as np
    code_records = []
    for relative, expected in FROZEN_CODE_SHA256.items():
        path = REPO / relative
        actual = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        if actual != expected:
            raise RuntimeError(f"pinned historical code differs: {path}")
        code_records.append({"path": str(path), "lf_sha256": actual})
    fp = FROZEN_ROOT / FREEZE_REL
    if not fp.is_file() or sha256(fp) != REQUIRED_FREEZE_SHA256:
        raise RuntimeError(f"required historical freeze unavailable or hash mismatch: {fp}")
    freeze = load_json(fp)
    if freeze.get("V02_NUMERICAL_STATE_FROZEN") is not True:
        raise RuntimeError("historical numerical state is not frozen")
    records = list(freeze["definitions"].values()) + [freeze["calibration_population"]]
    state = {x["branch"]: x for x in freeze["fitted_state"]}
    engine = private_engine()
    if engine.REQUIRED_FREEZE_SHA256 != REQUIRED_FREEZE_SHA256 or set(state) != set(BRANCHES):
        raise RuntimeError("historical engine/state inventory mismatch")
    for branch, item in state.items():
        records += [item["normalization"]["center"], item["normalization"]["scale"],
                    item["bandwidth"], item["rff"]["weights"], item["rff"]["bias"]]
        expected = engine.BRANCHES[branch]
        if item["input_width"] != expected["width"] or item["rff_features"] != 512 or item["rff"]["seed"] != expected["seed"]:
            raise RuntimeError("historical width/RFF seed drift")
    consumed = []
    for rec in records:
        p = (FROZEN_ROOT / rec["path"]).resolve()
        if not p.is_relative_to(FROZEN_ROOT.resolve()):
            raise RuntimeError("historical manifest path escape")
        verify_artifact({**rec, "path": str(p)})
        consumed.append({**rec, "path": str(p)})
    # Source-code provenance is LF-normalised only for text: no scientific
    # JSON/NPY is normalised. Windows checkout EOL does not alter code identity.
    for rec in freeze["source_code"]:
        p = REPO / rec["path"]
        raw = p.read_bytes().replace(b"\r\n", b"\n")
        variants = (raw, raw.replace(b"\n", b"\r\n"))
        if not any(hashlib.sha256(v).hexdigest() == rec["sha256"] and len(v) == rec["bytes"] for v in variants):
            raise RuntimeError(f"frozen implementation source mismatch: {p}")
    inventory = load_json(FROZEN_ROOT / freeze["calibration_population"]["path"])
    fit = [x for x in inventory["trajectories"] if x["used_for_fitting"]]
    if len(fit) != 4 or any(x["replica"] != 2 for x in fit):
        raise RuntimeError("historical calibration provenance is not R2-only")
    anchor = {"freeze": freeze, "state": state}
    for branch in BRANCHES:
        loaded = engine.load_frozen_state(anchor, branch)
        if not all(np.isfinite(loaded[k]).all() for k in ("center", "scale", "weights", "bias")):
            raise RuntimeError("non-finite frozen numeric state")
        if load_json(FROZEN_ROOT / state[branch]["bandwidth"]["path"])["bandwidth"] != engine.BRANCHES[branch]["bandwidth"]:
            raise RuntimeError("historical bandwidth drift")
    mapping = FROZEN_ROOT / MAPPING_REL
    # Mapping is verified against the graph's ordered identity, not rebuilt.
    with mapping.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    graph = load_json(FROZEN_ROOT / freeze["definitions"]["graph_definitions"]["path"])
    if len(rows) != N_RESIDUES or [int(x["embedding_index"]) for x in rows] != list(range(N_RESIDUES)):
        raise RuntimeError("frozen mapping order/count mismatch")
    nodes = graph["nodes"]
    for row, node in zip(rows, nodes):
        if (int(node["embedding_index"]) != int(row["embedding_index"])
                or node["site"] != row["aa"] + row["structure_resid"]
                or AA.get(row["resname"]) != row["aa"]):
            raise RuntimeError("mapping does not match frozen graph")
    anchor["receipt"] = {
        "freeze": artifact(fp), "consumed_assets": consumed, "mapping": artifact(mapping),
        "frozen_implementation": code_records,
        "verification_policy": "serialized-state hashes; no calibration replay/refit",
        "historical_calibration_replicas": [2], "prospective_replicas": "R1/R2/R3 evaluation only",
    }
    return engine, anchor


def context_map(candidate):
    cluster = CANDIDATE_CLUSTERS[candidate]
    return {"A": f"cluster{cluster}__apo", "P": f"cluster{cluster}__probe_only",
            "C": f"{candidate}__candidate_no_probe", "CP": f"{candidate}__candidate_probe"}


@dataclass(frozen=True)
class Job:
    candidate: str
    cluster: int
    context: str
    system: str
    replica: int
    seed: int
    topology: Path | None = None
    trajectory: Path | None = None
    sequence: str = ""
    selection: str = ""
    input_record: dict | None = None

    @property
    def key(self):
        return f"{self.system}__replica_{self.replica:02d}"

    def identity(self):
        return {k: getattr(self, k) for k in ("candidate", "cluster", "context", "system", "replica", "seed")}


def expected_jobs():
    return [Job(c, cluster, a, s, r, SEEDS[r])
            for c, cluster in CANDIDATE_CLUSTERS.items()
            for a, s in context_map(c).items() for r in REPLICAS]


def discover_inputs(job, source=SOURCE_ROOT):
    root = Path(source)
    top_dir = root / "systems" / job.system
    matches = sorted(top_dir.rglob("minimized.pdb")) if top_dir.exists() else []
    if len(matches) != 1:
        raise RuntimeError(f"expected one topology systems/{job.system}/**/minimized.pdb; found {len(matches)}")
    prod = root / "production" / job.system / f"replica_{job.replica:02d}"
    paths = {"topology": matches[0], "trajectory": prod / "trajectory.dcd",
             "state": prod / "state.csv", "progress": prod / "progress.json"}
    for k in ("topology", "trajectory", "state"):
        if not paths[k].is_file():
            raise FileNotFoundError(paths[k])
    return paths


def completion_acceptance(job, progress, final_step, frames, finite, topology_compatible):
    if final_step != FINAL_STEP or frames != RAW_FRAMES or not finite or not topology_compatible:
        raise RuntimeError(f"physical trajectory incomplete/incompatible: {job.key}")
    if progress is None:
        raise RuntimeError(f"missing progress identity/seed record: {job.key}")
    if progress.get("system") != job.system or progress.get("replica") != job.replica or progress.get("seed") != job.seed:
        raise RuntimeError(f"progress system/replica/seed mismatch: {job.key}")
    if progress.get("target_ns") != 10.0 or progress.get("timestep_fs") != 2.0 or progress.get("restraint_k_kj_mol_nm2") != 0.0:
        raise RuntimeError(f"production protocol mismatch: {job.key}")
    if progress.get("status") == "complete" and progress.get("completed_ns") == 10.0:
        return "physical_complete_and_progress_complete"
    known = (job.candidate == "PACER0010" and job.context == "CP" and job.replica == 3)
    if known and progress.get("status") == "running" and progress.get("completed_ns") == 9.9:
        return "KNOWN_STALE_PROGRESS_EXCEPTION_PACER0010_CP_R3_PHYSICAL_COMPLETE"
    raise RuntimeError(f"unexpected progress inconsistency: {job.key}")


def read_state(path):
    import numpy as np
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        rows = [{k.lstrip('#').strip('"'): v for k, v in r.items()} for r in reader]
    if len(rows) != RAW_FRAMES:
        raise RuntimeError(f"state.csv rows {len(rows)} != {RAW_FRAMES}: {path}")
    steps = [int(x["Step"]) for x in rows]
    if steps != list(range(5000, FINAL_STEP + 1, 5000)):
        raise RuntimeError(f"state.csv step cadence mismatch: {path}")
    times = np.asarray([float(x["Time (ps)"]) for x in rows])
    if not np.isfinite(times).all() or not np.allclose(np.diff(times), RAW_FRAME_SPACING_PS, rtol=0, atol=1e-4):
        raise RuntimeError(f"state.csv physical time spacing mismatch: {path}")
    for name in ("Potential Energy (kJ/mole)", "Kinetic Energy (kJ/mole)", "Temperature (K)", "Box Volume (nm^3)", "Density (g/mL)"):
        if not np.isfinite([float(x[name]) for x in rows]).all():
            raise RuntimeError(f"state.csv non-finite {name}: {path}")
    return {"last_step": steps[-1], "rows": len(rows), "first_time_ps": float(times[0]),
            "last_time_ps": float(times[-1]), "time_spacing_ps": RAW_FRAME_SPACING_PS}


def sequence_projection(source_sequence, target_sequence):
    """Unique ordered subsequence only; never a guessed alignment/reindexing."""
    # Count paths up to two without recursion: multiple receptor copies must
    # yield the explicit ambiguity error even for long protein topologies.
    n, m = len(source_sequence), len(target_sequence)
    counts = [[0] * m + [1] for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            counts[i][j] = min(2, counts[i + 1][j] +
                               (counts[i + 1][j + 1] if source_sequence[i] == target_sequence[j] else 0))
    if counts[0][0] != 1:
        raise RuntimeError("receptor sequence projection is missing or ambiguous; no mapping guessed")
    ids, j = [], 0
    for i in range(n):
        if j == m:
            break
        if source_sequence[i] == target_sequence[j] and counts[i + 1][j + 1]:
            ids.append(i)
            j += 1
    return tuple(ids)


def topology_selection(universe, sequence):
    """Project actual construct into the unchanged frozen 270-node order."""
    if len(sequence) != N_RESIDUES:
        raise RuntimeError("target sequence is not the frozen 270-residue contract")
    # OpenMM preserves construct order but splits it into independently capped
    # chains at the missing loop. Chain IDs, segids and restarted PDB resids are
    # provenance, not embedding identities. Project ALL standard protein
    # residues in topology order; a second possible projection must fail closed.
    residues = [r for r in universe.select_atoms("protein").residues if r.resname in AA]
    ids = sequence_projection("".join(AA[r.resname] for r in residues), sequence)
    selected = [residues[i] for i in ids]
    selection = "resindex " + " ".join(str(r.ix) for r in selected)
    mapping = [{"embedding_index": i, "source_resindex": int(r.ix), "source_resid": int(r.resid),
                "source_chain": "|".join(sorted(set(r.atoms.chainIDs))),
                "source_segid": r.segid, "resname": r.resname} for i, r in enumerate(selected)]
    return selection, {"receptor_residue_count": N_RESIDUES, "construct_standard_residues": len(residues),
                       "excluded_construct_resindices": [int(r.ix) for i, r in enumerate(residues) if i not in ids],
                       "sequence": sequence, "residue_order_compatible": True,
                       "policy": "unique topology-ordered sequence subsequence across construct fragments into fixed embedding order; no coordinate changes",
                       "mapping": mapping}


def sanitize_topology(source, system):
    # Original topology and atom order remain untouched; only CONECT records
    # are dropped from a separate parser input, as in the Stage-B template.
    data = b"".join(x for x in Path(source).read_bytes().splitlines(keepends=True) if not x.startswith(b"CONECT"))
    p = RESULT_ROOT / "input_audit/topology_sanitized" / f"{system}_minimized_no_conect.pdb"
    write_bytes(p, data)
    return p


def validate_selected_atoms(universe, selection, sequence):
    from project.pacer_dc_training.extract_geom2vec_atom14 import ATOM14
    residues = universe.select_atoms(selection).residues
    if len(residues) != N_RESIDUES:
        raise RuntimeError("selected receptor is not the frozen 270-residue order")
    total = 0
    for i, (r, aa) in enumerate(zip(residues, sequence)):
        names = ["CD1" if aa == "I" and atom.name == "CD" else atom.name for atom in r.atoms]
        if len(names) != len(set(names)):
            raise RuntimeError(f"duplicate receptor atom names at frozen residue {i}")
        missing = set(ATOM14[aa]) - set(names)
        if missing:
            raise RuntimeError(f"missing atom14 at frozen residue {i}: {sorted(missing)}")
        total += len(ATOM14[aa])
    return total


def input_availability(source=SOURCE_ROOT):
    records = []
    for job in expected_jobs():
        try:
            paths = discover_inputs(job, source)
            records.append({**job.identity(), "files_present": True, "paths": {k: str(v) for k, v in paths.items()}})
        except (RuntimeError, FileNotFoundError) as e:
            records.append({**job.identity(), "files_present": False, "reason": str(e)})
    return {"source_root": str(Path(source).resolve()), "expected_trajectories": 36,
            "file_sets_present": sum(x["files_present"] for x in records), "jobs": records}


def authenticate_inputs(source=SOURCE_ROOT):
    import numpy as np
    import MDAnalysis as mda
    _, anchor = verify_frozen_anchor()
    with (FROZEN_ROOT / MAPPING_REL).open(encoding="utf-8-sig", newline="") as f:
        seq = "".join(AA[x["resname"]] for x in csv.DictReader(f))
    availability = input_availability(source)
    if availability["file_sets_present"] != 36:
        missing = next(x["reason"] for x in availability["jobs"] if not x["files_present"])
        raise RuntimeError(f"Stage4 download not ready: {availability['file_sets_present']}/36 input sets; {missing}")
    jobs = []
    # The MD job matrix independently authenticates candidate/cluster/seed
    # provenance. Its bookkeeping status need not be a completion indicator.
    matrix = Path(source) / "job_matrix.csv"
    with matrix.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    keyed = {(r["system"], int(r["replica"])): r for r in rows}
    if len(rows) != 36 or len(keyed) != 36:
        raise RuntimeError("Stage4 job matrix must contain 36 unique jobs")
    for i, job in enumerate(expected_jobs(), 1):
        paths = discover_inputs(job, source)
        row = keyed.get((job.system, job.replica), {})
        expected_context = {"A": "apo", "P": "probe_only", "C": "candidate_no_probe", "CP": "candidate_probe"}[job.context]
        if (row.get("candidate_id") != job.candidate or int(row.get("receptor_cluster", -1)) != job.cluster
                or row.get("context") != expected_context
                or int(row.get("seed", -1)) != job.seed or int(row.get("production_steps", -1)) != FINAL_STEP
                or float(row.get("production_length_ns", -1)) != 10):
            raise RuntimeError(f"job matrix mismatch: {job.key}")
        print(f"authenticate {i}/36 {job.key}", flush=True)
        before = {k: (p.stat().st_size, p.stat().st_mtime_ns) for k, p in paths.items() if p.is_file()}
        state = read_state(paths["state"])
        progress = load_json(paths["progress"]) if paths["progress"].is_file() else None
        top = sanitize_topology(paths["topology"], job.system)
        u = mda.Universe(str(top), str(paths["trajectory"]))
        try:
            selection, mapping = topology_selection(u, seq)
            mapping["atom14_heavy_atoms"] = validate_selected_atoms(u, selection, seq)
            if len(u.trajectory) != RAW_FRAMES or not np.isclose(u.trajectory.dt, RAW_FRAME_SPACING_PS, rtol=0, atol=1e-3):
                raise RuntimeError(f"DCD frame count/physical dt mismatch: {job.key}")
            finite = True
            for ts in u.trajectory:
                if not np.isfinite(ts.positions).all():
                    finite = False
                    break
            reason = completion_acceptance(job, progress, state["last_step"], len(u.trajectory), finite, True)
        finally:
            u.trajectory.close()
        rec = {**job.identity(), "trajectory": artifact(paths["trajectory"]),
               "topology": artifact(paths["topology"]), "sanitized_topology": artifact(top),
               "state": artifact(paths["state"]), "state_final": state,
               "progress": artifact(paths["progress"]), "progress_status": progress["status"],
               "progress_completed_ns": progress["completed_ns"], "acceptance_reason": reason,
               "coordinates_finite": finite, "dcd_frames": RAW_FRAMES, "dcd_spacing_ps": RAW_FRAME_SPACING_PS,
               "receptor": mapping, "selection": selection}
        after = {k: (p.stat().st_size, p.stat().st_mtime_ns) for k, p in paths.items() if p.is_file()}
        if before != after:
            raise RuntimeError(f"source changed while authenticating (download still active): {job.key}")
        jobs.append(replace(job, topology=top, trajectory=paths["trajectory"], sequence=seq, selection=selection, input_record=rec))
    receipt = {"schema": "pacer.stage4.input_authentication.v1", "status": "36_OF_36_ACCEPTED",
               "source_root": str(Path(source).resolve()), "job_matrix": artifact(matrix),
               "jobs": [x.input_record for x in jobs], "historical_anchor": anchor["receipt"],
               "temporal_contract": TEMPORAL_CONTRACT, "implementation": code_identity(),
               "forbidden_operations_performed": FORBIDDEN}
    write_json(RESULT_ROOT / "input_audit/INPUT_AUTHENTICATION_v01.json", receipt)
    return jobs, receipt


def parser(description):
    p = argparse.ArgumentParser(description=description)
    g = p.add_mutually_exclusive_group(required=True)
    for name in ("smoke", "run", "verify", "preflight"):
        g.add_argument(f"--{name}", action="store_true")
    p.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    return p


def preflight(source=SOURCE_ROOT):
    _, anchor = verify_frozen_anchor()
    return {"historical_freeze_verified": True, "historical_anchor": anchor["receipt"],
            "mapping": {"candidates": 3, "systems": 12, "trajectories": 36},
            "temporal_contract": TEMPORAL_CONTRACT, "availability": input_availability(source)}
