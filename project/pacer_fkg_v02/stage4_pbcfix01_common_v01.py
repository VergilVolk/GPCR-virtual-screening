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
import functools
import threading
from contextlib import contextmanager
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
VERSION = "PACER_STAGE4_PBCFIX01_FKG_V02_v01"
SOURCE_ROOT = Path(os.environ.get("PACER_STAGE4_PBCFIX01_SOURCE_ROOT", "C:/projects/PACER_STAGE4_PBCFIX01_MD_backup"))
FROZEN_ROOT = Path(os.environ.get("PACER_STAGE4_FROZEN_ROOT", "C:/projects/GPCR-virtual-screening-fkg-v02"))
RESULT_ROOT = REPO / "project/results/pacer_stage4_pbcfix01_fkg_v02_v01"
CACHE_ROOT = REPO / "project/cache/pacer_stage4_pbcfix01_fkg_v02_v01"
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
CANDIDATE_CLUSTERS = {"PACERGEN02432": 1, "PACERGEN00123": 1, "PACERGEN00462": 1}
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
CLAIM = ("Prospective frozen PACER-FKG dynamic differential evaluation only. "
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
    for rel in ("project/results/pacer_stage4_pbcfix01_fkg_v02_v01",
                "project/cache/pacer_stage4_pbcfix01_fkg_v02_v01"):
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
        "stage4_pbcfix01_common_v01.py", "run_stage4_pbcfix01_inputs_v01.py", "run_stage4_pbcfix01_phase1_bs256_v01.py",
        "run_stage4_pbcfix01_phase2a_frozen_apply_v01.py",
        "run_stage4_pbcfix01_phase2b_graph_region_v01.py", "stage4_pbcfix01_deployment_v01.py")]


def private_engine():
    """Read-only historical functions with an explicit asset-root binding."""
    path = REPO / "project/pacer_fkg_v02/run_phase3_fkg.py"
    spec = importlib.util.spec_from_file_location("_stage4_pbcfix01_private_frozen_phase3", path)
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


PROTOCOL = "pacer_stage4_rerun_v01_pbcfix01"
SERVER_ROOT = f"/data/pacer_stage4_rerun_v01_deployment/{PROTOCOL}"
CONTEXT_NAMES = dict(zip(ALIASES, ("apo", "probe_only", "candidate_no_probe", "candidate_probe")))
SEGMENT_STEPS = 50_000
SEGMENT_ENDS = tuple(range(SEGMENT_STEPS, FINAL_STEP + 1, SEGMENT_STEPS))
TEMPORAL_CONTRACT.update(physical_first_time_ps=10, physical_last_time_ps=10000,
                         analysis_first_time_ps=10, analysis_last_time_ps=9960,
                         segments=100, frames_per_segment=10,
                         dcd_time_policy="receipt/state validated; chain frame i has time (i+1)*10 ps")


def context_map(candidate):
    cluster = CANDIDATE_CLUSTERS[candidate]
    return {a: f"{PROTOCOL}__{candidate}__cluster_{cluster:02d}__{name}"
            for a, name in CONTEXT_NAMES.items()}


@dataclass(frozen=True)
class Job:
    candidate: str
    cluster: int
    context: str
    system: str
    replica: int
    seed: int
    topology: Path | None = None
    trajectory: tuple[Path, ...] | None = None
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


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def member(root, relative):
    """Only portable, contained campaign paths; never follow an escaping link."""
    root = Path(root).resolve()
    name = str(relative)
    require(bool(name) and "\\" not in name and ":" not in name
            and not Path(name).is_absolute() and ".." not in Path(name).parts,
            f"unsafe campaign path: {name}")
    path = (root / name).resolve()
    require(path.is_relative_to(root) and path != root, f"campaign path escape: {name}")
    return path


def campaign_relative(value):
    """Translate only the declared server root; retain raw receipts unchanged."""
    value = str(value)
    require("\\" not in value and ":" not in value and ".." not in PurePosixPath(value).parts,
            f"unsafe server campaign path: {value}")
    if value.startswith("/"):
        require(value.startswith(SERVER_ROOT + "/"), f"unrecognized server root: {value}")
        value = value[len(SERVER_ROOT) + 1:]
    require(bool(value) and value != "." and not value.startswith("/"), "empty campaign member")
    return value


def input_scope(source):
    return hashlib.sha256(os.path.normcase(str(Path(source).resolve())).encode()).hexdigest()[:12]


_LOCK_STATE = threading.local()


@contextmanager
def analysis_lock():
    """OS advisory lock releases on process death; nested calls are reentrant."""
    if getattr(_LOCK_STATE, "held", False):
        yield
        return
    path = guard_output(RESULT_ROOT / ".analysis.lock")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as handle:
        if path.stat().st_size == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError("another PBCFIX01 analysis/recovery process is active") from exc
        _LOCK_STATE.held = True
        try:
            yield
        finally:
            _LOCK_STATE.held = False
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def exclusive_analysis(function):
    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        with analysis_lock():
            return function(*args, **kwargs)
    return wrapped


def quarantine_outputs(paths):
    """Preserve failed artifacts in this analysis version, with no overwrites."""
    import uuid
    token = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "_" + uuid.uuid4().hex
    plan = []
    for original in paths:
        original = guard_output(original)
        if not original.exists():
            continue
        root = next(r.resolve() for r in (RESULT_ROOT, CACHE_ROOT) if original.is_relative_to(r.resolve()))
        require(original != root and not original.is_symlink()
                and original.name != ".analysis.lock" and "_recovery" not in original.relative_to(root).parts,
                "unsafe recovery target")
        target = guard_output(root / "_recovery" / token / original.relative_to(root))
        require(not target.exists(), "recovery target already exists")
        plan.append((original, target))
    # All target paths have been resolved and checked before any move.
    for original, target in plan:
        target.parent.mkdir(parents=True, exist_ok=True)
        os.rename(original, target)
    receipt = {"status": "FAILED_ANALYSIS_ARTIFACTS_PRESERVED", "moves":
               [{"original": str(p), "preserved": str(q)} for p, q in plan]}
    write_json(RESULT_ROOT / "_recovery" / token / "RECOVERY_v01.json", receipt)
    return receipt


def csv_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def verify_sums(root, filename="SHA256SUMS_PBCFIX01.txt"):
    manifest = member(root, filename)
    lines = manifest.read_text(encoding="utf-8-sig").splitlines()
    require(bool(lines), "empty package integrity receipt")
    seen, records = set(), []
    for line in lines:
        digest, relative = line.split("  ", 1)
        require(relative not in seen and len(digest) == 64, "duplicate/invalid integrity entry")
        seen.add(relative)
        path = member(root, relative)
        require(sha256(path) == digest, f"package integrity mismatch: {relative}")
        records.append(artifact(path))
    for name in ("manifests/STAGE4_RERUN_JOB_MATRIX.csv", "manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv",
                 f"config/{PROTOCOL}.json", f"scripts/run_{PROTOCOL}.py"):
        require(name in seen, f"required package member not authenticated: {name}")
    return [artifact(manifest), *records]


def job_row(job):
    return {"job_id": f"{job.system}__R{job.replica}", "system_id": job.system,
            "candidate_id": job.candidate, "context": CONTEXT_NAMES[job.context],
            "replica": str(job.replica), "seed": str(job.seed), "production_ns": "10",
            "production_steps": str(FINAL_STEP), "timestep_fs": "2", "trajectory_ps": "10",
            "output_path": f"production/{job.system}/replica_{job.replica:02d}",
            "checkpoint_path": f"production/{job.system}/replica_{job.replica:02d}/segments/<end_step>/checkpoint.chk"}


def validate_matrix(source):
    """Validate all 36 identities even before coordinates have arrived."""
    source = Path(source).resolve()
    for output in (RESULT_ROOT.resolve(), CACHE_ROOT.resolve()):
        require(not source.is_relative_to(output) and not output.is_relative_to(source),
                "source and dedicated analysis output roots overlap")
    rows = csv_rows(source / "manifests/STAGE4_RERUN_JOB_MATRIX.csv")
    keyed = {r["job_id"]: r for r in rows}
    require(len(rows) == len(keyed) == 36, "job matrix must contain 36 unique jobs")
    for job in expected_jobs():
        expected = job_row(job)
        row = keyed.get(expected["job_id"], {})
        require(all((campaign_relative(row.get(k, "")) if k in ("output_path", "checkpoint_path") else row.get(k)) == v
                    for k, v in expected.items()), f"matrix identity/protocol mismatch: {job.key}")
        require(row.get("source_job_id") == expected["job_id"].replace(PROTOCOL, "pacer_stage4_rerun_v01"),
                f"wrong source job identity: {job.key}")
    systems = csv_rows(source / "manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv")
    by_system = {r["system_id"]: r for r in systems}
    require(len(systems) == len(by_system) == 12, "expected 12 independent systems")
    for job in expected_jobs()[::3]:
        row = by_system.get(job.system, {})
        expected = {"candidate_id": job.candidate, "cluster": "cluster_01", "alias": job.context,
                    "context": CONTEXT_NAMES[job.context], "system_path": f"systems/{job.system}",
                    "has_candidate": str(job.context in ("C", "CP")), "has_ACH": str(job.context in ("P", "CP"))}
        require(all((campaign_relative(row.get(k, "")) if k == "system_path" else row.get(k)) == v
                    for k, v in expected.items()), f"system manifest mismatch: {job.system}")
    cfg = load_json(source / f"config/{PROTOCOL}.json")
    expected = {"protocol_id": PROTOCOL, "production_steps": FINAL_STEP, "production_ns": 10,
                "timestep_fs": 2, "trajectory_ps": 10, "checkpoint_ps": 100,
                "seeds": list(SEEDS.values()), "restraint_expression": "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2"}
    require(all(cfg.get(k) == v for k, v in expected.items()), "production configuration drift")
    return keyed


def discover_inputs(job, source=SOURCE_ROOT, segment_count=100):
    root = Path(source)
    base = member(root, job_row(job)["output_path"])
    segments = [member(base, f"segments/{end:08d}") for end in SEGMENT_ENDS[:segment_count]]
    paths = {"topology": member(root, f"systems/{job.system}/minimized.pdb"),
             "progress": member(base, "progress.json"),
             "trajectory": tuple(member(d, "trajectory.dcd") for d in segments),
             "segments": segments}
    for p in [paths["topology"], paths["progress"], *paths["trajectory"],
              *(member(d, name) for d in segments for name in ("state.csv", "checkpoint.chk", "receipt.json"))]:
        if not p.is_file():
            raise FileNotFoundError(p)
    # Ignore archived .interrupted.* directories; reject additional canonical segments.
    found = {p.name for p in (base / "segments").iterdir() if p.is_dir() and p.name.isdigit()}
    require(found <= {f"{end:08d}" for end in SEGMENT_ENDS}, "unexpected canonical segment directory")
    return paths


def authenticated_outputs(directory, outputs, records):
    require(bool(outputs), f"empty outputs receipt: {directory}")
    for name, digest in outputs.items():
        path = member(directory, name)
        record = artifact(path)
        require(record["bytes"] > 0, f"empty committed output: {path}")
        require(record["sha256"] == digest, f"output SHA256 mismatch: {path}")
        records.append(record)


def production_identity(job, source, row):
    """Reconstruct the production fingerprint without loading OpenMM/checkpoints."""
    source = Path(source)
    records = []
    def read(relative):
        p = member(source, relative)
        record = artifact(p)
        obj = load_json(p)
        verify_artifact(record)
        records.append(record)
        return obj
    def digest(relative):
        rec = artifact(member(source, relative))
        records.append(rec)
        return rec["sha256"]
    sid = job.system
    build = read(f"systems/{sid}/build_audit.json")
    require(build.get("status") == "COMPLETE" and build.get("protocol_id") == PROTOCOL,
            "incomplete or wrong membrane build")
    authenticated_outputs(source / "systems" / sid, build["outputs"], records)
    require("minimized.pdb" in build["outputs"] and "system.xml" in build["outputs"] and "state.xml" in build["outputs"],
            "build receipt lacks topology/system/state")
    gate_rel = f"logs/server_preflight_{PROTOCOL}.json"
    gate = read(gate_rel)
    config_sha = digest(f"config/{PROTOCOL}.json")
    require(gate.get("status") == "PASS" and gate["config_sha256"] == config_sha
            and gate["builds"][sid] == digest(f"systems/{sid}/build_audit.json"), "production gate drift")
    approval = read("authorization.json")
    require(approval.get("protocol_id") == PROTOCOL and approval["preflight_sha256"] == digest(gate_rel)
            and approval.get("scientific_review_approved") is True
            and approval.get("user_explicit_authorization") is True, "production authorization receipt mismatch")
    authorization = dict(protocol_id=PROTOCOL, system_id=sid,
                         build_audit_sha256=digest(f"systems/{sid}/build_audit.json"),
                         config_sha256=config_sha, package_manifest_sha256=digest("SHA256SUMS_PBCFIX01.txt"),
                         system_sha256=digest(f"systems/{sid}/system.xml"),
                         topology_sha256=digest(f"systems/{sid}/minimized.pdb"),
                         state_sha256=digest(f"systems/{sid}/state.xml"),
                         source_snapshot_sha256=digest("evidence/source_snapshot.json"),
                         preflight_sha256=digest(gate_rel), authorization_sha256=digest("authorization.json"))
    smoke_sha = digest(f"smoke/{sid}/audit.json")
    require(gate["smokes"][sid] == smoke_sha, "smoke receipt drift")
    identity = dict(authorization_identity=authorization, job=row, config_sha256=config_sha,
                    system_sha256=authorization["system_sha256"], topology_sha256=authorization["topology_sha256"],
                    smoke_sha256=smoke_sha)
    prep_rel = f"equilibration/{sid}/replica_{job.replica:02d}"
    ready_rel = f"{prep_rel}/ready_{PROTOCOL}.json"
    ready = read(ready_rel)
    prep_identity = dict(identity, initial_state_sha256=authorization["state_sha256"])
    require(ready.get("passed") is True and ready["identity"] == prep_identity, "preparation identity drift")
    previous = authorization["state_sha256"]
    labels = ("nvt", "npt", "release_01", "release_02", "release_03", "release_04")
    require(len(ready["stage_receipts"]) == len(labels), "preparation stage inventory mismatch")
    for label, expected_sha in zip(labels, ready["stage_receipts"]):
        rel = f"{prep_rel}/{label}/receipt.json"
        require(digest(rel) == expected_sha, "preparation receipt SHA256 mismatch")
        receipt = read(rel)
        require(receipt.get("passed") is True and receipt["identity"] == dict(prep_identity, stage=label, previous_state_sha256=previous),
                "preparation chain mismatch")
        authenticated_outputs(source / prep_rel / label, receipt["outputs"], records)
        previous = receipt["outputs"]["state.xml"]
    require(campaign_relative(ready["start_state"]) == f"{prep_rel}/release_04/state.xml"
            and ready["start_state_sha256"] == previous, "wrong production starting state")
    identity.update(start_state_sha256=previous, preparation_receipt_sha256=digest(ready_rel))
    return identity, records


def read_segment_state(path, end):
    import numpy as np
    rows = [{k.lstrip('#').strip('"'): v for k, v in r.items()} for r in csv_rows(path)]
    steps = list(range(end - SEGMENT_STEPS + 5000, end + 1, 5000))
    require(len(rows) == 10 and [int(r["Step"]) for r in rows] == steps,
            f"segment state step/count mismatch: {path}")
    times = np.asarray([float(r["Time (ps)"]) for r in rows])
    require(np.allclose(times, np.asarray(steps) * .002, rtol=0, atol=1e-4), f"segment state time mismatch: {path}")
    # The new producer does NOT request Density in StateDataReporter.
    for name in ("Potential Energy (kJ/mole)", "Kinetic Energy (kJ/mole)", "Temperature (K)", "Box Volume (nm^3)"):
        require(np.isfinite([float(r[name]) for r in rows]).all(), f"nonfinite state {name}: {path}")
    require(all(float(r["Box Volume (nm^3)"]) > 0 for r in rows), "nonpositive box volume")
    return {"first_step": steps[0], "last_step": steps[-1], "rows": 10,
            "first_time_ps": float(times[0]), "last_time_ps": float(times[-1]),
            "volumes_nm3": [float(r["Box Volume (nm^3)"]) for r in rows]}


def canonical_time(ts):
    # OpenMM creates a fresh DCDReporter per segment, resetting its DCD header
    # clock. CSV and committed step ranges establish the absolute MD time.
    ts.time = (ts.frame + 1) * RAW_FRAME_SPACING_PS
    return ts


def attach_segments(universe, trajectories):
    # continuous=True would discard overlapping header times. Never use it.
    universe.load_new([str(p) for p in trajectories], continuous=False)
    universe.trajectory.add_transformations(canonical_time)


def inspect_segments(job, paths, identity, sequence):
    import numpy as np
    import MDAnalysis as mda
    topology_record = artifact(paths["topology"])
    top = sanitize_topology(paths["topology"], job.system)
    u = mda.Universe(str(top))
    try:
        selection, mapping = topology_selection(u, sequence)
        mapping["atom14_heavy_atoms"] = validate_selected_atoms(u, selection, sequence)
        segments = []
        for index, directory in enumerate(paths["segments"]):
            end = SEGMENT_ENDS[index]
            require(directory.name == f"{end:08d}" and paths["trajectory"][index].resolve() == (directory / "trajectory.dcd").resolve(),
                    "segment directory/order mismatch")
            receipt_path = member(directory, "receipt.json")
            receipt_record = artifact(receipt_path)
            receipt = load_json(receipt_path)
            require(receipt.get("identity") == identity and receipt.get("start_step") == end - SEGMENT_STEPS
                    and receipt.get("end_step") == end, f"segment identity/continuity mismatch: {directory}")
            require(set(receipt["outputs"]) == {"trajectory.dcd", "state.csv", "checkpoint.chk"}, "segment output inventory mismatch")
            records = []
            authenticated_outputs(directory, receipt["outputs"], records)
            state = read_segment_state(directory / "state.csv", end)
            u.load_new(str(directory / "trajectory.dcd"))
            try:
                require(len(u.trajectory) == 10 and np.isclose(u.trajectory.dt, 10, rtol=0, atol=1e-3),
                        f"segment DCD frame/dt mismatch: {directory}")
                require(u.trajectory.n_atoms == len(u.atoms), f"DCD/topology atom count mismatch: {directory}")
                times, fingerprints, boxes = [], [], []
                for ts in u.trajectory:
                    require(np.isfinite(ts.positions).all(), f"nonfinite DCD positions: {directory}")
                    require(ts.dimensions is not None and np.isfinite(ts.dimensions).all()
                            and np.all(ts.dimensions[:3] > 0) and np.all(ts.dimensions[3:] > 0)
                            and np.all(ts.dimensions[3:] < 180) and ts.volume > 0, f"invalid membrane box: {directory}")
                    require(np.isclose(ts.volume / 1000, state["volumes_nm3"][ts.frame], rtol=5e-5, atol=1e-3),
                            f"DCD/state box volume mismatch: {directory}")
                    times.append(ts.time)
                    fingerprints.append(hashlib.sha256(ts.positions.tobytes()).hexdigest())
                    boxes.append(ts.dimensions.tolist())
                local = np.arange(1, 11) * 10
                absolute = local + index * 100
                require(np.allclose(times, local, rtol=0, atol=.01) or np.allclose(times, absolute, rtol=0, atol=.01),
                        f"unexpected DCD header clock: {directory}")
            finally:
                u.trajectory.close()
            verify_artifact(receipt_record)
            segments.append({"start_step": end - SEGMENT_STEPS, "end_step": end,
                             "receipt": receipt_record, "outputs": records, "state": state,
                             "atoms": len(u.atoms), "coordinate_frame_sha256": fingerprints, "boxes_angstrom_degrees": boxes,
                             "dcd_header_first_ps": times[0], "dcd_header_last_ps": times[-1]})
        attach_segments(u, paths["trajectory"])
        try:
            require(len(u.trajectory) == 10 * len(segments), "chain lost/duplicated frames")
            # Prove every chained coordinate and box matches its original file.
            for index in range(len(u.trajectory)):
                ts = u.trajectory[index]
                require(np.isclose(ts.time, (index + 1) * 10), "chain time discontinuity")
                segment, local = segments[index // 10], index % 10
                require(hashlib.sha256(ts.positions.tobytes()).hexdigest() == segment["coordinate_frame_sha256"][local]
                        and np.array_equal(ts.dimensions, segment["boxes_angstrom_degrees"][local]),
                        "chain coordinate/box mismatch (lost, duplicated or reordered frame)")
        finally:
            u.trajectory.close()
        verify_artifact(topology_record)
        return {**job.identity(), "topology": topology_record, "sanitized_topology": artifact(top),
                "trajectory": [next(r for r in seg["outputs"] if Path(r["path"]).name == "trajectory.dcd") for seg in segments],
                "progress": artifact(paths["progress"]),
                "segments": segments, "selection": selection, "receptor": mapping,
                "atoms": len(u.atoms), "chain_coordinates_and_boxes_verified": True,
                "dcd_frames": 10 * len(segments), "dcd_spacing_ps": 10, "coordinates_finite": True}
    finally:
        u.trajectory.close()


def verify_input_record(record):
    seen = set()
    def check(asset):
        key = (asset["path"], asset["sha256"], asset["bytes"])
        if key not in seen:
            verify_artifact(asset)
            seen.add(key)
    for key in ("topology", "sanitized_topology", "progress"):
        check(record[key])
    for asset in [*record["trajectory"], *record.get("identity_assets", [])]:
        check(asset)
    for seg in record["segments"]:
        check(seg["receipt"])
        for asset in seg["outputs"]:
            check(asset)


def input_availability(source=SOURCE_ROOT):
    records = []
    for job in expected_jobs():
        try:
            discover_inputs(job, source)
            records.append({**job.identity(), "files_present": True})
        except (RuntimeError, FileNotFoundError) as exc:
            records.append({**job.identity(), "files_present": False, "reason": str(exc)})
    return {"source_root": str(Path(source).resolve()), "expected_trajectories": 36,
            "file_sets_present": sum(x["files_present"] for x in records), "jobs": records}


def frozen_sequence():
    return "".join(AA[r["resname"]] for r in csv_rows(FROZEN_ROOT / MAPPING_REL))


@exclusive_analysis
def authenticate_inputs(source=SOURCE_ROOT):
    _, anchor = verify_frozen_anchor()
    rows = validate_matrix(source)
    campaign = verify_sums(source)
    sequence = frozen_sequence()
    jobs = []
    for index, job in enumerate(expected_jobs(), 1):
        print(f"authenticate {index}/36 {job.key}", flush=True)
        paths = discover_inputs(job, source)
        identity, assets = production_identity(job, source, rows[job_row(job)["job_id"]])
        before = artifact(paths["progress"])
        progress = load_json(paths["progress"])
        require(progress.get("identity") == identity and progress.get("status") == "COMPLETE"
                and progress.get("completed_steps") == FINAL_STEP and progress.get("completed_ns") == 10,
                f"production not complete/progress identity mismatch: {job.key}")
        record = inspect_segments(job, paths, identity, sequence)
        record["identity_assets"] = assets
        record["production_identity"] = identity
        record["acceptance_reason"] = "100_COMMITTED_SEGMENTS_PROGRESS_COMPLETE"
        verify_artifact(before)
        verify_input_record(record)  # detect changed downloads during scan
        jobs.append(replace(job, topology=Path(record["sanitized_topology"]["path"]),
                            trajectory=paths["trajectory"], sequence=sequence,
                            selection=record["selection"], input_record=record))
    for asset in campaign:
        verify_artifact(asset)
    receipt = {"schema": "pacer.stage4.pbcfix01.input_authentication.v1", "status": "36_OF_36_ACCEPTED",
               "source_root": str(Path(source).resolve()),
               "job_matrix": artifact(Path(source) / "manifests/STAGE4_RERUN_JOB_MATRIX.csv"),
               "campaign_assets": campaign, "jobs": [j.input_record for j in jobs],
               "historical_anchor": anchor["receipt"], "temporal_contract": TEMPORAL_CONTRACT,
               "implementation": code_identity(), "forbidden_operations_performed": FORBIDDEN}
    write_json(RESULT_ROOT / "input_audit/INPUT_AUTHENTICATION_v01.json", receipt)
    return jobs, receipt


def parser(description):
    p = argparse.ArgumentParser(description=description)
    group = p.add_mutually_exclusive_group(required=True)
    for name in ("smoke", "run", "verify", "preflight"):
        group.add_argument(f"--{name}", action="store_true")
    p.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    return p


def preflight(source=SOURCE_ROOT):
    validate_matrix(source)
    assets = verify_sums(source)
    _, anchor = verify_frozen_anchor()
    return {"historical_freeze_verified": True, "historical_anchor": anchor["receipt"],
            "campaign_assets": assets, "temporal_contract": TEMPORAL_CONTRACT,
            "availability": input_availability(source), "full_authentication_performed": False}


