"""Read-only PACER-FKG v02 Phase-0 inventory/compatibility audit.

This program never imports or runs ViSNet.  Its only writes are guarded v02
JSON/CSV outputs; all MD and historical inputs are opened read-only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.coordinates.DCD import DCDReader

from path_guard import guard_write_path


REPO = Path(__file__).resolve().parents[2]
V02 = REPO / "project/pacer_fkg_v02"
PROV = REPO / "project/results/pacer_fkg_v02_longmd_v01/provenance"
BACKUP = Path(r"C:\projects\PACER_DC_MD_backup")
SHORT = Path(r"C:\projects\PACER_DC_MD_delivery_20260927")
BASELINE = "0a4b4f238c73add75a09a0875f698250a39bdf33"
SYSTEMS = ("apo", "probe_only", "compound110__candidate_no_probe", "compound110__candidate_probe")
CONTEXT = {
    "apo": "apo", "probe_only": "probe_only",
    "compound110__candidate_no_probe": "candidate_no_probe",
    "compound110__candidate_probe": "candidate_probe",
}
SEEDS = {1: 27101, 2: 38201, 3: 49301}
RESNAME_MAP = {"HSE": "HIS", "HSD": "HIS", "HSP": "HIS"}
AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}
ATOM14 = {
    "A":"N CA C O CB", "R":"N CA C O CB CG CD NE CZ NH1 NH2", "N":"N CA C O CB CG OD1 ND2",
    "D":"N CA C O CB CG OD1 OD2", "C":"N CA C O CB SG", "Q":"N CA C O CB CG CD OE1 NE2",
    "E":"N CA C O CB CG CD OE1 OE2", "G":"N CA C O", "H":"N CA C O CB CG ND1 CD2 CE1 NE2",
    "I":"N CA C O CB CG1 CG2 CD1", "L":"N CA C O CB CG CD1 CD2", "K":"N CA C O CB CG CD CE NZ",
    "M":"N CA C O CB CG SD CE", "F":"N CA C O CB CG CD1 CD2 CE1 CE2 CZ", "P":"N CA C O CB CG CD",
    "S":"N CA C O CB OG", "T":"N CA C O CB OG1 CG2", "W":"N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "Y":"N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "V":"N CA C O CB CG1 CG2",
}
ATOM14 = {k: v.split() for k, v in ATOM14.items()}


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path, block: int = 16 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, obj) -> None:
    safe = guard_write_path(path, REPO)
    safe.parent.mkdir(parents=True, exist_ok=True)
    safe.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=REPO, text=True, encoding="utf-8").strip()


def baseline_manifest() -> dict:
    # The commit determines membership; SHA256 protects the actual local
    # checkout bytes, including the platform's configured line endings.
    raw = subprocess.check_output(
        ["git", "ls-tree", "-r", "-z", "--format=%(objectname)%x09%(objectsize)%x09%(path)", BASELINE], cwd=REPO
    ).decode("utf-8")
    rows = []
    for record in raw.rstrip("\0").split("\0"):
        oid, size, path = record.split("\t", 2)
        content = (REPO / path).read_bytes()
        rows.append({"path": path.replace("\\", "/"), "sha256": hashlib.sha256(content).hexdigest(), "size": len(content)})
    return {
        "schema": "pacer_fkg_v02.baseline_tracked_tree.v1", "baseline_commit": BASELINE,
        "branch": git("branch", "--show-current"), "tag": "encoder-candidate-ah-frozen-20260929",
        "manifest_creation_time": now(), "hash_source": "local worktree bytes for paths tracked at baseline commit",
        "total_file_count": len(rows), "files": rows,
    }


def component_inventory() -> dict:
    tracked = set(git("ls-tree", "-r", "--name-only", BASELINE).splitlines())
    entries: list[dict] = []
    def add(path: str, classification: str, participation: list[str], note: str = "") -> None:
        p = REPO / path
        entries.append({"path": path, "classification": classification, "participation": participation,
                        "present_locally": p.exists(), "tracked_at_baseline": path in tracked, "note": note})
    source = {
      "project/pacer_dc_training/pacer_factorial_kernel_graph.py": ["kernel", "contrasts", "graph diffusion"],
      "project/pacer_dc_training/analyze_pacer_fkg_g1.py": ["G1", "regions", "contrast diagnostics"],
      "project/pacer_dc_training/audit_pacer_fkg_g2a_atom14_pbc.py": ["G2A", "atom14 mapping", "PBC"],
      "project/pacer_dc_training/audit_pacer_fkg_g2a_pbc.py": ["G2A PBC", "PBC"],
      "project/pacer_dc_training/analyze_pacer_fkg_g2b.py": ["G2B", "normalization", "bandwidth", "RFF", "contrasts"],
      "project/pacer_dc_training/audit_pacer_fkg_g2b_robustness.py": ["G2B robustness", "kernel", "RFF"],
      "project/pacer_dc_training/analyze_pacer_fkg_g2c.py": ["G2C", "graph", "RFF"],
      "project/pacer_dc_training/run_pacer_fkg_r2r3_audit.py": ["R2/R3 FKG audit", "kernel", "regions", "contrasts"],
      "project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py": ["authenticated short-data preprocessing", "atom14", "selection"],
      "project/pacer_dc_training/extract_geom2vec_atom14.py": ["Geom2Vec input contract", "atom14 heavy atoms"],
      "project/pacer_dc_training/build_m4_multistructure_graph.py": ["graph definitions", "regions"],
      "project/scripts/run_pacer_dc_production_md.py": ["production trajectory provenance", "PBC writer", "seeds"],
      "project/pacer_dc_training/FOUR_CONTEXT_WINDOW_CONTRACT.md": ["short-data contract", "contrasts"],
      "project/docs/PACER_FKG_IMPLEMENTATION_AUDIT_20260927.md": ["implementation audit"],
      "project/docs/PACER_FACTORIAL_KERNEL_GRAPH_METHOD.md": ["method definition"],
    }
    for p, roles in source.items(): add(p, "tracked source/code", roles)
    result_files = []
    roots = (
      "project/results/pacer_dc_fkg_G1_diagnostics_v01", "project/results/pacer_dc_fkg_G2A_v01",
      "project/results/pacer_dc_fkg_G2A_PBC_v01", "project/results/pacer_dc_fkg_G2B_v01",
      "project/results/pacer_dc_fkg_G2B_robustness_v01", "project/results/pacer_dc_fkg_G2C_v01",
      "project/results/pacer_dc_fkg_R2R3_v01", "project/results/pacer_dc_geom2vec_pilot_v01",
      "project/results/pacer_dc_fkg_G2A_v02", "project/results/pacer_dc_fkg_U2_erratum_v02",
      "project/results/pacer_dc_four_context_v01", "project/results/pacer_dc_geom2vec_R2R3_full_v01",
    )
    for path in sorted(tracked):
        if any(path == r or path.startswith(r + "/") for r in roots):
            result_files.append(path)
            roles = ["historical v01 result/calibration"]
            if "G2B_DIRECTION_AUDIT" in path: roles += ["normalization state", "bandwidth", "RFF state"]
            if "M4_MULTISTRUCTURE_GRAPH" in path: roles += ["graph definitions", "region definitions"]
            add(path, "tracked result/calibration artifact", roles)
    generated = (
      "project/results/pacer_dc_four_context_v01/compound110",
      "project/results/pacer_dc_geom2vec_R2R3_full_v01",
      "project/results/pacer_dc_G2A_source_v01",
      "project/results/pacer_dc_production_v01",
      "project/results/pacer_dc_membrane_reference_v01",
      "project/results/pacer_dc_restraint_release_v01",
      "project/results/pacer_dc_common_protein_v01",
    )
    for p in generated: add(p, "generated/ignored artifact", ["historical v01 input/state"], "Directory-level artifact; may be absent in this worktree.")
    add(str(SHORT), "external dependency", ["authenticated short-data mapping", "topologies", "atom14", "raw windows"], "Read-only local delivery with SHA256SUMS.txt and manifest.json.")
    add(r"C:\projects\geom2vec-source", "external dependency", ["Geom2Vec implementation dependency"])
    add("project/tools/geom2vec_checkpoints/visnet-base-1.pt", "external dependency", ["frozen ViSNet checkpoint"], "Path referenced by historical audit; file may be excluded/absent.")
    counts = Counter(x["classification"] for x in entries)
    return {"schema":"pacer_fkg_v02.v01_component_inventory.v1", "created_at":now(),
            "baseline_commit":BASELINE, "entries":entries, "counts":dict(counts),
            "scope_note":"Exact files are listed for tracked implementation/results; bulky generated state is identified by exact directory root."}


def parse_server_sums() -> dict[str, str]:
    out = {}
    for line in (BACKUP / "SERVER_SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        digest, rel = line.split(None, 1)
        out[rel.strip().removeprefix("./").replace("\\", "/")] = digest.lower()
    return out


def role(path: Path) -> str:
    n = path.name
    return {"trajectory.dcd":"trajectory", "progress.json":"completion_metadata", "state.csv":"thermodynamic_qc",
            "checkpoint.chk":"restart_checkpoint", "latest_state.xml":"final_state"}.get(n,
            "transfer_checksum_manifest" if n.startswith("SERVER_SHA256SUMS") else
            "backup_validation" if n.startswith("BACKUP_VALIDATION") else
            "transfer_log" if path.parent.name == "_rsync_logs" else "other")


def file_inventory() -> tuple[dict, bool]:
    expected_manifest = (BACKUP / "SERVER_SHA256SUMS.txt.sha256").read_text().split()[0].lower()
    actual_manifest = sha256(BACKUP / "SERVER_SHA256SUMS.txt")
    server = parse_server_sums()
    files, mismatch = [], False
    for p in sorted((x for x in BACKUP.rglob("*") if x.is_file()), key=lambda x: x.as_posix()):
        rel = p.relative_to(BACKUP).as_posix(); local = sha256(p); remote = server.get(rel)
        status = "MATCH" if remote == local else "NOT_RECORDED" if remote is None else "MISMATCH"
        mismatch |= status == "MISMATCH"
        files.append({"relative_path":rel, "role":role(p), "size":p.stat().st_size,
                      "local_sha256":local, "server_recorded_sha256":remote, "hash_match_status":status})
        print(f"HASH {status:12} {rel}", flush=True)
    required = [f"{s}/replica_{r:02d}/{n}" for s in SYSTEMS for r in (1,2,3)
                for n in ("trajectory.dcd","state.csv","checkpoint.chk","latest_state.xml","progress.json")]
    return ({"schema":"pacer_fkg_v02.long_md_file_inventory.v1", "created_at":now(),
             "backup_root":str(BACKUP), "server_manifest_companion_expected_sha256":expected_manifest,
             "server_manifest_local_sha256":actual_manifest, "server_manifest_companion_match":expected_manifest == actual_manifest,
             "server_record_count":len(server), "all_server_recorded_files_match":not mismatch,
             "expected_60_production_files_present":all((BACKUP/x).is_file() for x in required),
             "file_count":len(files), "files":files}, mismatch or expected_manifest != actual_manifest)


def short_reference() -> tuple[dict, dict[str, Path]]:
    manifest = {x["path"]: x for x in json.loads((SHORT/"manifest.json").read_text(encoding="utf-8-sig"))}
    topologies, evidence = {}, {}
    for system in SYSTEMS:
        c = CONTEXT[system]
        rel = f"compound110/replica_02/window_000/topology/{c}_w000_minimized.pdb"
        seqrel = f"compound110/replica_02/window_000/atom14/{c}_w000.csv"
        rawrel = f"compound110/replica_02/window_000/raw/{c}_w000.npz"
        paths = {"topology":SHORT/rel, "sequence":SHORT/seqrel, "raw":SHORT/rawrel}
        checks = {k:{"path":str(v), "sha256":sha256(v), "manifest_sha256":manifest[v.relative_to(SHORT).as_posix()]["sha256"]}
                  for k,v in paths.items()}
        for item in checks.values(): item["match"] = item["sha256"] == item["manifest_sha256"]
        evidence[system] = checks; topologies[system] = paths["topology"]
    return evidence, topologies


def min_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box[:3] * np.round(delta / box[:3])


def inspect_trajectories(topologies: dict[str, Path]) -> tuple[list[dict], list[dict]]:
    metadata, compat = [], []
    for system in SYSTEMS:
      for replica in (1,2,3):
        traj = BACKUP/system/f"replica_{replica:02d}"/"trajectory.dcd"; topo = topologies[system]
        with DCDReader(str(traj), convert_units=True) as reader:
            nframes, natoms, dt = len(reader), reader.n_atoms, float(reader.dt)
            sample_ids = sorted(set((0, nframes//4, nframes//2, (3*nframes)//4, nframes-1)))
            samples=[]; finite=True; boxes=[]
            for i in sample_ids:
                ts=reader[i]; finite &= bool(np.isfinite(ts.positions).all()); boxes.append(ts.dimensions.tolist() if ts.dimensions is not None else None)
                samples.append({"frame":i,"time_ps":float(ts.time),"finite":bool(np.isfinite(ts.positions).all()),
                                "box":ts.dimensions.tolist() if ts.dimensions is not None else None})
            first_time=float(reader[0].time); last_time=float(reader[-1].time)
        u=Universe(str(topo), str(traj)); receptor=u.select_atoms("chainID E and protein")
        residues=list(receptor.residues); names=np.asarray(receptor.names); ca=receptor.select_atoms("name CA")
        seq=[]; unsupported=[]; missing=[]; unexpected=[]; masks=[]; atom14_indices=[]
        for j,res in enumerate(residues):
            std=RESNAME_MAP.get(res.resname,res.resname); aa=AA3.get(std)
            if aa is None: unsupported.append({"index":j,"resid":int(res.resid),"resname":res.resname}); seq.append("X"); continue
            seq.append(aa); local={("CD1" if std=="ILE" and x.name=="CD" else x.name):x.index for x in res.atoms}
            expected=ATOM14[aa]; miss=[x for x in expected if x not in local]
            if miss: missing.append({"index":j,"resid":int(res.resid),"resname":res.resname,"atoms":miss})
            masks.append([x in local for x in expected] + [False]*(14-len(expected)))
            atom14_indices.extend(local[x] for x in expected if x in local)
            known=set(expected)|{"H","HN","HA","HA2","HA3","HB1","HB2","HB3","HG1","HG2","HG3","HD1","HD2","HD3","HE1","HE2","HE3","HZ","HZ1","HZ2","HZ3","HH","HH11","HH12","HH21","HH22","HT1","HT2","HT3","OXT"}
            extra=sorted({x.name for x in res.atoms if not x.name.startswith("H") and x.name not in known})
            if extra: unexpected.append({"index":j,"resid":int(res.resid),"resname":res.resname,"atoms":extra})
        sequence="".join(seq)
        short_csv=SHORT/f"compound110/replica_02/window_000/atom14/{CONTEXT[system]}_w000.csv"
        short_sequence=short_csv.read_text(encoding="utf-8").splitlines()[1].split(",",1)[1].strip()
        # Full-frame CA PBC audit. No transformations are attached to Universe.
        valid_box=0; invalid_box=0; split_pairs=0; severe_pairs=0; temporal_wrap=0; max_direct=0.; prev=None
        frame_finite=True
        for ts in u.trajectory:
            pos=ca.positions.astype(np.float64, copy=True); frame_finite &= bool(np.isfinite(pos).all())
            box=ts.dimensions
            if box is None or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0): invalid_box+=1; continue
            valid_box+=1; direct=np.diff(pos,axis=0); dd=np.linalg.norm(direct,axis=1); mi=np.linalg.norm(min_image(direct,box),axis=1)
            split_pairs += int(np.sum((dd>10.0)&(mi<5.0))); severe_pairs += int(np.sum((dd>10.0)&(mi>=5.0))); max_direct=max(max_direct,float(dd.max()))
            if prev is not None:
                jump=pos-prev; jd=np.linalg.norm(jump,axis=1); jm=np.linalg.norm(min_image(jump,box),axis=1)
                temporal_wrap += int(np.sum((jd>box[:3].min()/2)&(jm<5.0)))
            prev=pos
        state=list(csv.DictReader((BACKUP/system/f"replica_{replica:02d}"/"state.csv").open(newline="",encoding="utf-8")))
        state_times=[float(x['Time (ps)']) for x in state]
        progress=json.loads((BACKUP/system/f"replica_{replica:02d}"/"progress.json").read_text())
        rec={"system":system,"replica":replica,"topology_selected":str(topo),"trajectory_selected":str(traj),
             "atom_count":natoms,"topology_atom_count":len(u.atoms),"residue_count":len(u.residues),
             "receptor_atom_count":len(receptor),"receptor_residue_count":len(residues),"stored_frames":nframes,
             "first_trajectory_time_ps":first_time,"last_trajectory_time_ps":last_time,"dt_ps":dt,
             "state_csv_rows":len(state),"state_csv_first_time_ps":state_times[0],"state_csv_last_time_ps":state_times[-1],
             "state_csv_minus_dcd_first_time_ps":state_times[0]-first_time,
             "state_csv_minus_dcd_last_time_ps":state_times[-1]-last_time,
             "box_pbc_available":valid_box==nframes,"valid_box_frames":valid_box,"invalid_box_frames":invalid_box,
             "expected_50_ns_consistency":bool(nframes>1 and math.isclose(last_time,50000.,abs_tol=.1) and progress.get('completed_ns')==50.0),
             "finite_coordinate_sample":finite,"sample_frames":samples,"progress":progress,
             "pbc_full_scan":{"frames":nframes,"adjacent_ca_split_pairs":split_pairs,"adjacent_ca_non_pbc_severe_pairs":severe_pairs,
                              "temporal_ca_wrapping_events":temporal_wrap,"max_adjacent_ca_direct_A":max_direct}}
        metadata.append(rec)
        compat.append({"system":system,"replica":replica,"sequence":sequence,"short_reference_sequence":short_sequence,
          "sequence_matches_authenticated_short_data":sequence==short_sequence,"expected_270_residues":len(sequence)==270,
          "residue_order_consistent":sequence==short_sequence,"unsupported_or_unexpected_residue_names":unsupported,
          "atom14_slot_count":int(sum(len(ATOM14.get(x,[])) for x in sequence)),"atom14_mask_sha256":hashlib.sha256(np.asarray(masks,dtype=np.uint8).tobytes()).hexdigest(),
          "missing_expected_atoms":missing,"unexpected_heavy_atoms":unexpected,
          "backbone_N_CA_C_O_complete":not any(set(x["atoms"]) & {"N","CA","C","O"} for x in missing),
          "glycine_count":sequence.count("G"),"glycine_has_no_cb_by_contract":True,
          "terminal_atom_policy":"OXT and terminal hydrogens are ignored because they are not atom14 slots",
          "finite_receptor_ca_coordinates_full_trajectory":frame_finite,
          "finite_all_atom_coordinates_deterministic_sample":finite,
          "pbc_route_compatible":split_pairs==0 and severe_pairs==0 and temporal_wrap==0,
          "compatible":bool(len(u.atoms)==natoms and sequence==short_sequence and len(sequence)==270 and not unsupported and not missing and finite and frame_finite and split_pairs==0 and severe_pairs==0 and temporal_wrap==0)})
        print(f"SCAN {system} R{replica}: {nframes} frames, {natoms} atoms, compatible={compat[-1]['compatible']}",flush=True)
    return metadata,compat


def completion_qc(metadata: list[dict]) -> dict:
    rows=[]
    for x in metadata:
        p=x["progress"]; seed=p.get("seed"); expected=SEEDS[x["replica"]]
        rows.append({"system":x["system"],"replica":x["replica"],"status":p.get("status"),"target_ns":p.get("target_ns"),
          "completed_ns":p.get("completed_ns"),"seed_recorded_locally":seed,"expected_paired_seed":expected,
          "seed_locally_verified":seed==expected,"state_rows":x["state_csv_rows"],"trajectory_frames":x["stored_frames"],
          "completion_consistent":p.get("status")=="complete" and p.get("completed_ns")==50.0 and x["expected_50_ns_consistency"]})
    return {"schema":"pacer_fkg_v02.long_md_completion_qc.v1","created_at":now(),"expected_design":"4 systems x 3 replicas x 50 ns = 600 ns",
            "paired_seed_contract":SEEDS,"all_12_complete":len(rows)==12 and all(x["completion_consistent"] for x in rows),
            "all_seeds_locally_verified":all(x["seed_locally_verified"] for x in rows),"trajectories":rows}


def preprocessing_audit(short_evidence: dict, metadata: list[dict]) -> dict:
    return {"schema":"pacer_fkg_v02.long_md_encoder_preprocessing_audit.v1","created_at":now(),
      "established_route": ["production trajectory.dcd", "system minimized.pdb topology (CONECT-stripped copy)",
        "MDAnalysis selection: chainID E and protein", "extract selected coordinates without transform",
        "pacer_dc_atom14_v0.1 mapping", "flatten expected atom14 heavy-atom slots", "ViSNet (NOT RUN IN PHASE 0)"],
      "route_source":"project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py",
      "evidence_for_authenticated_route":["external delivery trajectory_windows.json files identify preprocessing_version pacer_dc_atom14_v0.1",
        "delivery raw/atom14/topology hashes are recorded by its manifest.json and SHA256SUMS.txt",
        "historical G2A PBC audit reports raw-coordinate identity and no short-trajectory receptor PBC splits"],
      "molecules_made_whole":False,"receptor_centering":False,"alignment":False,"receptor_selection":"chainID E and protein",
      "reference_structure":"per-system minimized.pdb from pacer_dc_membrane_reference_v01; read-only sanitized copies from PACER_DC_MD_delivery_20260927 used here",
      "missing_atom_behavior":"Expected atom14 slots absent from topology remain zero in historical construction; frozen Geom2Vec flattening nevertheless enumerates expected slots. Phase 0 rejects any missing expected heavy atom rather than silently repairing it.",
      "atom14_ordering":"mdgen restype_name_to_atom14_names; ILE CHARMM CD renamed CD1; HSE/HSD/HSP mapped to HIS",
      "residue_ordering":"MDAnalysis chainID E protein topology order; authenticated reference sequence length 270",
      "periodic_box_treatment":"DCDReporter(enforcePeriodicBox=False); box vectors retained in DCD; extractor applies no unwrap/wrap/whole transform and selects stored receptor coordinates directly",
      "long_md_full_scan_pbc_summary":{"all_frames_have_box":all(x['box_pbc_available'] for x in metadata),
        "adjacent_ca_split_pairs":sum(x['pbc_full_scan']['adjacent_ca_split_pairs'] for x in metadata),
        "non_pbc_severe_pairs":sum(x['pbc_full_scan']['adjacent_ca_non_pbc_severe_pairs'] for x in metadata),
        "temporal_ca_wrapping_events":sum(x['pbc_full_scan']['temporal_ca_wrapping_events'] for x in metadata)},
      "short_reference_authentication":short_evidence,
      "route_usable_without_change":all(x['pbc_full_scan']['adjacent_ca_split_pairs']==0 and x['pbc_full_scan']['adjacent_ca_non_pbc_severe_pairs']==0 and x['pbc_full_scan']['temporal_ca_wrapping_events']==0 for x in metadata),
      "no_new_pbc_procedure_invented":True}


def main() -> None:
    write_json(V02/"BASELINE_TRACKED_TREE_SHA256.json", baseline_manifest())
    write_json(V02/"V01_COMPONENT_INVENTORY.json", component_inventory())
    prior_inventory = PROV/"LONG_MD_FILE_INVENTORY.json"
    if prior_inventory.is_file():
        inventory = json.loads(prior_inventory.read_text(encoding="utf-8"))
        mismatch = (not inventory.get("server_manifest_companion_match", False) or
                    not inventory.get("all_server_recorded_files_match", False))
    else:
        inventory, mismatch=file_inventory(); write_json(prior_inventory,inventory)
    if mismatch: raise SystemExit("STOP: server/local checksum mismatch or manifest authentication failure")
    short_evidence,topologies=short_reference()
    metadata,compat=inspect_trajectories(topologies)
    write_json(PROV/"TRAJECTORY_METADATA_AUDIT.json", {"schema":"pacer_fkg_v02.trajectory_metadata_audit.v1","created_at":now(),"trajectories":metadata})
    write_json(PROV/"LONG_MD_COMPLETION_QC_AUDIT.json", completion_qc(metadata))
    write_json(PROV/"LONG_MD_ENCODER_PREPROCESSING_AUDIT.json", preprocessing_audit(short_evidence,metadata))
    write_json(PROV/"LONG_MD_ENCODER_COMPATIBILITY.json", {"schema":"pacer_fkg_v02.long_md_encoder_compatibility.v1","created_at":now(),
      "frozen_input_contract":"270 ordered receptor residues; deterministic atom14 expected-heavy-atom slots; finite coordinates; established no-transform PBC route",
      "all_12_compatible":len(compat)==12 and all(x['compatible'] for x in compat),"trajectories":compat})


if __name__ == "__main__":
    main()
