from pathlib import Path
import shutil
import sys

ROOT = Path(r"C:\projects\GPCR-virtual-screening")
SCRIPTS = ROOT / "project" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import run_pacer_xr_glide_rerun_v01 as g

ANCHORS = [
    "PACERGEN01755",  # PACER0010
    "PACERGEN00094",  # PACER0073
    "PACERGEN00123",  # PACER0027
]

# Load and verify the official frozen top200 first.
full = g.select_scope(False)
roster = full[full["candidate_id"].isin(ANCHORS)].copy()
roster = roster.sort_values("rerun_final_rank")

assert len(roster) == 3, roster
assert set(roster["candidate_id"]) == set(ANCHORS)

# Independent emergency output: never overwrite full200 summaries/channels.
EMERGENCY_ROOT = (
    ROOT / "project" / "results" / "pacer_xr_emergency_anchors3_v01"
)

# Reuse the already-completed full200 Glide preparation jobs.
src_jobs = (
    ROOT / "project" / "results" / "pacer_xr_rerun_v01"
    / "glide" / "jobs"
)
dst_jobs = EMERGENCY_ROOT / "glide" / "jobs"
dst_jobs.mkdir(parents=True, exist_ok=True)

for cid in ANCHORS:
    src = src_jobs / f"prepare_{cid}"
    dst = dst_jobs / f"prepare_{cid}"
    assert src.exists(), f"Missing cached Glide preparation: {src}"

    if not dst.exists():
        shutil.copytree(src, dst)

def emergency_scope(_acceptance=False):
    return roster.copy()

def emergency_output_root(_acceptance=False):
    return EMERGENCY_ROOT

# Patch this module.
g.select_scope = emergency_scope
g.output_root = emergency_output_root

# Also patch globals used by the optimized/legacy implementations.
for name in ("run_optimized", "run_glide_legacy"):
    fn = getattr(g, name, None)
    if fn is not None:
        fn.__globals__["select_scope"] = emergency_scope
        fn.__globals__["output_root"] = emergency_output_root

print("EMERGENCY GLIDE ANCHORS3")
print(roster[["candidate_id", "rerun_final_rank"]].to_string(index=False))
print("Output:", EMERGENCY_ROOT)
print("Starting 3 x 11 Glide evaluations ...")

g.run_glide_legacy(False, workers=2)