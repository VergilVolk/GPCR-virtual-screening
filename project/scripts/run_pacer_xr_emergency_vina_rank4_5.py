from pathlib import Path
import sys
ROOT = Path(r"C:\\projects\\GPCR-virtual-screening")
SCRIPTS = ROOT / "project" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import run_pacer_xr_vina_rerun_v01 as v
ANCHORS=[
    "PACERGEN01954",
    "PACERGEN01955"
]
full=v.select_scope(False)
roster=full[full["candidate_id"].isin(ANCHORS)].copy().sort_values("rerun_final_rank")
assert len(roster)==2, roster
EMERGENCY_ROOT=ROOT/"project"/"results"/"pacer_xr_emergency_rank4_5_v01"
def emergency_scope(_acceptance=False): return roster.copy()
def emergency_output_root(_acceptance=False): return EMERGENCY_ROOT
v.select_scope=emergency_scope
v.output_root=emergency_output_root
v.run_vina.__globals__["select_scope"]=emergency_scope
v.run_vina.__globals__["output_root"]=emergency_output_root
print("EMERGENCY VINA ANCHORS3")
print(roster[["candidate_id","rerun_final_rank"]].to_string(index=False))
print("Output:",EMERGENCY_ROOT)
print("Starting Vina anchor campaign ...")
v.run_vina(False,workers=12)
