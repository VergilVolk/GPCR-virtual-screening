"""One-command production orchestrator; preflight is read-only, --run is explicit."""
import argparse,json
from pacer_xr_production_v01 import *
from run_pacer_xr_glide_rerun_v01 import run_glide
from run_pacer_xr_vina_rerun_v01 import run_vina
from assemble_pacer_xr_six_channels_v01 import assemble_production,verify_channel_provenance
from run_pacer_xr_rerun_v01 import run_fusion

def preflight():
    roster=load_roster(); anchors=verify_anchors(roster); assets=load_asset_lock()
    verify_environment(assets,'glide'); verify_environment(assets,'vina'); refs,lock=load_references()
    from pacer_xr_glide_scheduler_v01 import science_contract,scheduler_contract
    science_contract(); scheduler=scheduler_contract()
    return dict(status='PASS',source_sha256=SOURCE_SHA,rows=len(roster),anchors=anchors,grids=len(assets['receptors']),glide_scheduler=scheduler,reference_channels={k:len(v) for k,v in refs.items()},reference_commit=lock['reference_commit'])

def stage_ready(root,engine,roster):
    try: verify_channel_provenance(root,engine,roster); return True
    except (ValueError,FileNotFoundError,KeyError,json.JSONDecodeError): return False

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True); mode.add_argument('--run',action='store_true'); mode.add_argument('--acceptance',action='store_true'); mode.add_argument('--preflight',action='store_true')
    ap.add_argument('--glide-workers',type=int,default=3); ap.add_argument('--vina-workers',type=int,default=8)
    args=ap.parse_args(); require(args.glide_workers>0 and args.vina_workers>0,'Workers must be positive')
    report=preflight(); print(json.dumps(report,indent=2),flush=True)
    if args.preflight: return
    acceptance=args.acceptance; roster=select_scope(acceptance); root=output_root(acceptance)
    plan=['verify_rerun_top200','Glide_PDB_and_ensemble','Vina_PDB_and_ensemble','six_channel_assembly','reference_verification','PACER_XR','cascade','full200','natural_top20','final_shortlist20']
    write_json(root/'orchestration_plan.json',dict(scope='acceptance_rank1_3' if acceptance else 'full200',stages=plan,acceptance_stops_after_raw_assembly=acceptance))
    try:
        for engine,runner,workers in [('glide',run_glide,args.glide_workers),('vina',run_vina,args.vina_workers)]:
            if stage_ready(root,engine,roster): print(engine+' complete valid stage RESUMED',flush=True)
            else: runner(acceptance,workers)
        master=JobStore(root/'orchestrator_jobs')
        inputs=dict(source=SOURCE_SHA,codes=code_hashes(),glide=sha(root/'glide/glide_channels.csv'),vina=sha(root/'vina/vina_channels.csv'))
        def assembly_operation(work):
            assemble_production(acceptance)
            return dict(stage='six_channel_assembly',rows=len(roster)),[root/'candidate_six_channel_raw_scores.csv',root/'six_channel_provenance.json']
        master.execute('assemble',inputs,assembly_operation)
        load_references()
        if not acceptance:
            def fusion_operation(work):
                result=run_fusion()
                return result,[OUTPUT/name for name in ['pacer_xr_rerun_full200_ranking.csv','pacer_xr_natural_top20.csv','pacer_xr_final_shortlist20.csv','fusion_audit.json']]
            master.execute('fusion',dict(raw=sha(root/'candidate_six_channel_raw_scores.csv'),refs=sha(P/'config/pacer_xr_references_v01.json'),codes=code_hashes()),fusion_operation)
        write_json(root/'orchestration_status.json',dict(status='PASS',scope='acceptance_rank1_3' if acceptance else 'full200',fusion_executed=not acceptance,final_selection_executed=not acceptance,resume_counts=master.counts))
    except Exception as exc:
        write_json(root/'orchestration_status.json',dict(status='FAILED',error=str(exc),scope='acceptance_rank1_3' if acceptance else 'full200')); raise

if __name__=='__main__': main()
