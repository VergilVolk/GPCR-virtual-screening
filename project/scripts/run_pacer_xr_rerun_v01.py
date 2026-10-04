"""Full200-only PACER-XR fusion, frozen cascade, and experimental shortlist."""
import argparse,json
import pandas as pd
from pacer_xr_production_v01 import *
from assemble_pacer_xr_six_channels_v01 import verify_channel_provenance

def run_fusion():
    roster=load_roster(); anchors=verify_anchors(roster); load_asset_lock()
    references,reference_lock=load_references()
    verify_channel_provenance(OUTPUT,'glide',roster); verify_channel_provenance(OUTPUT,'vina',roster)
    raw=OUTPUT/'candidate_six_channel_raw_scores.csv'; scores=pd.read_csv(raw)
    prov=json.loads((OUTPUT/'six_channel_provenance.json').read_text())
    require(prov['complete'] and prov['scope']=='full200' and prov['n_candidates']==200,'Fusion requires full200 production provenance')
    require(prov['raw_sha256']==sha(raw) and prov['source_sha256']==SOURCE_SHA and prov['codes']==code_hashes(),'Six-channel provenance mismatch')
    require(prov['glide_sha256']==sha(OUTPUT/'glide/glide_channels.csv') and prov['vina_sha256']==sha(OUTPUT/'vina/vina_channels.csv'),'Assembler is stale relative to engine outputs')
    ranked,threshold=fuse(scores,roster,references)
    require(len(ranked)==200,'Full ranking must contain exactly 200 candidates')
    natural=natural_top20(ranked); final=shortlist20(ranked)
    # All checks precede any final-ranking/selection writes.
    write_csv(OUTPUT/'pacer_xr_rerun_full200_ranking.csv',ranked)
    write_csv(OUTPUT/'pacer_xr_natural_top20.csv',natural)
    write_csv(OUTPUT/'pacer_xr_final_shortlist20.csv',final)
    audit=dict(source_sha256=SOURCE_SHA,n_candidates=200,reference_lock_sha256=sha(P/'config/pacer_xr_references_v01.json'),reference_channels=reference_lock['channels'],Glide_BEmin_official_top1pct_raw_threshold=threshold,cascade_head_candidates=int((ranked.cascade_stage==1).sum()),anchors=anchors,all_six_complete=True,normalization='Existing empirical midrank; lower energies yield higher ranks',cascade='Stable source-rank order on ties; stage1 Glide_BEmin rank descending; remaining PACER_XR descending',anchor_rule='Applied only to final shortlist; full200 scores and ranking unchanged',claim_boundary='Ranking hypothesis; not a functional PAM probability')
    write_json(OUTPUT/'fusion_audit.json',audit)
    return audit

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--run',action='store_true',required=True)
    ap.parse_args(); run_fusion()

if __name__=='__main__': main()
