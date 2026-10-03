"""Reuse Stage3D selection implementation with corrected wiring and identity gates."""
import ast
import contextlib
import datetime
import hashlib
import json
import runpy
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from identity_contract import validate_inputs, graph, state_parent, sha

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
OLD = ROOT/'project/results/pacer_stage3_v01'
DOCK = ROOT/'project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01'
ROUTER = ROOT/'project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv'

def rec(p):
    return dict(path=str(p.relative_to(ROOT)),sha256=sha(p))

def write(name,obj):
    p=OUT/name
    assert not p.exists(),p
    p.write_text(json.dumps(obj,indent=2)+'\n')

def pose_smiles(p):
    return [line[len('REMARK SMILES '):].strip() for line in p.read_text().splitlines() if line.startswith('REMARK SMILES ') and not line.startswith('REMARK SMILES IDX ')]

def main():
    completed_execution = (OUT/'md_shortlist_pose_manifest.csv').exists()
    identity=validate_inputs(ROOT,DOCK,ROUTER)
    preserved={str(p):sha(p) for p in OLD.rglob('*') if p.is_file()}
    preserved.update({str(ROUTER):sha(ROUTER),str(DOCK/'ledger.jsonl'):sha(DOCK/'ledger.jsonl')})
    originals=[OLD/'stage3d_c_gate.py',OLD/'stage3d_def_merge_diversity_pose.py']
    policy=json.loads((OLD/'STAGE3D_FREEZE_MANIFEST.json').read_text())
    assert policy['selection']['rule']=='structural gate PASS -> frozen M4-safe final_rank -> one representative per Murcko scaffold'
    assert policy['selection']['recommended_md_set_size']==3
    assert policy['structural_gate']['thresholds']==dict(cluster_coverage_fraction=0.8,median_native_pocket_coverage=0.5,median_n_contacted_residues=8,max_centroid_outlier_distance_A=20.0)
    amended=[]
    for source in originals:
        text=source.read_text()
        text=text.replace('OUT  = ROOT / "project/results/pacer_stage3_v01"','OUT  = ROOT / "project/results/pacer_stage3d_corrected_20261002_v01"')
        marker='DOCK = ROOT / "project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01"'
        text=text.replace(marker,marker+'\nfrom identity_contract import validate_inputs, graph\nROUTER_PATH = ROOT / "project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv"\nvalidate_inputs(ROOT, DOCK, ROUTER_PATH)')
        if source.name=='stage3d_def_merge_diversity_pose.py':
            text=text.replace('pd.read_csv(OUT / "pacer200_m4_safe_routed_ranking.csv")','pd.read_csv(ROUTER_PATH)')
            text=text.replace('m = router.merge(gate.drop(columns=["canonical_smiles"]), on="candidate_id", how="left", validate="1:1")','router["stable_graph_identifier"] = router.canonical_smiles.map(graph)\ngate["stable_graph_identifier"] = gate.canonical_smiles.map(graph)\nm = router.merge(gate.drop(columns=["canonical_smiles"]), on=["candidate_id", "stable_graph_identifier"], how="left", validate="1:1")')
            # Repair ID/index alignment of the reporting field; selection still uses final_rank.
            text=text.replace('m["stage3d_priority"] = pd.Series({cid: i for i, cid in enumerate(prio.pair_id, 1)})','m["stage3d_priority"] = m.candidate_id.map({cid: i for i, cid in enumerate(prio.pair_id, 1)})')
        ast.parse(text)
        target=OUT/source.name
        if target.exists():
            assert target.read_text()==text
        else:
            target.write_text(text)
        amended.append(target)
    if not completed_execution:
        with (OUT/'stage3d_stdout.log').open('x') as log,contextlib.redirect_stdout(log):
            for script in amended:runpy.run_path(str(script),run_name='__main__')
    merged=pd.read_csv(OUT/'pacer200_stage3d_merged.csv')
    div=pd.read_csv(OUT/'pacer_stage3d_diverse_shortlist.csv')
    man=pd.read_csv(OUT/'md_shortlist_pose_manifest.csv')
    routed=pd.read_csv(ROUTER)
    assert len(merged)==200 and merged.candidate_id.nunique()==200
    check=merged.set_index('candidate_id'); original=routed.set_index('candidate_id')
    assert check.final_rank.equals(original.loc[check.index,'final_rank'])
    assert np.allclose(check.pacer_binding_score,original.loc[check.index,'pacer_binding_score'])
    assert check.canonical_smiles.map(graph).equals(original.loc[check.index,'canonical_smiles'].map(graph))
    # Recomputed unchanged docking features must reproduce the historical gate table.
    old_gate=pd.read_csv(OLD/'pacer200_structural_gate.csv').sort_values('candidate_id').reset_index(drop=True)
    new_gate=pd.read_csv(OUT/'pacer200_structural_gate.csv').sort_values('candidate_id').reset_index(drop=True)
    pd.testing.assert_frame_equal(old_gate,new_gate,check_exact=False,rtol=1e-12,atol=1e-12)
    selected=div[div.md_set_member].copy()
    assert len(selected)==len(man)==3
    assert selected.murcko_scaffold.nunique()==3
    expected=[];seen=set()
    for row in merged[merged.structural_gate_pass].sort_values(['final_rank','pair_id']).itertuples():
        sc=MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(row.canonical_smiles))
        if sc not in seen:seen.add(sc);expected.append(row.candidate_id)
    assert selected.candidate_id.tolist()==expected[:3]
    ledger=pd.read_json(DOCK/'ledger.jsonl',lines=True)
    frozen=[]
    for row in man.itertuples():
        group=ledger[ledger.molecule_id==row.candidate_id].copy()
        xyz=group[['pose_centroid_x','pose_centroid_y','pose_centroid_z']].to_numpy(float)
        group['d']=np.linalg.norm(xyz-np.median(xyz,axis=0),axis=1)
        pick=group.sort_values(['d','native_pocket_coverage','vina_affinity','cluster'],ascending=[True,False,True,True],kind='mergesort').iloc[0]
        assert row.source_cluster==pick.cluster and row.source_state==pick.state_id
        src=ROOT/row.pose_source_file;dst=ROOT/row.pose_file
        assert sha(src)==sha(dst)==row.pose_file_sha256
        states=pose_smiles(src);assert states and all(graph(s)==graph(row.source_state_smiles) for s in states)
        assert state_parent(row.source_state_smiles)==state_parent(row.canonical_smiles)
        frozen.append(dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,stable_graph_identifier=graph(row.canonical_smiles),source_state=int(row.source_state),source_state_smiles=row.source_state_smiles,source_cluster=int(row.source_cluster),source_pose=rec(src),frozen_pose=rec(dst),router_final_rank=int(row.router_final_rank),diverse_rank=int(row.diverse_rank)))
    merged['murcko_scaffold']=merged.canonical_smiles.map(lambda s:MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s)))
    merged['selected']=merged.candidate_id.isin(selected.candidate_id)
    merged['stage3d_rank']=merged.stage3d_priority
    merged['selection_reason']=np.where(merged.selected,'structural gate PASS; top-3 distinct Murcko scaffolds by frozen M4-safe final_rank',merged.structural_gate_reason)
    selected=merged[merged.selected].sort_values('final_rank').copy()
    assert selected.candidate_id.tolist()==expected[:3]
    merged.to_csv(OUT/'stage3d_ranked_candidates.csv',index=False)
    selected.to_csv(OUT/'stage3d_selected_candidates.csv',index=False)
    validation=dict(status='PASS',identity=identity,merged_population=200,gate_unchanged=True,router_order_unchanged=True,selected_count=3,selected_candidates=selected.candidate_id.tolist(),scaffold_selection_matches_frozen_rule=True,pose_selection_matches_frozen_rule=True,selected_pose_graph_identity_pass=True,historical_drugclip_contamination=False)
    write('stage3d_validation.json',validation)
    manifest=dict(schema='pacer.stage3d.corrected.freeze.v1',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),canonical_input=rec(ROOT/'project/results/pacer_candidates_v01/predock_portfolio.csv'),drugclip_input=rec(ROUTER),docking_inputs=[rec(DOCK/n) for n in ['metadata.json','ledger.jsonl','framewise_scores.csv','statewise_scores.csv']],identity=identity,structural_gate=policy['structural_gate']['thresholds'],selection_rule=policy['selection']['rule'],selected_count=3,pose_selection_rule=policy['selection']['pose_freeze_rule'],frozen_poses=frozen,policy_sources=[rec(p) for p in originals]+[rec(OLD/'STAGE3D_FREEZE_MANIFEST.json')],policy_source_note='Historical Stage3D artifacts used only as frozen policy/implementation evidence, never as DrugCLIP ranking input.',historical_drugclip_ranking_used=False)
    write('stage3d_selection_manifest.json',manifest)
    assert all(sha(Path(p))==h for p,h in preserved.items())
    for p in OUT.glob('*.json'):
        assert '3df43a327d085bbe883840810f9876bddd403d799cb3f381d26294f54592bb29' not in p.read_text()
    write('stage3d_provenance_audit.json',dict(status='PASS',canonical_input=manifest['canonical_input'],drugclip_input=manifest['drugclip_input'],docking_inputs=manifest['docking_inputs'],identity=identity,selection_policy=manifest['selection_rule'],selected_count=3,config='No separate config; frozen constants preserved in copied implementation',implementation=[rec(p) for p in [Path(__file__),OUT/'identity_contract.py',*amended]],source_implementation=manifest['policy_sources'],corrections=['corrected input and output wiring','full molecular graph plus candidate ID identity gate','reporting priority mapped by candidate identity instead of integer dataframe index'],historical_drugclip_bundle_used=False,old_stage3d_files_unchanged=True,outputs=[rec(p) for p in sorted(OUT.rglob('*')) if p.is_file()],training=False,redocking=False,drugclip_inference=False,router_rerun=False))
    print(json.dumps(validation,indent=2))

if __name__=='__main__':
    main()
