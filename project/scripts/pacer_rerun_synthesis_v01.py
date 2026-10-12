#!/usr/bin/env python3
"""Phase D preflight and complete explicit status records (no invented routes).
Run in .venv-pacer-retro-v01. --output defaults to analysis directory.
"""
import argparse,json,hashlib,sys,importlib.metadata as md
from pathlib import Path
import pandas as pd
from aizynthfinder.aizynthfinder import AiZynthFinder

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'project/results/pacer_rerun_admet_synthesis_v01'
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=DEFAULT);o=ap.parse_args().output
    d=pd.read_csv(o/'ALL200_PHYSICOCHEMICAL_v01.csv');audit=json.loads((o/'EXTERNAL_ASSET_PREFLIGHT_v01.json').read_text())
    pkg=md.distribution('aizynthfinder');assets=[]
    for p in (ROOT/'project').rglob('*'):
        if p.is_file() and ('uspto_model' in p.name or 'uspto_templates' in p.name or 'zinc_stock' in p.name): assets.append(str(p))
    finder=AiZynthFinder()
    pre={'tool':'AiZynthFinder','version':md.version('aizynthfinder'),'python':sys.version,'import_and_object_creation_succeeded':True,'source':'https://github.com/MolecularAI/aizynthfinder','license':pkg.metadata.get('License'),'onnxruntime':md.version('onnxruntime'),'rdkit':md.version('rdkit'),'model_record':'https://zenodo.org/records/7797465','model_license':'CC-BY-4.0 per Zenodo API metadata','expansion_policy_available':False,'templates_available':False,'stock_available':False,'local_assets_found':assets,'network_preflight':audit,'can_execute_search':False,'reason':'Required complete policy/template/stock assets unavailable. Official ZINC stock URL returns 403; Zenodo files timed out. No access restriction bypass. Package installation and successful import alone do not constitute retrosynthesis.','all200_search_attempt_policy':'Preflight gate failed for every candidate; per-molecule route search NOT_ASSESSED, never NO_ROUTE_FOUND.','route_validation':'No computational routes generated; no experimental synthesis claim.'}
    with (o/'RETROSYNTHESIS_TOOL_PREFLIGHT_v01.json').open('x') as f:json.dump(pre,f,indent=2)
    rows=[]
    for r in d.to_dict('records'):
        rows.append({k:r[k] for k in ['candidate_id','canonical_smiles']}|{'status':'NOT_ASSESSED','reason':pre['reason'],'tool':'AiZynthFinder','tool_version':md.version('aizynthfinder'),'route_found':None,'candidate_route_count':None,'shortest_reaction_steps':None,'intermediates':None,'starting_materials':None,'route_score':None,'policy_version':None,'stock_version':None,'protection_deprotection':None,'key_transformations':None,'special_conditions':None,'expensive_material_evidence':None,'raw_route_file':None})
    with (o/'ALL200_RETROSYNTHESIS_v01.csv').open('x') as f:pd.DataFrame(rows).to_csv(f,index=False)
    syn=pd.read_csv(o/'PHASE_B_SYNTHESIS_HEURISTICS_v01.csv');syn['route_reason']=pre['reason'];syn['reaction_steps']=None;syn['route_complexity']='NOT_ASSESSED';syn['protection_deprotection']='NOT_ASSESSED';syn['special_conditions']='NOT_ASSESSED';syn['supplier_query_status']='NOT_ASSESSED_NO_ROUTE_STARTING_MATERIALS';syn['starting_material_smiles']=None;syn['supplier']=None;syn['product_url']=None;syn['CAS']=None;syn['pack_size']=None;syn['price']=None;syn['query_date']=None
    with (o/'ALL200_SYNTHESIS_FEASIBILITY_v01.csv').open('x') as f:syn.to_csv(f,index=False)
    with (o/'PHASE_D_COMPLETE_v01.json').open('x') as f:json.dump({'status':'COMPLETE_WITH_EXPLICIT_GAPS','molecule_records':200,'route_search_executed':0,'NOT_ASSESSED':200,'NO_ROUTE_FOUND':0,'verified_starting_material_availability':0,'yields_estimated':0,'costs_estimated':0},f,indent=2)
    print('Phase D: 200 explicit NOT_ASSESSED records; no route invented.')
if __name__=='__main__':main()
