"""Fail closed on canonical/router/docking parent identity before Stage3D joins."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

CANONICAL_HASH = '0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b'
ROUTER_HASH = 'bda22c9354a9ebf4a263ce92ccf1257c8b27e8e3e747404b0f4fa97c3fb34d89'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def graph(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError('Invalid SMILES')
    return Chem.MolToSmiles(rdMolStandardize.Cleanup(mol), canonical=True, isomericSmiles=True)

def state_parent(smiles):
    mol = rdMolStandardize.Cleanup(Chem.MolFromSmiles(str(smiles)))
    return Chem.MolToSmiles(rdMolStandardize.Uncharger().uncharge(mol), canonical=True, isomericSmiles=True)

def validate_inputs(root, dock, router_path):
    canonical = root / 'project/results/pacer_candidates_v01/predock_portfolio.csv'
    assert sha(canonical) == CANONICAL_HASH and sha(router_path) == ROUTER_HASH
    table = pd.read_csv(canonical)
    router = pd.read_csv(router_path)
    table['_graph'] = table.canonical_smiles.map(graph)
    router['_graph'] = router.canonical_smiles.map(graph)
    for df in [table, router]:
        assert len(df) == df.candidate_id.nunique() == df._graph.nunique() == 200
    assert set(table._graph) == set(router._graph)
    pairs = table[['candidate_id','_graph']].merge(router[['candidate_id','_graph']], on='_graph', validate='one_to_one', suffixes=('_canonical','_router'))
    assert pairs.candidate_id_canonical.eq(pairs.candidate_id_router).all()
    assert router.pair_id.eq(router.candidate_id).all()
    meta = json.loads((dock/'metadata.json').read_text())
    assert meta['source'].replace('\\','/') == 'project/results/pacer_candidates_v01/predock_portfolio.csv'
    assert meta['n_molecules'] == 200 and meta['clusters'] == list(range(10))
    ledger = pd.read_json(dock/'ledger.jsonl', lines=True)
    assert len(ledger) == 2000 and ledger.molecule_id.nunique() == 200
    assert set(ledger.molecule_id) == set(table.candidate_id)
    assert (ledger.error.isna() | ledger.error.astype(str).str.strip().eq('')).all()
    assert not ledger.duplicated(['molecule_id','cluster']).any()
    assert ledger.groupby('molecule_id').cluster.apply(lambda x:set(x)==set(range(10))).all()
    assert np.isfinite(ledger[['vina_affinity','pose_centroid_x','pose_centroid_y','pose_centroid_z','native_pocket_coverage','n_contacted_residues']].to_numpy(float)).all()
    parents = dict(zip(table.candidate_id, table._graph))
    for row in ledger[['molecule_id','state_smiles']].drop_duplicates().itertuples(index=False):
        assert state_parent(row.state_smiles) == state_parent(parents[row.molecule_id]), row.molecule_id
    return dict(canonical_population=200,docking_population=200,drugclip_population=200,canonical_docking_graph_overlap=200,canonical_drugclip_graph_overlap=200,docking_drugclip_graph_overlap=200,missing=0,extra=0,duplicates=0,docking_tasks_expected=2000,docking_tasks_completed=2000,receptor_clusters=10,state_parent_mapping='Full isomeric graph after explicit protonation neutralization; no truncated InChIKey matching',canonical_input_sha256=sha(canonical),drugclip_input_sha256=sha(router_path))
