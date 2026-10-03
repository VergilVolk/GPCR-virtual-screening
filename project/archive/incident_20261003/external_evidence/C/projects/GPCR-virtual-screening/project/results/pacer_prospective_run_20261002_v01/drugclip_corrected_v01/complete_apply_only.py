"""Apply frozen DrugCLIP models to the canonical docked population; no fitting."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path('/mnt/c/projects/GPCR-virtual-screening')
HISTORY = Path('/mnt/c/projects/GPCR-virtual-screening-final')
OUT = ROOT / 'project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01'
STAGE = HISTORY / 'project/results/pacer_prospective_run_20261002_v01/stage3'
INPUT = ROOT / 'project/results/pacer_candidates_v01/predock_portfolio.csv'
OLD = ROOT / 'project/artifacts/drugclip2023_m4_loto_pacer200_v01'
SEEDS = [20260925, 20260926, 20260927]
EXPECTED_INPUT = '0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b'
EXPECTED_REPS = 'f0a71918b707c13053444a7de9f41d47487201b92f8b79f5ccc5940f1fdcef3a'
EXPECTED_2026 = '997d343b0ba4b9436a30d76b7b3a0281f99d7e7f78f2e1fbdef5dc3bb0c67a31'
EXPECTED_HEADS = ['d33b31d69c4437327b41c939c57521e27f56abb7491f4b75d0ef86d1533e58fd', '76d21a31436bdd51583b4f160ccca8ce7383af3848232bf14925a996a555c9a7', '9cb96a2ce47d04180c526760791abfd19df3bf2f2db69a07dc306aea0f478a8b']

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def save(name, data):
    p = OUT / name
    if p.exists():
        raise RuntimeError('Refusing overwrite: ' + str(p))
    p.write_text(json.dumps(data, indent=2) + '\n')

def graph(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None:
        raise ValueError('Invalid molecule')
    return Chem.MolToSmiles(rdMolStandardize.Cleanup(m), canonical=True, isomericSmiles=True)

def validate(frame):
    assert len(frame) == 200
    assert frame.candidate_id.nunique() == 200 and frame.canonical_smiles.nunique() == 200
    gs = frame.canonical_smiles.map(graph)
    assert gs.nunique() == 200 and set(gs) == set(TABLE.canonical_smiles.map(graph))
    original = TABLE.set_index('candidate_id').canonical_smiles.map(graph)
    assert all(original[cid] == g for cid, g in zip(frame.candidate_id, gs))
    return dict(rows=200, graph_overlap=200, missing=0, extra=0, duplicates=0)

def ordered_reps(archive):
    ids = list(map(str, archive['molecule_ids']))
    keys = [graph(s) for s in ids]
    assert len(keys) == len(set(keys)) == 200
    lookup = dict(zip(keys, range(200)))
    assert set(keys) == set(TABLE.canonical_smiles.map(graph))
    values = archive['molecule_representations'][[lookup[graph(s)] for s in TABLE.canonical_smiles]]
    assert values.shape == (200, 512) and np.isfinite(values).all()
    return values.astype(np.float32)

def record(p):
    return dict(path=str(p), sha256=sha(p))

def run(command, logfile):
    with (OUT / logfile).open('x') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('Command failed; see ' + logfile)

assert sha(INPUT) == EXPECTED_INPUT
TABLE = pd.read_csv(INPUT)
validate(TABLE)
CANON = record(INPUT)
torch.set_num_threads(4)

def task_a():
    reps_path = STAGE / 'embeddings.npz'
    assert sha(reps_path) == EXPECTED_REPS
    base = HISTORY / 'project/tools/DrugCLIP/artifacts/checkpoint_best.pt'
    assert sha(base) == 'dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e'
    manifest = json.loads((OLD / 'bundle_manifest.json').read_text())
    frozen = OLD / manifest['representations_file']
    assert sha(frozen) == manifest['representations_sha256']
    with np.load(frozen, allow_pickle=False) as z:
        pocket = z['m4_pocket_representation'].copy().astype(np.float32)
    with np.load(reps_path, allow_pickle=False) as z:
        reps = ordered_reps(z)
    assert pocket.shape == (512,) and np.isfinite(pocket).all()
    adapters = []
    for entry, seed in zip(manifest['adapters'], SEEDS):
        path = OLD / entry['file']
        assert sha(path) == entry['sha256']
        payload = torch.load(path, map_location='cpu')
        assert payload['base_seed'] == seed and payload['held_target'] == 'M4R'
        assert payload['mol_project']['linear1.weight'].shape[1] == 512
        adapters.append(dict(file=str(path), sha256=sha(path)))
    archive_path = OUT / 'prospective200_2023_frozen_representations.npz'
    assert not archive_path.exists()
    np.savez_compressed(archive_path, candidate_ids=TABLE.candidate_id.to_numpy(dtype=str), canonical_smiles=TABLE.canonical_smiles.to_numpy(dtype=str), molecule_representations=reps, m4_pocket_representation=pocket)
    pocket_hash = hashlib.sha256(pocket.tobytes(order='C')).hexdigest()
    save('bundle_manifest.json', dict(protocol_id=manifest['protocol_id'], representations_file=archive_path.name, representations_sha256=sha(archive_path), adapters=adapters, candidate_rows=200, canonical_input=CANON, candidate_representations=record(reps_path), pocket_source=dict(**record(frozen), array_key='m4_pocket_representation', pocket_id='M4R_cluster0', float32_c_order_sha256=pocket_hash)))
    scores = OUT / 'prospective200_2023_gpcr_loto_scores.csv'
    run([sys.executable, '-B', str(ROOT / 'project/scripts/score_pacer200_drugclip2023_m4_loto.py'), '--bundle-dir', str(OUT), '--output', str(scores)], '2023_scoring_stdout.log')
    df = pd.read_csv(scores).rename(columns={'pair_id': 'candidate_id'})
    validate(df)
    assert np.isfinite(df[[f'seed{s}' for s in SEEDS] + ['score_2023_gpcr_loto']]).all().all()
    df.to_csv(scores, index=False)
    audit = dict(canonical_input=CANON, embeddings=record(reps_path), backbone_checkpoint=record(base), adapters=adapters, pocket_representation_source=record(frozen), pocket_representation_sha256=pocket_hash, pocket_hash_encoding='float32 C-order bytes', protocol_id=manifest['protocol_id'], seeds=SEEDS, aggregation='mean of three normalized projected cosine scores; seed sd ddof=0', validation=validate(df), finite=True, score_range=[float(df.score_2023_gpcr_loto.min()), float(df.score_2023_gpcr_loto.max())], output=record(scores))
    # The stock scoring audit predates the candidate_id column rename; refresh its output hash.
    stock = scores.with_suffix('.audit.json')
    prior = json.loads(stock.read_text()); prior['output_sha256'] = sha(scores); prior['canonical_input'] = CANON
    stock.write_text(json.dumps(prior, indent=2) + '\n')
    save('2023_provenance_audit.json', audit)
    print('2023_GPCR_LOTO_COMPLETE', flush=True)

def task_b():
    checkpoint = ROOT / 'project/data/external/drugclip_science2026/litpcba_identity_90.pt'
    assert sha(checkpoint) == EXPECTED_2026
    pockets = [STAGE / f'{p}_pocket.lmdb' for p in ['7TRQ', '7TRP', '7TRS']]
    mols = STAGE / 'molecules.lmdb'
    env = lmdb.open(str(mols), subdir=False, readonly=True, lock=False)
    with env.begin() as txn:
        assert txn.stat()['entries'] == 200
    env.close()
    checkout = HISTORY / 'project/tools/DrugCLIP'
    unicore = HISTORY / 'project/tools/Uni-Core'
    assert (checkout / '.git/HEAD').read_text().strip() == '7a3a3fa33673f8668c811790f2e4681c98af44ef'
    assert (unicore / '.git/HEAD').read_text().strip() == '44f6386f4dcd7137fc1e5d5e768117d635d64a26'
    extraction = ROOT / 'project/scripts/extract_drugclip_embeddings_cpu.py'
    archive_path = OUT / 'prospective200_science2026_representations.npz'
    assert not archive_path.exists()
    command = [sys.executable, '-B', str(extraction), '--drugclip', str(checkout), '--unicore', str(unicore), '--checkpoint', str(checkpoint), '--molecules', str(mols), '--pockets', *map(str, pockets), '--output', str(archive_path), '--batch-size', '8', '--progress-every', '5']
    print('SCIENCE2026_BACKBONE_START', flush=True)
    run(command, '2026_extraction_stdout.log')
    extraction_audit = json.loads(archive_path.with_suffix('.audit.json').read_text())
    assert not extraction_audit['missing_checkpoint_keys'] and not extraction_audit['unexpected_checkpoint_keys']
    with np.load(archive_path, allow_pickle=False) as z:
        reps = torch.as_tensor(ordered_reps(z))
        assert list(map(str,z['pocket_ids'])) == ['7TRQ_M4_allosteric','7TRP_M4_allosteric','7TRS_M4_allosteric']
        pocket = torch.as_tensor(z['pocket_representations'].astype(np.float32))
    assert pocket.shape == (3,512) and torch.isfinite(pocket).all()
    sys.path.insert(0, str(ROOT / 'project/scripts'))
    from score_pacer200_drugclip2023_m4_loto import project
    values, heads = [], []
    for seed, expected in zip(SEEDS, EXPECTED_HEADS):
        path = ROOT / f'project/results/drugclip_science2026/family_aug_v01/science2026_13target_ep80_famaug.seed{seed}.projection.pt'
        assert sha(path) == expected
        payload = torch.load(path, map_location='cpu')
        assert payload['seed'] == seed and 'science2026_90.projection.pt' in payload['base_projection']
        with torch.inference_mode():
            values.append((project(reps, payload['mol_project']) @ project(pocket, payload['pocket_project']).T).max(dim=1).values.numpy())
        heads.append(record(path))
    matrix = np.stack(values, axis=1)
    assert np.isfinite(matrix).all()
    df = TABLE[['candidate_id','canonical_smiles']].copy()
    for i, seed in enumerate(SEEDS): df[f'seed{seed}'] = matrix[:,i]
    df['score_2026_family_aug'] = matrix.mean(axis=1)
    df['score_2026_seed_sd'] = matrix.std(axis=1,ddof=0)
    df = df.sort_values('score_2026_family_aug', ascending=False, kind='mergesort')
    df['rank_2026_family_aug'] = range(1,201)
    output = OUT / 'prospective200_2026_family_aug_scores.csv'
    assert not output.exists(); df.to_csv(output,index=False)
    save('2026_provenance_audit.json',dict(canonical_input=CANON, backbone_checkpoint=record(checkpoint), representations=record(archive_path), extraction_code=record(extraction), extraction_command=command, molecule_lmdb=record(mols), pocket_lmdbs=[record(p) for p in pockets], checkout_commit=(checkout/'.git/HEAD').read_text().strip(), unicore_commit=(unicore/'.git/HEAD').read_text().strip(), dictionaries=[record(checkout/'data/dict_mol.txt'),record(checkout/'data/dict_pkt.txt')], heads=heads, seeds=SEEDS, protocol='Science2026 frozen backbone + frozen family-aug projections; updated prospective three crystallographic M4 pockets', pocket_ids=['7TRQ_M4_allosteric','7TRP_M4_allosteric','7TRS_M4_allosteric'], pocket_population_note='Uses the three pockets declared in the updated prospective mainline, rather than reproducing the historical ten GaMD pocket population.', aggregation='maximum normalized projected cosine over declared pockets per seed, then three-seed arithmetic mean', validation=validate(df), finite=True, score_range=[float(df.score_2026_family_aug.min()),float(df.score_2026_family_aug.max())], output=record(output)))
    print('2026_FAMILY_AUG_COMPLETE', flush=True)

def task_cd():
    a = pd.read_csv(OUT/'prospective200_2023_gpcr_loto_scores.csv')
    b = pd.read_csv(OUT/'prospective200_2026_family_aug_scores.csv')
    validate(a);validate(b)
    a['_graph'] = a.canonical_smiles.map(graph); b['_graph'] = b.canonical_smiles.map(graph)
    merged = a.merge(b[['_graph','candidate_id','score_2026_family_aug']],on='_graph',validate='one_to_one',suffixes=('','_2026'))
    assert (merged.candidate_id == merged.candidate_id_2026).all()
    # Keep canonical input order for the router's frozen method='first' tie rule.
    df = TABLE[['candidate_id','canonical_smiles']].merge(merged[['candidate_id','score_2023_gpcr_loto','score_2026_family_aug']],on='candidate_id',validate='one_to_one')
    df['pair_id'] = df.candidate_id; df['target'] = 'M4R'
    df['score_2026_13t_famaug'] = df.score_2026_family_aug
    validate(df)
    inp = OUT/'prospective200_two_model_scores.csv'; assert not inp.exists();df.to_csv(inp,index=False)
    script = ROOT/'project/scripts/pacer_drugclip_router.py'
    assert sha(script) == '22ec38f18ed96e20bf726602e3997cc1aad8e1cb13127ddb1279d459f11a290e'
    output = OUT/'prospective200_m4_safe_routed_ranking.csv'
    run([sys.executable,'-B',str(script),'--input',str(inp),'--output',str(output)],'router_stdout.log')
    routed = pd.read_csv(output); validate(routed)
    assert routed.route.eq('m4_2023_gpcr_loto').all() and routed.frozen_protocol.all()
    assert np.isfinite(routed[['pacer_binding_score','score_2023_gpcr_loto','score_2026_13t_famaug']]).all().all()
    save('router_provenance_audit.json',dict(canonical_input=CANON,input=record(inp),output=record(output),router_script=record(script),protocol='pacer_drugclip_m4_safe_router_v01',weights=dict(old=0.5,new=0.5),m4_policy='2023 percentile only',route_counts=routed.route.value_counts().to_dict(),finite=True,validation=validate(routed),source_scores=[record(OUT/'prospective200_2023_gpcr_loto_scores.csv'),record(OUT/'prospective200_2026_family_aug_scores.csv')]))
    save('COMPLETION_RECEIPT.json',dict(status='COMPLETE',canonical_input=CANON,backbone_2023_reused=True,validation=validate(routed),outputs=[record(OUT/n) for n in ['prospective200_2023_gpcr_loto_scores.csv','prospective200_2026_family_aug_scores.csv','prospective200_two_model_scores.csv','prospective200_m4_safe_routed_ranking.csv']],training=False,redocking=False,stage3d_executed=False))
    print('M4_SAFE_ROUTER_COMPLETE',flush=True)

if __name__ == '__main__':
    try:
        {'a':task_a,'b':task_b,'cd':task_cd}[sys.argv[1]]()
    except Exception as exc:
        save('FAILURE_'+sys.argv[1]+'.json',dict(status='BLOCKED',error=repr(exc),canonical_input=CANON))
        raise
