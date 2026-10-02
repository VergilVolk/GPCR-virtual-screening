"""Prospective input preparation and monitoring; launch is explicitly user-only.

No membrane builder is substituted for the frozen GaMD receptors. Missing
cluster-specific systems are blockers, never aliases to the historical 7TRS box.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PROJECT = REPO / 'project'
FROZEN = PROJECT / 'results/pacer_stage3d_corrected_20261002_v01'
SEEDS = [27101, 38201, 49301]
SOURCE_HASHES = {
    'stage3d_selected_candidates.csv': '5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3',
    'md_shortlist_pose_manifest.csv': '8ef3c4be8ff3800543d09e3b37c11331671e961c7bd3f29cc70850b787e311a6',
}
RUNNER = PROJECT / 'scripts/run_cm00734_stage_b_production_md_v01.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def inputs():
    for name, digest in SOURCE_HASHES.items():
        assert sha(FROZEN / name) == digest, f'Frozen input changed: {name}'
    with (FROZEN / 'md_shortlist_pose_manifest.csv').open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    assert [r['candidate_id'] for r in rows] == ['PACER0010', 'PACER0073', 'PACER0027']
    for row in rows:
        assert sha(REPO / row['pose_file']) == row['pose_file_sha256']
    return rows


def prepare():
    from rdkit import Chem
    manifest_path=HERE/'master_manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()).get('unique_systems')==12:
        raise SystemExit('Use build_pacer_stage4_prospective_membrane.py --all --dry-run; legacy prepare-inputs must not reset the new system namespace')
    if any((HERE/'production').rglob('progress.json')) or (HERE/'preflight_receipt.json').exists():
        raise SystemExit('Refusing to overwrite campaign inputs after preflight/production receipts exist')
    rows = inputs()
    ligand_dir = HERE / 'ligand_inputs'
    ligand_dir.mkdir(exist_ok=True)
    records = []
    jobs = []
    for row in rows:
        # MODEL 1 is the Stage3D selected pose, not a new pose-selection rule.
        lines = (REPO / row['pose_file']).read_text().split('ENDMDL', 1)[0].splitlines()
        smiles = next(x[len('REMARK SMILES '):] for x in lines
                      if x.startswith('REMARK SMILES ') and not x.startswith('REMARK SMILES IDX '))
        molecule = Chem.MolFromSmiles(smiles)
        expected = Chem.MolFromSmiles(row['source_state_smiles'])
        canonical = lambda m: Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
        assert canonical(molecule) == canonical(expected), 'Frozen microstate mismatch'
        pairs = []
        for line in lines:
            if line.startswith('REMARK SMILES IDX '):
                values = list(map(int, line[len('REMARK SMILES IDX '):].split()))
                pairs.extend(zip(values[::2], values[1::2]))
        coordinates = {int(x[6:11]): tuple(float(x[a:b]) for a, b in [(30,38),(38,46),(46,54)])
                       for x in lines if x.startswith(('ATOM', 'HETATM'))}
        assert len(pairs) == molecule.GetNumAtoms()
        assert {a for a, _ in pairs} == set(range(1, molecule.GetNumAtoms()+1))
        conf = Chem.Conformer(molecule.GetNumAtoms())
        for a, serial in pairs:
            conf.SetAtomPosition(a-1, coordinates[serial])
        molecule.AddConformer(conf)
        explicit = Chem.AddHs(molecule, addCoords=True)
        # Existing polar H coordinates are kept; only omitted hydrogens are constructed.
        parents = []
        for line in lines:
            if line.startswith('REMARK H PARENT '):
                values = list(map(int, line[len('REMARK H PARENT '):].split()))
                parents.extend(zip(values[::2], values[1::2]))
        consumed = set()
        for parent, serial in parents:
            hs = [a.GetIdx() for a in explicit.GetAtomWithIdx(parent-1).GetNeighbors()
                  if a.GetAtomicNum() == 1 and a.GetIdx() not in consumed]
            assert hs, 'Polar hydrogen-parent mapping failed'
            explicit.GetConformer().SetAtomPosition(hs[0], coordinates[serial])
            consumed.add(hs[0])
        for a, serial in pairs:
            assert tuple(explicit.GetConformer().GetAtomPosition(a-1)) == coordinates[serial]
        explicit.SetProp('_Name', row['candidate_id'])
        explicit.SetProp('source_pose_sha256', row['pose_file_sha256'])
        explicit.SetProp('source_model', '1')
        explicit.SetProp('source_state_smiles', smiles)
        sdf = ligand_dir / (row['candidate_id'] + '_frozen_state.sdf')
        writer = Chem.SDWriter(str(sdf)); writer.write(explicit); writer.close()
        loaded = next(iter(Chem.SDMolSupplier(str(sdf), removeHs=False)))
        # SDF readers infer chirality from 3D coordinates. The source SMILES
        # has no defined stereo; do not promote that inference to identity.
        assert not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in molecule.GetAtoms())
        assert not any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in molecule.GetBonds())
        Chem.RemoveStereochemistry(loaded)
        assert canonical(Chem.RemoveHs(loaded)) == canonical(molecule)
        for a, serial in pairs:
            observed = tuple(loaded.GetConformer().GetAtomPosition(a-1))
            assert max(abs(x-y) for x,y in zip(observed,coordinates[serial])) < 0.0001
        receptor = PROJECT / f"results/m4_gamd_ensemble/receptors/cluster_{int(row['source_cluster']):02d}_receptor.pdb"
        records.append(dict(row, receptor_file=str(receptor.relative_to(REPO)).replace('\\','/'),
                            receptor_sha256=sha(receptor), sdf_file=str(sdf.relative_to(REPO)).replace('\\','/'),
                            sdf_sha256=sha(sdf), state_formal_charge=Chem.GetFormalCharge(molecule),
                            state_inchikey=Chem.MolToInchiKey(molecule),
                            stereo_identity_source='source_state_smiles; SDF geometry does not define unspecified stereo',
                            heavy_coordinates_unchanged=True, parameterization_status='NOT_STARTED'))
        for context in ['apo', 'probe_only', 'candidate_no_probe', 'candidate_probe']:
            system = f"{row['candidate_id']}__c{int(row['source_cluster']):02d}__{context}"
            for replica, seed in enumerate(SEEDS,1):
                output = f'production/{system}/replica_{replica:02d}'
                jobs.append(dict(candidate_id=row['candidate_id'], receptor_cluster=row['source_cluster'],
                    context=context, system=system, replica=replica, seed=seed,
                    production_length_ns=10, production_steps=5000000,
                    pose_sha256=row['pose_file_sha256'] if context.startswith('candidate') else '',
                    system_path=f'systems/{system}', output_path=output,
                    checkpoint_path=output+'/checkpoint.chk', status='BLOCKED_SYSTEM_BUILD'))
    assert len(jobs) == 36
    with (HERE/'job_matrix.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(jobs[0])); writer.writeheader(); writer.writerows(jobs)
    evidence = [RUNNER, PROJECT/'scripts/build_pacer_dc_membrane_reference.py',
        PROJECT/'scripts/build_cm00734_stage_b_membrane_v01.py',
        PROJECT/'scripts/run_pacer_dc_short_equilibration.py',
        PROJECT/'scripts/run_pacer_dc_restraint_release.py',
        PROJECT/'docs/PACER_DC_CLOSE_LOOP_CONTRACT_v01.md',
        PROJECT/'docs/PACER_DC_CLOSE_LOOP_20NS_EXECUTION_PLAN_v01.md',
        PROJECT/'environment_pacer_dc_md.yml']
    manifest = dict(protocol_id=HERE.name, execution_platform='Linux GPU server',
        server_repository='/data/GPCR-virtual-screening', conda_env='pacer-dc-md',
        production_length_ns=10, timestep_fs=2, production_steps=5000000,
        replicas_per_context=3, paired_seeds=SEEDS, trajectory_ps=50, checkpoint_ps=100,
        temperature_K=300, pressure_bar=1, thermostat='LangevinMiddle', friction_per_ps=1,
        forcefields=['amber14/protein.ff14SB.xml','amber14/lipid17.xml','amber14/tip3p.xml',
                     'openff_unconstrained-2.2.1.offxml'], ligand_charge_policy='AM1-BCC via AmberTools',
        membrane='POPC', ionic_strength_molar=0.15, membrane_minimum_padding_nm=1,
        nonbonded=dict(method='PME',cutoff_nm=1,ewald_error_tolerance=0.0005),
        constraints='HBonds', rigid_water=True,
        membrane_barostat=dict(mode='XYIsotropic/ZFree',frequency=25,surface_tension_bar_nm=0),
        protein_heavy_restraint_k_initial=1000, production_restraint_k=0,
        historical_equilibration=dict(smoke_steps=500,smoke_timestep_fs=0.5,
            nvt_ps=1,npt_ps=1,timestep_fs=0.5,restraint_schedule=[500,100,10,0],
            restraint_stage_ps=0.5,restraint_timestep_fs=1),
        probe_smiles='CC(=O)OCC[N+](C)(C)C',
        control_policy='Cluster-matched four contexts; historical 7TRS controls are not eligible',
        control_policy_evidence='CLOSE_LOOP_CONTRACT_v01 section 3: incompatible lineage requires own four contexts',
        apo_reuse='NO',probe_only_reuse='NO',expected_job_count=36,
        source_hashes=SOURCE_HASHES,candidates=records,
        job_matrix_sha256=sha(HERE/'job_matrix.csv'),
        evidence={str(p.relative_to(REPO)).replace('\\','/'):sha(p) for p in evidence},
        system_build_status='PARTIAL_INPUTS_ONLY',preflight_status='NOT_RUN',
        ready_to_launch=False,production_started=False,
        blockers=['No validated builder accepting the frozen GaMD 274-residue construct (historical builder expects 270 receptor residues plus other protein chains)',
                  'Cluster-to-OPM residue/atom mapping and rigid transform are not frozen; 7TRS_OPM.pdb is absent in this checkout',
                  'Cluster-matched ACh probe coordinates/atom mapping are not frozen; no probe insertion by guesswork',
                  '12 membrane systems, 36 matched equilibration/release receipts and Linux-server CUDA/system preflight are missing'])
    dump(HERE/'master_manifest.json',manifest)
    print(json.dumps(dict(ligand_inputs=3,jobs=36,ready_to_launch=False)))


def status():
    with (HERE/'job_matrix.csv').open(newline='',encoding='utf-8') as f:
        jobs=list(csv.DictReader(f))
    counts={}; results=[]
    for row in jobs:
        out=HERE/row['output_path']; ledger=out/'progress.json'
        progress=json.loads(ledger.read_text()) if ledger.exists() else {}
        state=progress.get('status',row['status'])
        if (out/'failure.json').exists(): state='failed'
        counts[state]=counts.get(state,0)+1
        ns=float(progress.get('completed_ns',0))
        chk=out/'checkpoint.chk'
        results.append(dict(system=row['system'],replica=int(row['replica']),status=state,
            completed_ns=ns,step=round(ns*500000),checkpoint_exists=chk.exists(),
            last_checkpoint_mtime=chk.stat().st_mtime if chk.exists() else None,
            note='running is last ledger state, not a live process check'))
    print(json.dumps(dict(expected_jobs=len(jobs),counts=counts,jobs=results),indent=2))


def preflight():
    inputs()
    manifest=json.loads((HERE/'master_manifest.json').read_text(encoding='utf-8'))
    assert sha(HERE/'job_matrix.csv') == manifest['job_matrix_sha256']
    # Do not replace this gate with file existence or a tiny CUDA-only test.
    raise SystemExit('BLOCKED: '+ '; '.join(manifest['blockers']))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('action',choices=['prepare-inputs','status','preflight'])
    action=parser.parse_args().action
    {'prepare-inputs':prepare,'status':status,'preflight':preflight}[action]()
