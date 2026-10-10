"""Build an isolated, authenticated Stage4 rerun deployment artifact."""
import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROJECT = REPO/'project'
OUT = PROJECT/'results/pacer_stage4_rerun_v01'
ROOT = OUT/'server_bundle/pacer_stage4_rerun_v01'
STAGE3 = PROJECT/'results/pacer_xr_rerun_v01'
RANK_HASH = '6fff7e34d94466ff288b52281b1c3a64fef4ea6cd9da84ca45bc703429a04559'
IDS = ['PACERGEN02432', 'PACERGEN00123', 'PACERGEN00462']
PROTECTED = ['pacer_xr_rerun_v01', 'pacer_stage4_prospective_10ns_v01',
             'pacer_stage4_prospective_fkg_v02_v01', 'pacer_dc_production_v01']


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()


def dump(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def readcsv(p):
    with p.open(encoding='utf-8-sig', newline='') as f: return list(csv.DictReader(f))


def writecsv(p, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def snapshot():
    return {str(p.relative_to(REPO)).replace('\\', '/'): [p.stat().st_size, p.stat().st_mtime_ns]
            for d in PROTECTED for p in (PROJECT/'results'/d).rglob('*') if p.is_file()}


def copy(src, dest):
    assert src.is_file(), f'Missing required resource: {src}'
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)
    assert sha(src) == sha(dest)


def initialize():
    assert not OUT.exists(), 'Refuse to overwrite existing rerun preparation'
    assert sha(STAGE3/'pacer_xr_rerun_full200_ranking.csv') == RANK_HASH
    rows = sorted(readcsv(STAGE3/'pacer_xr_rerun_full200_ranking.csv'), key=lambda r: int(r['candidate_cascade_rank']))[:3]
    assert [r['candidate_id'] for r in rows] == IDS
    assert [int(r['candidate_cascade_rank']) for r in rows] == [1, 2, 3]
    dump(OUT/'logs/protected_before.json', snapshot())
    for d in ['config', 'inputs/receptors', 'inputs/ligands', 'inputs/poses', 'inputs/probe',
              'inputs/forcefield', 'scripts/legacy', 'manifests', 'logs', 'production', 'systems']:
        (ROOT/d).mkdir(parents=True, exist_ok=True)
    requests = json.loads((PROJECT/'results/pacer_xr_stage3_ablation_v01/glide_pose_requests.json').read_text())
    cluster_scores = readcsv(STAGE3/'glide/glide_cluster_scores.csv')
    selected = []
    for row in rows:
        cid = row['candidate_id']
        scores = [r for r in cluster_scores if r['candidate_id'] == cid]
        assert len(scores) == 10 and len({r['cluster_id'] for r in scores}) == 10
        best = min(scores, key=lambda r: float(r['BE_i']))
        assert int(best['cluster_id']) == 1
        assert abs(float(best['BE_i'])-float(row['Glide_BEmin_raw'])) < 1e-8
        req = [r for r in requests if r['candidate_id'] == cid and r['context'] == 'BEmin']
        assert len(req) == 1
        req = req[0]
        assert req['receptor_id'] == 'cluster_01' and req['expected_smiles'] == row['canonical_smiles']
        assert req['selected_variant_id'] == best['selected_variant_id']
        assert req['entry_index'] == int(best['selected_entry_index'])
        pv = Path(req['poseviewer']); assert pv.is_file()
        selected.append(dict(candidate_id=cid, cascade_rank=int(row['candidate_cascade_rank']),
            canonical_smiles=row['canonical_smiles'], canonical_smiles_sha256=hashlib.sha256(row['canonical_smiles'].encode()).hexdigest(),
            source_poseviewer=str(pv), source_poseviewer_sha256=sha(pv),
            entry_index=req['entry_index'], selected_variant_id=req['selected_variant_id'],
            glide_gscore=float(best['glide_gscore']), BEmin=float(best['BE_i']), pmf=float(best['pmf']),
            cluster='cluster_01'))
    assert len({r['source_poseviewer'] for r in selected}) == 1
    copy(Path(selected[0]['source_poseviewer']), ROOT/'inputs/poses/original_cluster_01_pv.maegz')
    dump(ROOT/'manifests/extraction_requests.json', selected)
    copy(STAGE3/'pacer_xr_rerun_full200_ranking.csv', ROOT/'manifests/frozen_ranking.csv')
    copy(STAGE3/'glide/glide_cluster_scores.csv', ROOT/'manifests/glide_cluster_scores.csv')
    copy(STAGE3/'glide/glide_all_poses.csv', ROOT/'manifests/glide_all_poses.csv')
    assets = json.loads((PROJECT/'config/pacer_xr_assets_v01.json').read_text())['receptors']['cluster_01']
    copy(PROJECT/assets['source_receptor'], ROOT/'inputs/receptors/cluster_01_receptor.pdb')
    copy(PROJECT/assets['prepared_receptor'], ROOT/'inputs/receptors/docking_prepared_cluster_01.mae')
    archive = PROJECT/'results/pacer_stage4_prospective_10ns_v01/server_bundle_v2/project'
    copy(archive/'data/pdb/7TRS_OPM.pdb', ROOT/'inputs/receptors/7TRS_OPM.pdb')
    copy(archive/'data/pdb/m4_ligands/ACH_ideal.sdf', ROOT/'inputs/probe/ACH_ideal.sdf')
    copy(PROJECT/'data/pdb/7TRS.pdb', ROOT/'inputs/probe/7TRS.pdb')
    # Historical OPM-oriented ACh coordinates and the complete original mapping audit.
    backup = Path('C:/projects/PACER_STAGE4_MD_backup')
    copy(backup/'orientation/cluster0/ACH.sdf', ROOT/'inputs/probe/historical_OPM_ACH.sdf')
    copy(backup/'orientation/cluster0/alignment_audit.json', ROOT/'manifests/historical_ACH_alignment_audit.json')
    legacy = ['pacer_stage4_prospective.py', 'build_pacer_stage4_prospective_membrane.py',
              'build_pacer_dc_membrane_reference.py', 'build_pacer_dc_openmm_reference_systems.py',
              'run_pacer_stage4_membrane_smoke.py', 'run_pacer_dc_short_equilibration.py',
              'run_pacer_dc_restraint_release.py', 'run_pacer_dc_production_md.py',
              'run_cm00734_stage_b_production_md_v01.py']
    for name in legacy: copy(PROJECT/'scripts'/name, ROOT/'scripts/legacy'/name)
    copy(PROJECT/'docs/PACER_STAGE4_FINAL_HANDOFF_v01.md', ROOT/'manifests/PACER_STAGE4_FINAL_HANDOFF_v01.md')
    copy(PROJECT/'environment_pacer_dc_md.yml', ROOT/'environment-linux.yml')
    for name in ['build_pacer_stage4_rerun_v01.py', 'run_pacer_stage4_rerun_v01.py',
                 'schedule_pacer_stage4_rerun_v01.py', 'preflight_pacer_stage4_rerun_v01.py']:
        copy(PROJECT/'scripts'/name, ROOT/'scripts'/name)
    dump(ROOT/'config/pacer_stage4_rerun_v01.json', dict(protocol_id='pacer_stage4_rerun_v01',
        production_steps=5000000, production_ns=10, timestep_fs=2, trajectory_ps=10, checkpoint_ps=100,
        seeds=[27101,38201,49301], contexts={'A':[], 'P':['ACH'], 'C':['candidate'], 'CP':['candidate','ACH']},
        temperature_K=300, friction_per_ps=1, pressure_bar=1, initial_restraint_k=1000,
        forcefields=['amber14/protein.ff14SB.xml','amber14/lipid17.xml','amber14/tip3p.xml','openff_unconstrained-2.2.1.offxml'],
        charge_method='AM1-BCC', membrane='POPC', padding_nm=1, ionic_strength_M=0.15,
        PME_cutoff_nm=1, ewald_error_tolerance=5e-4, constraints='HBonds', rigid_water=True,
        barostat=dict(XY='XYIsotropic',Z='ZFree',frequency=25,surface_tension_bar_nm=0),
        smoke=dict(steps=500,timestep_fs=0.5,heating_K=[50,100,200,300]),
        equilibration=dict(nvt_ps=1,npt_ps=1,timestep_fs=0.5),
        release=dict(schedule=[500,100,10,0],stage_ps=0.5,timestep_fs=1),
        control_policy='12 distinct systems / 36 independently executed tasks; common cluster_01 receptor baseline; no shared trajectories',
        production_started=False))
    print('Initialized isolated rerun; extract original MAEGZ entries next')


def finalize():
    sys.dont_write_bytecode=True
    environment=ROOT/'environment-linux.yml'
    environment.write_text((PROJECT/'environment_pacer_dc_md.yml').read_text(encoding='utf-8').replace('name: pacer-dc-md','name: pacer-stage4-rerun-v01'),encoding='utf-8')
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from scipy.spatial import cKDTree
    sys.path.insert(0, str(ROOT/'scripts/legacy'))
    import pacer_stage4_prospective as old
    extraction=ROOT/'manifests/extracted_pose_records.json'
    if extraction.exists(): copy(extraction,OUT/'logs/extracted_pose_records.json')
    records = json.loads((OUT/'logs/extracted_pose_records.json').read_text())
    for name in ['build_pacer_stage4_rerun_v01.py','run_pacer_stage4_rerun_v01.py',
                 'schedule_pacer_stage4_rerun_v01.py','preflight_pacer_stage4_rerun_v01.py']:
        copy(PROJECT/'scripts'/name,ROOT/'scripts'/name)
    transform = old.alignment(old.ca_rows(ROOT/'inputs/receptors/cluster_01_receptor.pdb',' '),
                              old.ca_rows(ROOT/'inputs/receptors/7TRS_OPM.pdb','R'))
    source = (ROOT/'inputs/receptors/cluster_01_receptor.pdb').read_text().splitlines()
    oriented = []
    for line in source:
        if line.startswith(('ATOM','HETATM')):
            xyz = old.transform_xyz([[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]], transform)[0]
            line = line[:30]+''.join(f'{x:8.3f}' for x in xyz)+line[54:]
        oriented.append(line)
    normalized = old.normalized_receptor_lines(oriented)
    (ROOT/'inputs/receptors/oriented_cluster_01.pdb').write_text('\n'.join(normalized)+'\n')
    assert len(old.ca_rows(ROOT/'inputs/receptors/cluster_01_receptor.pdb',' ')) == 274
    # Recover the historical QC ACh pose into cluster_01 using the frozen two-transform chain.
    prior = json.loads((ROOT/'manifests/historical_ACH_alignment_audit.json').read_text())
    assert sha(ROOT/'inputs/probe/historical_OPM_ACH.sdf')==prior['transformed_probe_sha256']
    assert sha(ROOT/'inputs/receptors/7TRS_OPM.pdb')==prior['reference_sha256']
    ach = Chem.SDMolSupplier(str(ROOT/'inputs/probe/historical_OPM_ACH.sdf'), removeHs=False)[0]
    assert ach is not None and Chem.GetFormalCharge(ach) == 1
    assert old.graph(ach) == old.graph(Chem.MolFromSmiles(old.PROBE_SMILES))
    oldmap, oldopm = prior['probe_mapping'], prior['alignment']
    def inverse(xyz, t): return (xyz-np.asarray(t['translation_A'])) @ np.asarray(t['rotation']).T
    # old OPM -> old cluster -> historical QC receptor; then QC -> new cluster -> OPM.
    qc_xyz = inverse(inverse(ach.GetConformer().GetPositions(), oldopm), oldmap)
    prior_ca = prior['probe_reference']
    qc_path = Path(prior_ca)
    if not qc_path.is_file():
        alternatives = [ROOT/'inputs/probe/7TRS.pdb', Path('C:/projects/PACER_STAGE4_MD_backup')/'reference_assets/probe_only/minimized.pdb']
        qc_path = next((p for p in alternatives if p.is_file()), qc_path)
    assert qc_path.is_file(), 'BLOCKED: historical ACh QC receptor required for exact remapping'
    assert sha(qc_path) == prior['probe_reference_sha256'], 'Historical ACh QC reference hash mismatch'
    copy(qc_path, ROOT/'inputs/probe/historical_QC_probe_reference.pdb')
    cluster_ca = old.ca_rows(ROOT/'inputs/receptors/cluster_01_receptor.pdb',' ')
    matched = {r['mobile_resid'] for r in transform['correspondence']}
    probe_map = old.alignment(old.ca_rows(qc_path, prior['probe_reference_receptor_chain']),
                              [r for r in cluster_ca if r['resid'] in matched])
    xyz = old.transform_xyz(old.transform_xyz(qc_xyz, probe_map), transform)
    for i, p in enumerate(xyz): ach.GetConformer().SetAtomPosition(i, p)
    with Chem.SDWriter(str(ROOT/'inputs/probe/ACH_cluster_01_OPM.sdf')) as w: w.write(ach)
    dump(ROOT/'manifests/orientation_pacer_stage4_rerun_v01.json', dict(cluster_to_OPM=transform,
        QC_to_cluster=probe_map, historical_inverse_chain=[oldopm,oldmap],
        convention='row vectors: xyz_OPM = xyz_cluster @ rotation + translation_A',
        native_ACH_removed_atoms=sum(l.startswith(('ATOM','HETATM')) and l[17:20]=='ACH' for l in source),
        source_receptor_sha256=sha(ROOT/'inputs/receptors/cluster_01_receptor.pdb'),
        OPM_sha256=sha(ROOT/'inputs/receptors/7TRS_OPM.pdb'),
        probe_reference_sha256=sha(qc_path)))
    receptor_xyz = np.array([[float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for l in normalized
                            if l.startswith(('ATOM','HETATM')) and l[76:78].strip()!='H'])
    tree = cKDTree(receptor_xyz)
    probe_heavy=[a.GetIdx() for a in ach.GetAtoms() if a.GetAtomicNum()>1]
    probe_distances=tree.query(ach.GetConformer().GetPositions()[probe_heavy])[0]
    assert 1.0<min(probe_distances)<5.0, 'BLOCKED: ACh/protein severe clash or coordinate-frame mismatch'
    dump(ROOT/'manifests/probe_pacer_stage4_rerun_v01.json',dict(identity='CC(=O)OCC[N+](C)(C)C',
        formal_charge=1,atom_count=ach.GetNumAtoms(),sha256=sha(ROOT/'inputs/probe/ACH_cluster_01_OPM.sdf'),
        historical_pose_sha256=prior['transformed_probe_sha256'],min_protein_heavy_distance_A=float(min(probe_distances)),
        atom_mapping=[dict(historical_index=i+1,OPM_cluster_01_index=i+1,element=a.GetSymbol(),formal_charge=a.GetFormalCharge()) for i,a in enumerate(ach.GetAtoms())],
        contexts=['probe_only','candidate_probe'],native_receptor_ACh_retained=False,status='PASS'))
    audits=[]
    for r in records:
        cid=r['candidate_id']; sdf=ROOT/f'inputs/ligands/{cid}_glide_cluster_01.sdf'
        mol=Chem.SDMolSupplier(str(sdf), removeHs=False)[0]; assert mol is not None
        assert mol.GetNumAtoms()==len(r['atoms'])
        assert Chem.GetFormalCharge(mol)==r['formal_charge']
        raw_xyz=np.asarray([a['xyz'] for a in r['atoms']])
        assert np.max(np.abs(mol.GetConformer().GetPositions()-raw_xyz))<0.00011
        assert [a.GetAtomicNum() for a in mol.GetAtoms()]==[a['atomic_number'] for a in r['atoms']]
        assert [a.GetFormalCharge() for a in mol.GetAtoms()]==[a['formal_charge'] for a in r['atoms']]
        source_bonds={(min(b['a'],b['b']),max(b['a'],b['b'])):b['order'] for b in r['bonds']}
        sdf_bonds={(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx())+1,max(b.GetBeginAtomIdx(),b.GetEndAtomIdx())+1):b.GetBondTypeAsDouble() for b in mol.GetBonds()}
        assert set(source_bonds)==set(sdf_bonds), 'Source atom/bond mapping changed'
        # Maestro aromatic bond flags and SDF Kekule orders may differ in representation; sanitized chemical graph is verified below.
        state=Chem.MolFromSmiles(r['state_smiles']); parent=Chem.MolFromSmiles(r['canonical_smiles'])
        assert old.parent(state)==old.parent(parent), 'Parent connectivity/stereo mismatch'
        # Compare pose chemistry with prepared state; unspecified stereo in the state is not inferred from 3D.
        observed=Chem.RemoveHs(mol); template=Chem.RemoveHs(state)
        match=observed.GetSubstructMatch(template, useChirality=True)
        assert len(match)==template.GetNumAtoms()==observed.GetNumAtoms(), 'State/pose graph mismatch'
        assert sorted(b.GetBondTypeAsDouble() for b in observed.GetBonds())==sorted(b.GetBondTypeAsDouble() for b in template.GetBonds())
        assert Chem.GetFormalCharge(state)==Chem.GetFormalCharge(mol)
        # Preserve all source atomic indices explicitly, including hydrogens.
        for i,a in enumerate(mol.GetAtoms()): a.SetAtomMapNum(i+1)
        for i,p in enumerate(old.transform_xyz(raw_xyz,transform)): mol.GetConformer().SetAtomPosition(i,p)
        oriented_sdf=ROOT/f'inputs/ligands/{cid}_OPM.sdf'
        with Chem.SDWriter(str(oriented_sdf)) as w: w.write(mol)
        heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1]
        distances=tree.query(mol.GetConformer().GetPositions()[heavy])[0]
        assert np.min(distances)>1.0 and np.min(distances)<5.0, 'BLOCKED: ligand/receptor frame or severe overlap'
        achheavy=[a.GetIdx() for a in ach.GetAtoms() if a.GetAtomicNum()>1]
        pair=cKDTree(ach.GetConformer().GetPositions()[achheavy]).query(mol.GetConformer().GetPositions()[heavy])[0]
        assert min(pair)>1.0, 'BLOCKED: candidate/ACh severe overlap'
        safe={k:v for k,v in r.items() if k not in ['source_poseviewer','properties','atoms','bonds']}
        safe.update(source_poseviewer='inputs/poses/original_cluster_01_pv.maegz',
            raw_pose=f'inputs/poses/{cid}_glide_cluster_01.mae', raw_pose_sha256=sha(ROOT/f'inputs/poses/{cid}_glide_cluster_01.mae'),
            sdf=f'inputs/ligands/{cid}_glide_cluster_01.sdf',sdf_sha256=sha(sdf),
            oriented_sdf=f'inputs/ligands/{cid}_OPM.sdf',oriented_sdf_sha256=sha(oriented_sdf),
            atom_count=mol.GetNumAtoms(),atom_mapping=[dict(source_index=i+1,sdf_index=i+1,opm_index=i+1,atomic_number=a['atomic_number']) for i,a in enumerate(r['atoms'])],
            parent_microstate_identity='PASS', stereo_policy='Prepared Glide isomeric state authoritative; retain pose stereo, do not invent unspecified parent stereo',
            max_SDF_roundtrip_error_A=float(np.max(np.abs(Chem.SDMolSupplier(str(sdf),removeHs=False)[0].GetConformer().GetPositions()-raw_xyz))),
            min_candidate_protein_heavy_A=float(min(distances)),min_candidate_ACH_heavy_A=float(min(pair)), status='PASS')
        audits.append(safe)
    # Verify the returned docking receptor is in the same coordinate frame as the frozen ensemble source.
    glide_ca=old.ca_rows(ROOT/'inputs/receptors/glide_cluster_01.pdb',' ')
    frame=old.alignment(cluster_ca,glide_ca,expected_count=len(glide_ca))
    assert np.max(np.abs(np.array([r['xyz'] for r in cluster_ca])-np.array([r['xyz'] for r in glide_ca])))<0.2, 'Returned docking receptor differs from source baseline'
    assert np.max(np.abs(np.array(frame['rotation'])-np.eye(3)))<0.01 and np.linalg.norm(frame['translation_A'])<0.2, 'Poseviewer receptor frame mismatch'
    pose_audit=dict(status='PASS',candidate_poses=audits,poseviewer_frame_check=frame,
                    coordinates_generated_or_optimized=False,production_started=False)
    dump(OUT/'STAGE4_RERUN_POSE_AUDIT.json',pose_audit)
    systems=[]; jobs=[]
    for r in audits:
        for alias,context in [('A','apo'),('P','probe_only'),('C','candidate_no_probe'),('CP','candidate_probe')]:
            sid=f'pacer_stage4_rerun_v01__{r["candidate_id"]}__cluster_01__{context}'
            systems.append(dict(system_id=sid,candidate_id=r['candidate_id'],cascade_rank=r['cascade_rank'],cluster='cluster_01',
                context=context,alias=alias,has_candidate=alias in ['C','CP'],has_ACH=alias in ['P','CP'],
                receptor='inputs/receptors/oriented_cluster_01.pdb',receptor_sha256=sha(ROOT/'inputs/receptors/oriented_cluster_01.pdb'),
                candidate_sdf=r['oriented_sdf'] if alias in ['C','CP'] else '',
                probe_sdf='inputs/probe/ACH_cluster_01_OPM.sdf' if alias in ['P','CP'] else '',
                system_path=f'systems/{sid}',status='NOT_BUILT'))
            for rep,seed in enumerate([27101,38201,49301],1):
                job_id=f'{sid}__R{rep}'; output=f'production/{sid}/replica_{rep:02d}'
                jobs.append(dict(job_id=job_id,system_id=sid,candidate_id=r['candidate_id'],context=context,replica=rep,seed=seed,
                    production_ns=10,production_steps=5000000,timestep_fs=2,trajectory_ps=10,
                    output_path=output,checkpoint_path=output+'/segments/<end_step>/checkpoint.chk',status='AWAITING_SERVER_BUILD'))
    assert len(systems)==12 and len(jobs)==len({j['job_id'] for j in jobs})==len({j['output_path'] for j in jobs})==36
    writecsv(OUT/'STAGE4_RERUN_SYSTEM_MANIFEST.csv',systems); writecsv(OUT/'STAGE4_RERUN_JOB_MATRIX.csv',jobs)
    before=json.loads((OUT/'logs/protected_before.json').read_text()); after=snapshot()
    assert before==after, 'Protected historical directories changed during preparation'
    dump(OUT/'logs/protected_after.json',after)
    dependencies={m:importlib.util.find_spec(m) is not None for m in ['rdkit','openmm','openff','pdbfixer','openmmforcefields']}
    blockers=[]
    if not all(dependencies.values()): blockers.append('Local MD stack unavailable; complete parameterization, 12 membrane builds and topology/overlap checks on Linux')
    manifest=dict(protocol_id='pacer_stage4_rerun_v01',status='READY_FOR_SERVER_PREFLIGHT',ready_for_server_preflight=True,ready_for_production=False,
        frozen_ranking_sha256=RANK_HASH,candidates=audits,local_dependencies=dependencies,
        systems_built=0,systems=12,jobs=36,total_production_ns=360,production_started=False,
        blockers=blockers+['Both RTX 5090 devices / driver / OpenMM CUDA runtime and 12 system smokes unverified',
                          '36 matched equilibration/release QC receipts required before production'],
        protected_history_check='PASS: complete file size/mtime inventory unchanged; source input hashes verified',
        control_sharing='A/P physically shareable for cluster_01, but no trajectories shared; all 36 tasks execute separately',
        assets={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in (ROOT/'inputs').rglob('*') if p.is_file()})
    dump(OUT/'STAGE4_RERUN_INPUT_MANIFEST.json',manifest)
    for name in ['STAGE4_RERUN_INPUT_MANIFEST.json','STAGE4_RERUN_POSE_AUDIT.json','STAGE4_RERUN_SYSTEM_MANIFEST.csv','STAGE4_RERUN_JOB_MATRIX.csv']:
        copy(OUT/name,ROOT/'manifests'/name)
    report='''# PACER Stage4 Rerun v01 preflight

Status: READY_FOR_SERVER_PREFLIGHT; NOT READY_FOR_PRODUCTION. Production not started.

PASS: frozen ranking SHA256; fixed Cascade ranks 1/2/3; cluster_01 BEmin and scoring provenance; unique source MAEGZ entries; prepared-state/parent chemical identity; formal charges; full source atom mapping; original coordinate preservation; poseviewer receptor coordinate frame; OPM rigid transform and 270-CA correspondence; 274-residue capped receptor baseline; explicit historical ACh remapping and +1 charge; native ACh exclusion; 12 contexts and 36 unique paired-seed task paths; all packaged input files; historical directory inventory unchanged.

Initial candidate/protein and candidate/ACh heavy-atom overlap screen passed at 1 A severe-clash threshold. This is not full system overlap QC. Lipid/water/ion topology, charge sums, disulfides, force-field template completeness and finite energies MUST pass the Linux builder and CUDA system smoke before simulation eligibility.

Actual constructed membrane systems: 0/12. Windows preparation Python has RDKit but lacks OpenMM/OpenFF/PDBFixer/openmmforcefields. No parameterization or successful system build is claimed. AM1-BCC and all 12 builds remain server work.

All candidates use the same cluster_01 receptor. A/P controls are physically shareable, but this package retains 12 named systems and 36 separately executed tasks. Identical controls across candidate families must not be pooled as nine independent receptor conditions.

Scientific parameters follow the supplied historical implementations: ff14SB/lipid17/TIP3P/OpenFF 2.2.1, AM1-BCC, POPC, 0.15 M, 1 nm padding/cutoff, PME, 300 K, 1 bar, HBonds, 2 fs production. Smoke 500 steps/0.5 fs; NVT/NPT each 1 ps/0.5 fs; release 500/100/10/0 each 0.5 ps/1 fs. These are inherited short equilibration stages, not convergence evidence. Trajectory output is explicitly 10 ps, following the user and final handoff correction of historical 50 ps metadata.

Pending server gates: actual RTX 5090 driver/NVRTC compatibility on BOTH devices; dependency solve; AM1-BCC; 12 builds; full topology/disulfide/clash validation; CUDA 500-step system smokes; 36 per-replica equilibration/release audits. Any failure stops the workflow. No automatic scientific parameter changes.

Concrete charge note: the winning Glide states for PACERGEN00123 and PACERGEN00462 have formal charge +2; this is a valid microstate relation to the frozen neutral canonical parent. ACh is +1. The inherited builder neutralizes the membrane baseline before inserting ligands, so final systems can have nonzero net charge. The new builder records total system charge and ion counts; it does not silently add counterions or change the inherited protocol. Any proposed change to neutralization/ion conditions must be reported for scientific review before running it.

Commands: see STAGE4_RERUN_DEPLOYMENT_GUIDE.md. Production requires an explicit manual --authorize-production flag after all server gates pass.
'''
    (OUT/'STAGE4_RERUN_PREFLIGHT_REPORT.md').write_text(report,encoding='utf-8')
    guide='''# PACER Stage4 Rerun v01 deployment

Replace USER, HOST and REMOTE_PARENT with your actual SSH destination; none is assumed or configured by this preparation.

```powershell
scp C:/projects/GPCR-virtual-screening/project/results/pacer_stage4_rerun_v01/pacer_stage4_rerun_v01_server_bundle.tar.gz USER@HOST:REMOTE_PARENT/
```

On Linux, use a new directory. Never extract into an existing campaign.

```bash
mkdir pacer_stage4_rerun_v01_deployment
cd pacer_stage4_rerun_v01_deployment
tar -xzf ../pacer_stage4_rerun_v01_server_bundle.tar.gz
cd pacer_stage4_rerun_v01
sha256sum -c SHA256SUMS.txt
conda env create -n pacer-stage4-rerun-v01 -f environment-linux.yml
conda activate pacer-stage4-rerun-v01
python scripts/preflight_pacer_stage4_rerun_v01.py --environment-only
python scripts/build_pacer_stage4_rerun_v01.py --all
python scripts/preflight_pacer_stage4_rerun_v01.py --systems
python scripts/schedule_pacer_stage4_rerun_v01.py --phase prepare
```

The inherited environment file pins CUDA NVRTC 13.2. This is an installation recipe, not a claim that the server driver supports it. Both GPU contexts must execute finite-force integration before continuing. An incompatible solver/driver result is a blocker; do not silently change the scientific stack.

Only after reviewing all 36 preparation receipts, the operator may manually start production:

```bash
python scripts/schedule_pacer_stage4_rerun_v01.py --phase production --authorize-production
```

No production command was executed locally. Scheduler uses two workers, physical GPU indices 0 and 1, with CUDA_VISIBLE_DEVICES required unset, and an exclusive scheduler lock. Each task holds a separate lock. Re-run the same command to resume. A completed task is authenticated and skipped. Resume verifies job, config, system XML, PDB topology, starting state and checkpoint hashes. Production saves immutable 100 ps segments, each with its own DCD, CSV, checkpoint and committed receipt; an interrupted segment is preserved as an orphan and replayed from the last committed checkpoint, never appended to a potentially longer DCD. Concatenate only committed segments in numeric order. No XML fallback is silently used for incompatible binary checkpoints.

Outputs: systems/SYSTEM, smoke/SYSTEM, equilibration/SYSTEM/replica_NN, production/SYSTEM/replica_NN/segments/END_STEP. logs contains GPU preflight, scheduler/task logs and readiness receipts. All paths resolve from package root. No repository checkout or external private assets are required. AM1-BCC resources come from packaged conda dependencies; force-field asset files and hashes are recorded at build time.
'''
    (OUT/'STAGE4_RERUN_DEPLOYMENT_GUIDE.md').write_text(guide,encoding='utf-8')
    copy(OUT/'STAGE4_RERUN_DEPLOYMENT_GUIDE.md',ROOT/'README.md')
    copy(OUT/'STAGE4_RERUN_PREFLIGHT_REPORT.md',ROOT/'manifests/STAGE4_RERUN_PREFLIGHT_REPORT.md')
    # Windows provenance stays outside the portable package; extraction has already completed.
    for temp in ['extraction_requests.json','extracted_pose_records.json']:
        if (ROOT/'manifests'/temp).exists(): (ROOT/'manifests'/temp).unlink()
    # Discard only the interpreter cache generated by our local import; no bytecode in portable delivery.
    cache=ROOT/'scripts/legacy/__pycache__'
    if cache.exists():
        for generated in cache.glob('pacer_stage4_prospective.*.pyc'): generated.unlink()
        if not any(cache.iterdir()): cache.rmdir()
    checks={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt'}
    (ROOT/'SHA256SUMS.txt').write_text(''.join(f'{digest}  {name}\n' for name,digest in sorted(checks.items())),encoding='utf-8')
    archive=OUT/'pacer_stage4_rerun_v01_server_bundle.tar.gz'
    with tarfile.open(archive,'w:gz') as tar: tar.add(ROOT,arcname=ROOT.name)
    with tarfile.open(archive,'r:gz') as tar:
        for name,digest in checks.items():
            assert hashlib.sha256(tar.extractfile(ROOT.name+'/'+name).read()).hexdigest()==digest
    assert snapshot()==before and sha(STAGE3/'pacer_xr_rerun_full200_ranking.csv')==RANK_HASH
    dump(OUT/'logs/package_verification.json',dict(status='PASS',members_verified=len(checks),
        archive=str(archive),bytes=archive.stat().st_size,sha256=sha(archive),production_started=False))
    archive.with_suffix(archive.suffix+'.sha256').write_text(sha(archive)+'  '+archive.name+'\n',encoding='utf-8')
    print(json.dumps(dict(status=manifest['status'],systems_built=0,jobs=len(jobs),archive=str(archive),bytes=archive.stat().st_size,sha256=sha(archive)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['initialize','finalize']); a=p.parse_args()
    {'initialize':initialize,'finalize':finalize}[a.action]()
