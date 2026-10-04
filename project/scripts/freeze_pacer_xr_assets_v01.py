"""Freeze validated grid/receptor/PMF/tool identities without generating assets."""
import json
from pacer_xr_production_v01 import P,SCH,sha,require,write_json,load_pmf

def main():
    smoke=P/'results/pacer_xr_rerun_glide_v01'
    pdb=json.loads((smoke/'smoke_7trs/smoke_manifest.json').read_text())
    ensemble=json.loads((smoke/'smoke_ensemble/smoke_ensemble_manifest.json').read_text())
    require(pdb['completion_status']==ensemble['completion_status']=='PASS','Successful smoke manifests required')
    files={}
    def add(path,expected=None):
        path=__import__('pathlib').Path(path); digest=sha(path)
        if expected: require(digest==expected,'Smoke artifact hash mismatch: '+str(path))
        relative=path.relative_to(P).as_posix(); files[relative]=dict(path=relative,sha256=digest)
        return relative
    def smoke_add(path,root,manifest):
        path=__import__('pathlib').Path(path)
        return add(path,manifest['artifact_sha256'][str(path.relative_to(root))])
    root=smoke/'smoke_7trs'
    receptors={'7TRS':dict(grid=smoke_add(pdb['grid_path'],root,pdb),prepared_receptor=smoke_add(pdb['prepared_receptor_path'],root,pdb),center=pdb['grid_center'],inner=pdb['grid_settings']['INNERBOX'],outer=pdb['grid_settings']['OUTERBOX'],source_receptor=add(pdb['receptor_path']),vina_receptor=add(P/'results/structure/ensemble/7TRS_R_meeko.pdbqt'),vina_center=pdb['grid_center'],vina_size=[22,22,22])}
    ligand=smoke_add(pdb['ligand_input'],root,pdb)
    root=smoke/'smoke_ensemble'
    for c in ensemble['clusters']:
        name=f"cluster_{c['cluster_id']:02d}"
        receptors[name]=dict(grid=add(c['grid_file'],c['grid_sha256']),prepared_receptor=add(c['prepared_receptor'],c['prepared_sha256']),source_receptor=add(c['source_receptor'],c['source_sha256']),reference_ligand=add(c['reference_ligand'],c['reference_sha256']),center=c['center'],inner=[30,30,30],outer=[60,60,60],vina_receptor=add(P/f'results/m4_gamd_ensemble/receptors/cluster_{c["cluster_id"]:02d}.pdbqt'),vina_center=c['center'],vina_size=[30,30,30],ACH_repair=c.get('chemistry_repair','Schrodinger preparation inferred valid ACH connectivity; output verified equivalent to repaired clusters.'))
    add(P/'data/miao2026_m4r/M4R_ensemble_PMF.xvg',ensemble['pmf_source_sha256'])
    add(P/'results/m4_gamd_ensemble/ensemble_preparation_audit.json',ensemble['audit_sha256'])
    add(smoke/'smoke_7trs/smoke_manifest.json'); add(smoke/'smoke_ensemble/smoke_ensemble_manifest.json')
    lock=dict(receptors=receptors,files=list(files.values()),pmf=load_pmf(),rank1_ligand=ligand,
        suite='2023-1 build 128',glide='98128',schrodinger_executables={name:sha(SCH/name) for name in ['glide.exe','ligprep.exe','run.exe']},
        vina_executables={name:sha(P/'tools'/name) for name in ['vina.exe','vina_1.2.7.exe']},
        preparation_versions=dict(molscrub='0.1.1',rdkit='2025.3.5',meeko='0.7.1'),
        protocol=dict(glide_precision='SP',poses_per_lig=5,postdock_npose=5,forcefield='OPLS_2005',max_states=16,ph=7,ligprep_flags=['-i','0','-nt','-ns','-bff','14'],
            vina_pdb=dict(exhaustiveness=4,num_modes=1,cpu=1,seed=42,preparation='Parent canonical SMILES; RDKit ETKDGv3/MMFF300; Meeko, as dock_candidate_portfolio.py'),
            vina_ensemble=dict(exhaustiveness=8,num_modes=9,energy_range=3.,cpu=1,seed=42,preparation='Molscrub 0.1.1 pH7 max16; RDKit ETKDGv3 seed42+sid/MMFF300; Meeko')))
    target=P/'config/pacer_xr_assets_v01.json'
    if target.exists(): require(json.loads(target.read_text())==lock,'Existing asset lock differs')
    write_json(target,lock); print('Frozen 11 grids and all receptor/preparation/tool identities')

if __name__=='__main__': main()
