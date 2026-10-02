#!/usr/bin/env python
"""Build only: frozen GaMD conformation, inherited PACER-DC membrane protocol."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
import pacer_stage4_prospective as s4

def arguments():
    p=argparse.ArgumentParser(description=__doc__)
    choose=p.add_mutually_exclusive_group(required=True)
    choose.add_argument('--candidate',choices=list(s4.EXPECTED)); choose.add_argument('--all',action='store_true')
    p.add_argument('--dry-run',action='store_true',help='Validate identities and report asset readiness; no build or dynamics')
    p.add_argument('--prepare-only',action='store_true',help='Write rigid orientation/ligand inputs; no membrane or dynamics')
    p.add_argument('--opm-reference',type=Path,default=s4.PROJECT/'data/pdb/7TRS_OPM.pdb')
    p.add_argument('--opm-chain',default='R')
    p.add_argument('--probe-reference',type=Path,default=s4.PROJECT/'results/pacer_dc_openmm_reference_qc_v01/probe_only/minimized.pdb')
    p.add_argument('--probe-reference-chain',default='E')
    return p.parse_args()

def assets(args):
    return {'opm':args.opm_reference,'probe_reference':args.probe_reference,
            '7TRS':s4.PROJECT/'data/pdb/7TRS.pdb','ACH_template':s4.PROJECT/'data/pdb/m4_ligands/ACH_ideal.sdf'}

def orient(row,args,history,qc):
    from openmm import app,unit
    from openff.units import unit as off_unit
    source=Path(row['receptor_file']); c=int(row['source_cluster']); d=s4.OUT/'orientation'/f'cluster{c}'
    d.mkdir(parents=True,exist_ok=True)
    transform=s4.alignment(s4.ca_rows(source,' '),s4.ca_rows(args.opm_reference,args.opm_chain))
    s4.require(len(s4.ca_rows(source,' '))==274,'Frozen construct changed')
    lines=[]
    excluded=[]
    for line in source.read_text().splitlines():
        if line.startswith(('ATOM','HETATM')):
            if line[17:20]=='ACH':
                excluded.append(line); continue
            xyz=s4.transform_xyz([[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]],transform)[0]
            line=line[:30]+''.join(f'{x:8.3f}' for x in xyz)+line[54:]
        lines.append(line)
    oriented=d/'oriented_receptor.pdb'; oriented.write_text('\n'.join(lines)+'\n')
    candidate=qc.openff_molecule(s4.candidate_molecule(row),row['candidate_id'])
    candidate.conformers[0]=s4.transform_xyz(candidate.conformers[0].m_as(off_unit.angstrom),transform)*off_unit.angstrom
    raw_probe=qc.posed_from_pdb(s4.PROJECT/'data/pdb/7TRS.pdb',qc.LIGAND_DIR/'ACH_ideal.sdf',residue_name='ACH')
    s4.require(s4.graph(raw_probe)==s4.graph(Chem.MolFromSmiles(s4.PROBE_SMILES)),'ACh graph mismatch')
    probe=qc.openff_molecule(raw_probe,'ACH')
    probe=history.replace_conformer_from_chain(probe,args.probe_reference,'F')
    # QC reference receptor -> frozen cluster, then identical cluster -> OPM transform.
    probe_map=s4.alignment(s4.ca_rows(args.probe_reference,args.probe_reference_chain),
                          [r for r in s4.ca_rows(source,' ') if r['resid'] in {x['mobile_resid'] for x in transform['correspondence']}])
    probe.conformers[0]=s4.transform_xyz(s4.transform_xyz(probe.conformers[0].m_as(off_unit.angstrom),probe_map),transform)*off_unit.angstrom
    s4.require(int(probe.total_charge.m)==1,'ACh formal charge changed')
    for mol in [candidate,probe]:
        mol.generate_unique_atom_names()
        mol.to_file(str(d/(mol.name+'.sdf')),file_format='SDF')
    audit=dict(cluster=c,reference=str(args.opm_reference),reference_sha256=s4.sha(args.opm_reference),
        alignment_atom_selection='CA; full normalized-sequence correspondence',
        opm_reference_chain=args.opm_chain,mobile_chain=' ',
        mobile=str(source),mobile_sha256=s4.sha(source),alignment=transform,
        output=str(oriented),output_sha256=s4.sha(oriented),
        candidate_pose=row['pose_file'],candidate_pose_sha256=row['pose_file_sha256'],
        transformed_candidate_sha256=s4.sha(d/(candidate.name+'.sdf')),
        excluded_native_ach_atoms=len(excluded),native_ach_source_sha256=s4.sha(source),
        probe_policy='Remove native source ACH; add explicit historical 7TRS/QC ACh only in probe contexts',
        probe_reference=str(args.probe_reference),probe_reference_sha256=s4.sha(args.probe_reference),
        probe_reference_receptor_chain=args.probe_reference_chain,probe_reference_ligand_chain='F',
        probe_mapping=probe_map,transformed_probe_sha256=s4.sha(d/'ACH.sdf'),
        protocol_source_sha256=s4.sha(s4.PROJECT/'scripts/build_pacer_dc_membrane_reference.py'))
    s4.dump(d/'alignment_audit.json',audit)
    return oriented,candidate,probe,audit

def build(row,args):
    from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat,Platform,XmlSerializer,app,unit
    from pdbfixer import PDBFixer
    import build_pacer_dc_membrane_reference as history
    import build_pacer_dc_openmm_reference_systems as qc
    oriented,candidate,probe,alignment=orient(row,args,history,qc)
    if args.prepare_only: return
    # Normalize Amber residue naming on a copy only. Keep caps and all residues.
    normalized=oriented.with_name('normalized_receptor.pdb')
    lines=oriented.read_text().splitlines()
    normalized.write_text('\n'.join(s4.normalized_receptor_lines(lines))+'\n')
    fixer=PDBFixer(filename=str(normalized)); fixer.findMissingResidues(); fixer.missingResidues={}
    fixer.findMissingAtoms(); missing=sum(len(x) for x in fixer.missingAtoms.values())
    before={(a.residue.id,a.name):np.asarray(fixer.positions[a.index].value_in_unit(unit.nanometer))
            for a in fixer.topology.atoms() if a.element.symbol!='H'}
    fixer.addMissingAtoms()
    for a in fixer.topology.atoms():
        key=(a.residue.id,a.name)
        if key in before: s4.require(np.allclose(before[key],fixer.positions[a.index].value_in_unit(unit.nanometer),atol=1e-8),'Existing heavy atom moved during repair')
    protein=app.Modeller(fixer.topology,fixer.positions)
    protein.topology.createDisulfideBonds(protein.positions)
    cyx={line[22:26].strip() for line in oriented.read_text().splitlines()
         if line.startswith('ATOM') and line[17:20]=='CYX'}
    paired={atom.residue.id for a,b in protein.topology.bonds()
            if a.name==b.name=='SG' for atom in (a,b)}
    s4.require(cyx.issubset(paired),'Frozen CYX disulfide connectivity not preserved')
    s4.require(sum(r.name not in ['ACE','NME'] for r in protein.topology.residues())==274,
               'Receptor residue count changed during repair')
    base,hqc,attempts=history.rebuild_protein_hydrogens(protein,protein.positions,history.forcefield())
    after_heavy={(a.residue.id,a.name):np.asarray(base.positions[a.index].value_in_unit(unit.nanometer))
                 for a in base.topology.atoms() if a.element.symbol!='H'}
    for key,xyz in before.items():
        s4.require(key in after_heavy and np.allclose(xyz,after_heavy[key],atol=1e-8),
                   'Hydrogen preparation moved existing protein heavy atoms')
    protein_atoms=base.topology.getNumAtoms()
    base.topology.createDisulfideBonds(base.positions)
    base.addMembrane(history.forcefield(),lipidType='POPC',minimumPadding=1*unit.nanometer,
                     ionicStrength=0.15*unit.molar,neutralize=True)
    shared=s4.OUT/'shared_bases'/f"cluster{row['source_cluster']}"; shared.mkdir(parents=True,exist_ok=True)
    base_path=shared/'shared_OPM_POPC_water.pdb'
    with base_path.open('w') as f: app.PDBFile.writeFile(base.topology,base.positions,f,keepIds=True)
    for molecule in [candidate,probe]:
        before_xyz=molecule.conformers[0].m.copy()
        molecule.assign_partial_charges('am1bcc')
        s4.require(np.array_equal(before_xyz,molecule.conformers[0].m),'Parameterization changed pose coordinates')
    defs=s4.systems_for(row)
    for name,names in defs.items():
        output=s4.OUT/'systems'/name
        s4.require(not output.exists(),'Refusing to overwrite system: '+str(output))
        output.mkdir(parents=True)
        ligands=[copy.deepcopy(probe if n=='ACH' else candidate) for n in names]
        modeller=app.Modeller(base.topology,base.positions)
        for mol in ligands: modeller.add(mol.to_topology().to_openmm(),mol.conformers[0].to_openmm())
        removed=history.remove_ligand_solvent_clashes(modeller,sum(m.n_atoms for m in ligands))
        system=history.forcefield(ligands).createSystem(modeller.topology,nonbondedMethod=app.PME,
            nonbondedCutoff=1*unit.nanometer,constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=5e-4)
        history.restrain_protein(system,modeller.topology,modeller.positions,protein_atoms,k=1000)
        system.addForce(MonteCarloMembraneBarostat(1*unit.bar,0*unit.bar*unit.nanometer,300*unit.kelvin,
            MonteCarloMembraneBarostat.XYIsotropic,MonteCarloMembraneBarostat.ZFree,25))
        integrator=LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,0.002*unit.picoseconds)
        sim=app.Simulation(modeller.topology,system,integrator,Platform.getPlatformByName('CPU'))
        sim.context.setPositions(modeller.positions)
        initial=sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        s4.require(np.isfinite(initial),'Nonfinite initial energy')
        sim.minimizeEnergy(maxIterations=100)
        final=sim.context.getState(getEnergy=True,getPositions=True)
        energy=final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        s4.require(np.isfinite(energy) and np.isfinite(final.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all(),'Nonfinite minimized system')
        for filename,positions in [('input.pdb',modeller.positions),('minimized.pdb',final.getPositions())]:
            with (output/filename).open('w') as f: app.PDBFile.writeFile(modeller.topology,positions,f)
        (output/'system.xml').write_text(XmlSerializer.serialize(system))
        (output/'state.xml').write_text(XmlSerializer.serialize(final))
        charges={m.name:dict(forcefield='openff_unconstrained-2.2.1.offxml',method='AM1-BCC',
            formal_charge=int(m.total_charge.m),partial_charges=m.partial_charges.m.tolist(),
            atom_names=[a.name for a in m.atoms],n_atoms=m.n_atoms) for m in ligands}
        s4.dump(output/'parameters.json',charges)
        s4.require(all(abs(sum(v['partial_charges'])-v['formal_charge'])<1e-5 for v in charges.values()),
                   'Parameter charge sum mismatch')
        audit=dict(status='COMPLETE',candidate=row['candidate_id'],cluster=int(row['source_cluster']),
            context=name.split('__')[-1],canonical_smiles=row['canonical_smiles'],state_smiles=row['source_state_smiles'],
            formal_state_charge=Chem.GetFormalCharge(Chem.MolFromSmiles(row['source_state_smiles'])),
            source_assets={key:dict(path=str(path),sha256=s4.sha(path)) for key,path in assets(args).items()},
            alignment=alignment,shared_base_sha256=s4.sha(base_path),ligands=names,
            protein_atoms=protein_atoms,candidate_atoms=candidate.n_atoms if candidate.name in names else 0,
            probe_atoms=probe.n_atoms if 'ACH' in names else 0,total_particles=system.getNumParticles(),
            missing_protein_atoms_added=missing,hydrogen_qc=hqc,hydrogen_attempts=attempts,
            cap_topology_policy='Existing ACE/NME caps retained; TER placed after NME, preserving both receptor fragments',
            removed_clashing_solvent_residues=removed,initial_energy_kj_mol=initial,
            minimized_energy_kj_mol=energy,finite=True,
            outputs={f:s4.sha(output/f) for f in ['input.pdb','minimized.pdb','system.xml','state.xml','parameters.json']},
            production_started=False)
        s4.dump(output/'build_audit.json',audit)
        (output/'build.log').write_text(json.dumps(dict(status='COMPLETE',initial_energy_kj_mol=initial,
                                                       minimized_energy_kj_mol=energy,finite=True))+'\n')
        print(json.dumps(dict(system=name,status='COMPLETE')),flush=True)
        del sim,integrator

def main():
    args=arguments(); rows=s4.frozen_rows(); chosen=rows if args.all else [r for r in rows if r['candidate_id']==args.candidate]
    paths=assets(args); missing=[str(p) for p in paths.values() if not p.is_file()]
    report=dict(stage3d_input_validation='PASS',candidates=[r['candidate_id'] for r in rows],
        selected=[r['candidate_id'] for r in chosen],system_count=12,missing_assets=missing,
        assets={name:dict(path=str(path),exists=path.is_file(),sha256=s4.sha(path) if path.is_file() else None) for name,path in paths.items()},
        production_started=False,dry_run=args.dry_run,
        identity=dict(parent_graph='PASS',microstate_charge='PASS',heavy_atom_mapping='PASS',
                      stereo='Full source SMILES identity; unspecified stereo is not inferred from coordinates'))
    if args.dry_run:
        s4.update_manifest(rows); s4.dump(s4.OUT/'builder_input_validation.json',report)
        print(json.dumps(report,indent=2)); return
    s4.require(not missing,'Missing formal assets: '+', '.join(missing))
    # Complete clusters can be reused without rebuilding their membrane base.
    pending=[]
    for row in chosen:
        names=s4.systems_for(row)
        existing=[name for name in names if (s4.OUT/'systems'/name).exists()]
        if existing:
            s4.require(len(existing)==4,'Partial cluster build exists; refusing to mix shared membrane bases')
            for name in existing:
                directory=s4.OUT/'systems'/name
                audit=json.loads((directory/'build_audit.json').read_text())
                s4.require(audit['status']=='COMPLETE','Existing build not complete')
                s4.require(audit['alignment']['mobile_sha256']==s4.sha(row['receptor_file']) and
                           audit['alignment']['candidate_pose_sha256']==row['pose_file_sha256'],'Existing build identity mismatch')
                for key,path in paths.items(): s4.require(audit['source_assets'][key]['sha256']==s4.sha(path),'Existing reference changed')
                for file,digest in audit['outputs'].items(): s4.require(s4.sha(directory/file)==digest,'Existing built asset changed')
            print(json.dumps(dict(candidate=row['candidate_id'],status='REUSED_COMPLETE_CLUSTER')),flush=True)
        else: pending.append(row)
    for row in pending:
        build(row,args); s4.update_manifest(rows)
    s4.update_manifest(rows)

if __name__=='__main__': main()
