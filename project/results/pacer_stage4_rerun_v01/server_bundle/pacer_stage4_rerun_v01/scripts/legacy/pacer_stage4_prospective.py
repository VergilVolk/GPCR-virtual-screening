"""Identity and geometry helpers for the prospective Stage4 preparation only."""
import csv
import hashlib
import json
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

PROJECT = Path(__file__).resolve().parents[1]
REPO = PROJECT.parent
OUT = PROJECT/'results/pacer_stage4_prospective_10ns_v01'
FROZEN = PROJECT/'results/pacer_stage3d_corrected_20261002_v01'
FREEZE_HASHES = {
 'stage3d_selected_candidates.csv':'5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3',
 'md_shortlist_pose_manifest.csv':'8ef3c4be8ff3800543d09e3b37c11331671e961c7bd3f29cc70850b787e311a6'}
EXPECTED = {'PACER0010':(9,'836a68d0cced23dc6e345a55d43e66119c27f26b9d49803635f97aed56e1e9de'),
 'PACER0073':(0,'3c503fc44a107cca2612dc71d9b0d529577bb5447eb2ea329dbe3adc32bcc37a'),
 'PACER0027':(4,'88cd9f7b228a25f65e82c9d82a49f65653cd67eb2abe53c3280146cc001d737e')}
SEEDS = [27101,38201,49301]
PROBE_SMILES = 'CC(=O)OCC[N+](C)(C)C'
NORMAL = {'CYX':'CYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

def require(condition, message):
    if not condition: raise ValueError(message)

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): digest.update(chunk)
    return digest.hexdigest()

def dump(path, data):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    tmp.replace(path)

def graph(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol),canonical=True,isomericSmiles=True)

def parent(mol):
    return graph(rdMolStandardize.Uncharger().uncharge(rdMolStandardize.Cleanup(Chem.RemoveHs(mol))))

def candidate_molecule(row):
    """Recover MODEL 1 with full Meeko SMILES atom mapping, never pose search."""
    lines=(REPO/row['pose_file']).read_text().split('ENDMDL',1)[0].splitlines()
    require(lines[0].strip()=='MODEL 1','Expected frozen MODEL 1')
    smiles=next(x[14:] for x in lines if x.startswith('REMARK SMILES ') and not x.startswith('REMARK SMILES IDX '))
    mol=Chem.MolFromSmiles(smiles); state=Chem.MolFromSmiles(row['source_state_smiles'])
    canonical=Chem.MolFromSmiles(row['canonical_smiles'])
    require(all(x is not None for x in [mol,state,canonical]),'Invalid frozen graph')
    require(graph(mol)==graph(state),'Pose/state isomeric graph mismatch')
    require(parent(mol)==parent(canonical),'State/parent graph mismatch')
    require(Chem.GetFormalCharge(mol)==Chem.GetFormalCharge(state),'State charge mismatch')
    coords={int(x[6:11]):[float(x[a:b]) for a,b in [(30,38),(38,46),(46,54)]]
            for x in lines if x.startswith(('ATOM','HETATM'))}
    atom_types={int(x[6:11]):x.split()[-1] for x in lines if x.startswith(('ATOM','HETATM'))}
    pairs=[]; hs=[]
    for x in lines:
        for label,target in [('REMARK SMILES IDX ',pairs),('REMARK H PARENT ',hs)]:
            if x.startswith(label):
                nums=list(map(int,x[len(label):].split()))
                require(len(nums)%2==0,'Odd atom-map entry')
                target.extend(zip(nums[::2],nums[1::2]))
    require(len(pairs)==mol.GetNumAtoms()==canonical.GetNumHeavyAtoms(),'Heavy-atom count mismatch')
    require({a for a,b in pairs}==set(range(1,mol.GetNumAtoms()+1)),'Incomplete SMILES mapping')
    require(len({b for a,b in pairs})==len(pairs),'Non-bijective atom mapping')
    elements={'A':'C','C':'C','N':'N','NA':'N','NS':'N','O':'O','OA':'O','OS':'O',
              'S':'S','SA':'S','P':'P','F':'F','Cl':'Cl','Br':'Br','I':'I','H':'H','HD':'H','HS':'H'}
    require({b for a,b in pairs}=={serial for serial,t in atom_types.items() if elements.get(t)!='H'},
            'Pose heavy-atom count/mapping mismatch')
    for a,b in pairs: require(elements.get(atom_types[b])==mol.GetAtomWithIdx(a-1).GetSymbol(),'Pose atom element mismatch')
    conf=Chem.Conformer(mol.GetNumAtoms())
    for a,b in pairs: conf.SetAtomPosition(a-1,coords[b])
    mol.AddConformer(conf)
    mol=Chem.AddHs(mol,addCoords=True); consumed=set()
    for a,b in hs:
        options=[n.GetIdx() for n in mol.GetAtomWithIdx(a-1).GetNeighbors()
                 if n.GetAtomicNum()==1 and n.GetIdx() not in consumed]
        require(bool(options),'Invalid H-parent mapping')
        mol.GetConformer().SetAtomPosition(options[0],coords[b]); consumed.add(options[0])
    for a,b in pairs:
        require(np.array_equal(np.asarray(mol.GetConformer().GetAtomPosition(a-1)),coords[b]),'Frozen coordinate changed')
    return mol

def frozen_rows():
    for name,digest in FREEZE_HASHES.items(): require(sha(FROZEN/name)==digest,'Freeze hash mismatch: '+name)
    def read(name):
        with (FROZEN/name).open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
    rows=read('md_shortlist_pose_manifest.csv'); selected=read('stage3d_selected_candidates.csv')
    require(len(rows)==3 and {r['candidate_id'] for r in rows}==set(EXPECTED),'Candidate set mismatch')
    require(len(selected)==3 and {r['candidate_id'] for r in selected}==set(EXPECTED),'Selection mismatch')
    known={r['candidate_id']:r for r in selected}
    manifest=json.loads((OUT/'master_manifest.json').read_text(encoding='utf-8'))
    historical={r['candidate_id']:r for r in manifest['candidates']}
    for row in rows:
        cid=row['candidate_id']; cluster,digest=EXPECTED[cid]
        require(int(row['source_cluster'])==cluster and int(row['source_state'])==0,'Cluster/state mismatch')
        require(row['pose_file_sha256']==digest==sha(REPO/row['pose_file']),'Pose hash mismatch')
        require(row['canonical_smiles']==known[cid]['canonical_smiles'],'Selection graph mismatch')
        receptor=PROJECT/f'results/m4_gamd_ensemble/receptors/cluster_{cluster:02d}_receptor.pdb'
        require(sha(receptor)==historical[cid]['receptor_sha256'],'Frozen receptor hash mismatch')
        row['receptor_file']=str(receptor); candidate_molecule(row)
    return rows

def ca_rows(path, chain):
    result=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('ATOM') and line[21]==chain and line[12:16].strip()=='CA' and line[16] in (' ','A'):
            result.append(dict(residue=line[17:20],resid=line[22:26].strip(),chain=chain,
                xyz=[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]))
    require(bool(result),'No alignment CA atoms: '+str(path))
    require(len({(x['resid']) for x in result})==len(result),'Duplicate CA residue IDs')
    return result

def alignment(mobile, reference, expected_count=270):
    """Require complete reference sequence correspondence; normalize names only."""
    norm=lambda r:NORMAL.get(r['residue'],r['residue'])
    blocks=SequenceMatcher(a=[norm(r) for r in mobile],b=[norm(r) for r in reference],autojunk=False).get_matching_blocks()
    pairs=[(i+k,j+k) for i,j,n in blocks for k in range(n)]
    require(len(reference)==len(pairs)==expected_count,'Incomplete reference CA correspondence')
    p=np.array([mobile[i]['xyz'] for i,j in pairs]); q=np.array([reference[j]['xyz'] for i,j in pairs])
    pc=p.mean(0); qc=q.mean(0); u,_,vt=np.linalg.svd((p-pc).T@(q-qc)); rotation=u@vt
    if np.linalg.det(rotation)<0: u[:,-1]*=-1; rotation=u@vt
    translation=qc-pc@rotation
    require(np.allclose(rotation.T@rotation,np.eye(3)) and np.linalg.det(rotation)>0,'Invalid rigid transform')
    return dict(ca_count=len(pairs),rmsd_A=float(np.sqrt(np.mean(np.sum((p@rotation+translation-q)**2,axis=1)))),
        rotation=rotation.tolist(),translation_A=translation.tolist(),source_center_A=pc.tolist(),target_center_A=qc.tolist(),
        correspondence=[dict(mobile_resid=mobile[i]['resid'],reference_resid=reference[j]['resid'],
            mobile_residue=mobile[i]['residue'],reference_residue=reference[j]['residue']) for i,j in pairs])

def transform_xyz(xyz, audit):
    return np.asarray(xyz)@np.array(audit['rotation'])+np.array(audit['translation_A'])

def normalized_receptor_lines(lines):
    """Keep capped fragments bonded: source TER precedes, rather than follows, NME.

    Move only topology separators after the already present NME caps. No new
    residues/bonds across the missing loop, and no coordinates are changed.
    """
    groups=[]
    for line in lines:
        if line.startswith(('ATOM','HETATM')):
            # Docking receptor files retain the native co-simulated ACh. The
            # four-context MD protein base excludes it; probe is added explicitly.
            if line[17:20]=='ACH': continue
            key=(line[21:27],line[17:20])
            if not groups or groups[-1][0]!=key: groups.append((key,[]))
            groups[-1][1].append(line)
    result=[]
    for (_,name),group in groups:
        result.extend(line[:17]+NORMAL.get(name,name)+line[20:] for line in group)
        if name=='NME': result.append('TER')
    require(sum(key[1]=='ACE' for key,group in groups)==2 and
            sum(key[1]=='NME' for key,group in groups)==2,'Expected two ACE/NME-capped receptor fragments')
    result.append('END')
    return result

def systems_for(row):
    cid=row['candidate_id']; c=int(row['source_cluster'])
    return {f'cluster{c}__apo':[],f'cluster{c}__probe_only':['ACH'],
            f'{cid}__candidate_no_probe':[cid],f'{cid}__candidate_probe':['ACH',cid]}

def update_manifest(rows):
    master_path=OUT/'master_manifest.json'; m=json.loads(master_path.read_text(encoding='utf-8'))
    m.update(unique_systems=12,expected_job_count=36,replicas=3,production_launch_status='NOT_STARTED',
             ready_to_launch=False,production_started=False)
    systems=[]; jobs=[]
    for row in rows:
        for name,ligands in systems_for(row).items():
            context=name.split('__')[-1]; audit=OUT/'systems'/name/'build_audit.json'
            build=json.loads(audit.read_text()) if audit.exists() else {}
            systems.append(dict(system=name,candidate=row['candidate_id'],cluster=int(row['source_cluster']),
                context=context,ligands=ligands,build_status=build.get('status','NOT_BUILT'),
                build_audit_sha256=sha(audit) if audit.exists() else None))
            for replica,seed in enumerate(SEEDS,1):
                output=f'production/{name}/replica_{replica:02d}'
                jobs.append(dict(candidate_id=row['candidate_id'],receptor_cluster=row['source_cluster'],context=context,
                    system=name,replica=replica,seed=seed,production_length_ns=10,production_steps=5000000,
                    pose_sha256=row['pose_file_sha256'] if row['candidate_id'] in ligands else '',
                    system_path=f'systems/{name}',output_path=output,checkpoint_path=output+'/checkpoint.chk',
                    status='AWAITING_EQUILIBRATION' if build.get('status')=='COMPLETE' else 'NOT_BUILT'))
    require(not any((OUT/'production').rglob('progress.json')),'Do not rewrite job matrix after production starts')
    with (OUT/'job_matrix.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(jobs[0])); writer.writeheader(); writer.writerows(jobs)
    m['systems']=systems; m['job_matrix_sha256']=sha(OUT/'job_matrix.csv')
    m['system_build_status']='COMPLETE' if all(x['build_status']=='COMPLETE' for x in systems) else 'PARTIAL'
    m['equilibration_status']='NOT_RUN'; m['preflight_status']='NOT_RUN'
    smoke=[]
    for item in systems:
        smoke_path=OUT/'smoke'/item['system']/'audit.json'
        result=json.loads(smoke_path.read_text()) if smoke_path.exists() else {}
        passed=result.get('finite',False) and result.get('steps')==500 and result.get('platform')=='CUDA'
        passed=passed and result.get('build_audit_sha256')==item['build_audit_sha256'] and item['build_status']=='COMPLETE'
        smoke.append(dict(system=item['system'],status='PASS' if passed else 'NOT_RUN'))
    m['smoke_systems']=smoke
    m['cuda_smoke_status']='PASS' if all(x['status']=='PASS' for x in smoke) else 'NOT_RUN'
    m['blockers']=['Equilibration/release and production are outside this implementation task',
                   'Actual membrane construction and CUDA smoke require formal OPM/ACh assets and server execution']
    m['implementation_sources']={name:sha(PROJECT/'scripts'/name) for name in
         ['pacer_stage4_prospective.py','build_pacer_stage4_prospective_membrane.py','run_pacer_stage4_membrane_smoke.py']}
    dump(master_path,m)
