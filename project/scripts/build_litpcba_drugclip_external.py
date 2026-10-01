#!/usr/bin/env python3
"""Build a label-blind DrugCLIP external benchmark from the LIT-PCBA 9-target subset."""
from __future__ import annotations
import argparse, json, pickle, re
from pathlib import Path
import lmdb
import numpy as np
import pandas as pd
from rdkit import Chem


def blocks(path):
    text=path.read_text(encoding='utf-8',errors='replace');parts=text.split('@<TRIPOS>MOLECULE')
    for part in parts[1:]:yield '@<TRIPOS>MOLECULE'+part


def write(records,path):
    path.parent.mkdir(parents=True,exist_ok=True);env=lmdb.open(str(path),subdir=False,map_size=2**34)
    with env.begin(write=True) as txn:
        for i,record in enumerate(records):txn.put(str(i).encode(),pickle.dumps(record))
    env.sync();env.close()


def origin(conf):
    match=re.search(r'^origin\s*=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',conf.read_text(errors='replace'),re.M)
    if not match:raise ValueError(f'No GOLD origin in {conf}')
    return np.asarray(list(map(float,match.groups())),dtype=np.float32)


def pocket_record(folder,radius,max_atoms):
    center=origin(folder/'gold.conf');rows=[]
    for line in (folder/'xtal.pdb').read_text(errors='replace').splitlines():
        if not line.startswith('ATOM'):continue
        try:xyz=np.asarray([float(line[30:38]),float(line[38:46]),float(line[46:54])],dtype=np.float32)
        except ValueError:continue
        distance=float(np.linalg.norm(xyz-center))
        if distance<=radius:rows.append((distance,line[12:16].strip(),xyz))
    rows=sorted(rows,key=lambda row:row[0])[:max_atoms]
    if not rows:raise ValueError(f'Empty pocket for {folder.name}')
    return {'pocket':folder.name.rstrip('_'),'pocket_atoms':[row[1] for row in rows],
            'pocket_coordinates':[row[2] for row in rows]},center,len(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--radius',type=float,default=8.0);p.add_argument('--max-atoms',type=int,default=256);a=p.parse_args()
    molecule_records={};manifest=[];failures=[];pocket_audit={}
    for folder in sorted(path for path in a.root.iterdir() if path.is_dir()):
        target=folder.name.rstrip('_');seen_names=set()
        for block in blocks(folder/'ligdecs.mol2'):
            name=block.splitlines()[1].strip();seen_names.add(name);mol=Chem.MolFromMol2Block(block,sanitize=True,removeHs=False)
            if mol is None:
                failures.append({'target':target,'name':name,'reason':'rdkit_mol2_parse'});continue
            try:mol=Chem.RemoveHs(mol);Chem.SanitizeMol(mol);smiles=Chem.MolToSmiles(mol,canonical=True)
            except Exception:
                failures.append({'target':target,'name':name,'reason':'sanitize'});continue
            label=int(name.lower().startswith('active'));manifest.append({'pair_id':f'LITPCBA:{target}:{name}','target':target,'molecule_name':name,'canonical_smiles':smiles,'label':label})
            if smiles not in molecule_records:
                molecule_records[smiles]={'atoms':[atom.GetSymbol() for atom in mol.GetAtoms()],
                    'coordinates':[np.asarray(mol.GetConformer().GetPositions(),dtype=np.float32)],
                    'smi':smiles,'mol':mol,'label':smiles}
        record,center,n=pocket_record(folder,a.radius,a.max_atoms);path=a.output_dir/f'{target}.pocket.lmdb';write([record],path);pocket_audit[target]={'origin':center.tolist(),'atoms':n,'path':str(path),'source':'GOLD origin; protein ATOM records within fixed radius'}
    table=pd.DataFrame(manifest).drop_duplicates(['target','canonical_smiles','label']).reset_index(drop=True)
    conflicts=table.groupby(['target','canonical_smiles']).label.nunique();bad=set(conflicts[conflicts>1].index);table=table[~table.apply(lambda row:(row.target,row.canonical_smiles) in bad,axis=1)].copy()
    molecules=a.output_dir/'molecules.lmdb';write(list(molecule_records.values()),molecules);table.to_csv(a.output_dir/'pairs.csv',index=False)
    audit={'targets':sorted(table.target.unique()),'pairs':len(table),'unique_molecules':len(molecule_records),'actives_by_target':table.groupby('target').label.sum().astype(int).to_dict(),'rows_by_target':table.groupby('target').size().astype(int).to_dict(),'parse_failures':failures,'label_rule':'molecule name starts with active; copied before model inference','pockets':pocket_audit,'claim_boundary':'Independent non-GPCR active/inactive retrieval benchmark; not affinity or PAM efficacy.'}
    (a.output_dir/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
