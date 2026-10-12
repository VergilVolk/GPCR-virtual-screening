#!/usr/bin/env python3
"""Independent post-hoc analysis. Never changes input or frozen ranking.

Usage: .venv-pacer-admet-v01/bin/python project/scripts/pacer_rerun_medchem_v01.py
Use --output NEW_DIRECTORY to reproduce; existing outputs are never overwritten.
"""
import argparse, hashlib, importlib.metadata as md, json, platform, subprocess, sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors as D, Lipinski as L, QED, rdMolDescriptors as R
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Contrib.SA_Score import sascorer

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT/'project/results/pacer_xr_rerun_v01/pacer_xr_rerun_full200_ranking.csv'
DEFAULT = ROOT/'project/results/pacer_rerun_admet_synthesis_v01'
TOP3 = ['PACERGEN02432','PACERGEN00123','PACERGEN00462']
KEY = ['candidate_id','canonical_smiles']
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def git(*args): return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()
def dump(out,name,obj):
    p=out/name
    with p.open('x') as f: json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False)
def csv(out,name,rows):
    df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    with (out/name).open('x') as f: df.to_csv(f,index=False)
    return df
def identity(row): return {k:row[k] for k in KEY}
def catalog(kind):
    p=FilterCatalogParams(); p.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs,kind)); return FilterCatalog(p)
def alerts(cat,m):
    records=[]
    for e in cat.GetMatches(m):
        matches=[]
        for fm in e.GetFilterMatches(m):
            atoms=sorted({int(p[1]) for p in fm.atomPairs})
            matches.append({'atom_indices_0based':atoms,'fragment_smiles':Chem.MolFragmentToSmiles(m,atoms) if atoms else None})
        records.append({'rule_id':e.GetDescription(),'catalog_properties':{k:e.GetProp(k) for k in e.GetPropList()},'matches':matches})
    return records
def encode(x): return json.dumps(x,ensure_ascii=False)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=DEFAULT);args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(INPUT,keep_default_na=False); mols=[Chem.MolFromSmiles(s) for s in df.canonical_smiles]
    canonical=[Chem.MolToSmiles(m,isomericSmiles=True) if m is not None else None for m in mols]
    discrepancies=[{'candidate_id':r.candidate_id,'input':r.canonical_smiles,'rdkit_canonical':s} for r,s in zip(df.itertuples(),canonical) if r.canonical_smiles!=s]
    validation={'input_file':str(INPUT),'sha256':sha(INPUT),'git_commit':git('rev-parse','HEAD'),'git_branch':git('branch','--show-current'),'parent_result_branch':'results/pacer-stage4-pbcfix01-fkg-v02-v01','initial_worktree_clean':True,'n_rows':len(df),'unique_ids':int(df.candidate_id.nunique()),'unique_input_smiles':int(df.canonical_smiles.nunique()),'unique_rdkit_canonical_smiles':len(set(canonical)),'all_parse':all(m is not None for m in mols),'id_smiles_bijection':df.candidate_id.nunique()==df.canonical_smiles.nunique()==len(df),'all_ids_PACERGEN':bool(df.candidate_id.str.fullmatch(r'PACERGEN\d+').all()),'top3_present':all(x in set(df.candidate_id) for x in TOP3),'input_smiles_column_matches_canonical':bool((df.smiles==df.canonical_smiles).all()),'canonical_serialization_differences':discrepancies,'identity_policy':'Original strings preserved verbatim; canonical serialization used only to check identity. No historical PACER mapping.'}
    dump(out,'INPUT_VALIDATION_v01.json',validation)
    assert len(df)==200 and validation['unique_ids']==validation['unique_rdkit_canonical_smiles']==200 and validation['all_parse'] and validation['top3_present'] and validation['all_ids_PACERGEN'], 'Identity validation failed'
    pkgs={}
    for pkg in ['rdkit','pandas','numpy','matplotlib','torch','chemprop','admet-ai','aizynthfinder','deepchem','openbabel']:
        try: pkgs[pkg]=md.version(pkg)
        except md.PackageNotFoundError: pkgs[pkg]=None
    dump(out,'ENVIRONMENT_v01.json',{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'packages':pkgs,'rdkit':rdBase.rdkitVersion,'default_python_before_setup':'3.14.7; RDKit/pandas/matplotlib/torch/ADMET/AiZynthFinder absent','environment_policy':'isolated .venv-pacer-admet-v01; original project environments untouched'})
    pains,brenk=catalog('PAINS'),catalog('BRENK');phys=[];drug=[];stereo=[];alertrows=[];syn=[]
    ion_queries={'amine_nonamide':'[NX3;!$(N-[C,S,P]=[O,S,N]);!$(N~[!#6;!#1])]', 'aromatic_nitrogen':'[nH0;+0]', 'carboxylic_acid':'C(=O)[O;H1,-1]', 'phenol':'[OX2H]c', 'sulfonamide':'S(=O)(=O)[NX3]'}
    for row,m in zip(df.to_dict('records'),mols):
        k=identity(row); mw=D.MolWt(m);lp=D.MolLogP(m);psa=D.TPSA(m);hbd=L.NumHDonors(m);hba=L.NumHAcceptors(m);rot=L.NumRotatableBonds(m);rings=R.CalcNumRings(m);mr=D.MolMR(m)
        carb=sum(a.GetAtomicNum()==6 for a in m.GetAtoms());het=sum(a.GetAtomicNum() not in [1,6] for a in m.GetAtoms());atoms_h=m.GetNumAtoms()+sum(a.GetTotalNumHs() for a in m.GetAtoms())
        phys.append({**k,'molecular_weight':mw,'exact_molecular_weight':D.ExactMolWt(m),'MolLogP':lp,'TPSA':psa,'HBD':hbd,'HBA':hba,'rotatable_bonds':rot,'aromatic_ring_count':R.CalcNumAromaticRings(m),'ring_count':rings,'heavy_atom_count':m.GetNumHeavyAtoms(),'FractionCSP3':R.CalcFractionCSP3(m),'formal_charge':Chem.GetFormalCharge(m),'molar_refractivity':mr,'atom_count_including_H':atoms_h,'carbon_count':carb,'heteroatom_count':het,'structure_valid':True,'sanitized_valence_valid':not bool(Chem.DetectChemistryProblems(m)),'component_count':len(Chem.GetMolFrags(m)),'evidence_type':'COMPUTED_DESCRIPTOR'})
        rules={'Lipinski':{'MW<=500':mw<=500,'MolLogP<=5':lp<=5,'HBD<=5':hbd<=5,'HBA<=10':hba<=10},'Veber':{'rotatable_bonds<=10':rot<=10,'TPSA<=140':psa<=140},'Ghose':{'160<=MW<=480':160<=mw<=480,'-0.4<=MolLogP<=5.6':-.4<=lp<=5.6,'40<=MR<=130':40<=mr<=130,'20<=atoms_including_H<=70':20<=atoms_h<=70},'Egan':{'MolLogP<=5.88':lp<=5.88,'TPSA<=131.6':psa<=131.6},'Muegge':{'200<=MW<=600':200<=mw<=600,'-2<=MolLogP<=5':-2<=lp<=5,'TPSA<=150':psa<=150,'ring_count<=7':rings<=7,'carbon_count>4':carb>4,'heteroatom_count>1':het>1,'rotatable_bonds<=15':rot<=15,'HBD<=5':hbd<=5,'HBA<=10':hba<=10}}
        dr={**k,'QED':QED.qed(m),'evidence_type':'EMPIRICAL_RULE_HEURISTIC','descriptor_convention':'RDKit MolLogP used consistently; Muegge uses MolLogP proxy rather than original XlogP; Egan original classifier approximated by common rectangular thresholds'}
        for name,checks in rules.items():
            fails=[cond for cond,passed in checks.items() if not passed];dr.update({name+'_status':'FAIL' if fails else 'PASS',name+'_violations':encode(fails),name+'_violation_count':len(fails)})
        dr['Lipinski_at_most_one_violation']=dr['Lipinski_violation_count']<=1
        drug.append(dr)
        centers=Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False)
        sinfo=list(Chem.FindPotentialStereo(m));details=[{'type':str(s.type),'centered_on':int(s.centeredOn),'specified':str(s.specified),'descriptor':str(s.descriptor)} for s in sinfo]
        nd=sum(str(s.type)=='Bond_Double' and str(s.specified)=='Unspecified' for s in sinfo)
        ns=sum(v=='?' for _,v in centers)
        # Enumerate on copies; unique=True handles symmetry. Cap is disclosed, never reported as exact if reached.
        isos=list(EnumerateStereoisomers(m,options=StereoEnumerationOptions(onlyUnassigned=True,unique=True,maxIsomers=4096,rand=0)))
        ion={name:[list(x) for x in m.GetSubstructMatches(Chem.MolFromSmarts(q))] for name,q in ion_queries.items()}
        stereo.append({**k,'specified_tetrahedral_centers':sum(v!='?' for _,v in centers),'unspecified_tetrahedral_centers':ns,'tetrahedral_atom_labels':encode(centers),'unspecified_double_bond_stereo':nd,'specified_double_bond_stereo':sum(str(s.type)=='Bond_Double' and str(s.specified)=='Specified' for s in sinfo),'potential_stereo_details':encode(details),'stereoisomer_count_preserving_specified':len(isos),'enumeration_cap':4096,'enumeration_capped':len(isos)>=4096,'independent_binary_upper_bound':2**(ns+nd),'ionizable_motif_matches_0based':encode(ion),'formal_charge':Chem.GetFormalCharge(m),'multicomponent':len(Chem.GetMolFrags(m))>1,'salt_form_assessment':'MULTICOMPONENT_REVIEW' if len(Chem.GetMolFrags(m))>1 else 'SINGLE_COMPONENT_NOT_A_SALT_PROOF','protonation_status':'NOT_ASSESSED','protonation_pH':None,'best_salt_form':'NOT_ASSESSED','interpretation':'Motif matches are structural only, not pKa or pH-specific speciation predictions; stereoisomer enumeration is theoretical, without conformational/experimental feasibility.'})
        pa,br=alerts(pains,m),alerts(brenk,m)
        alertrows.append({**k,'PAINS_hit':bool(pa),'PAINS_rule_count':len(pa),'PAINS_matches':encode(pa),'Brenk_hit':bool(br),'Brenk_rule_count':len(br),'Brenk_matches':encode(br),'status':'STRUCTURAL_ALERT' if pa or br else 'SCREENED_NO_MATCH','evidence_type':'STRUCTURAL_ALERT','atom_index_convention':'0-based input RDKit molecule, no atom reordering','interpretation':'No match is not proof of safety; Brenk is a broad medicinal-chemistry alert catalog, not a dedicated reactive-metabolite model.'})
        syn.append({**k,'SA_score':sascorer.calculateScore(m),'SA_model':'Ertl-Schuffenhauer RDKit Contrib SA_Score','SA_scale':'1 easier to 10 more difficult; heuristic, not route or price','SA_status':'COMPUTED_HEURISTIC','ring_count':rings,'unspecified_stereocenters':ns,'route_status':'NOT_ASSESSED','starting_material_availability':'AVAILABILITY_NOT_VERIFIED','yield_status':'NOT_ESTIMATED','yield':None,'cost_status':'NOT_ESTIMATED','cost':None,'currency':None,'price_date':None,'purification_included':None,'labor_included':None,'equipment_consumables_included':None,'yield_considered':None})
    csv(out,'ALL200_PHYSICOCHEMICAL_v01.csv',phys);csv(out,'ALL200_DRUGLIKENESS_v01.csv',drug);csv(out,'ALL200_STEREOCHEMISTRY_v01.csv',stereo);csv(out,'ALL200_STRUCTURAL_ALERTS_v01.csv',alertrows);csv(out,'PHASE_B_SYNTHESIS_HEURISTICS_v01.csv',syn)
    dump(out,'PHASE_AB_COMPLETE_v01.json',{'phase':'A,B','status':'COMPLETE','n':len(df),'input_sha256_after':sha(INPUT),'input_unchanged':sha(INPUT)==validation['sha256'],'script_sha256':sha(__file__)})
    print('PHASE A/B complete:',len(df),'molecules; canonical serialization differences:',len(discrepancies))
if __name__=='__main__': main()
