import importlib.util
import pathlib
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

MODULE=pathlib.Path(__file__).with_name('audit_pacer_fkg_g2a_atom14_pbc.py')
spec=importlib.util.spec_from_file_location('g2a',MODULE)
g2a=importlib.util.module_from_spec(spec);spec.loader.exec_module(g2a)


def synthetic():
    seq='AGI'
    names=[];rid=[];res=['ALA','GLY','ILE']
    for r,aa in enumerate(seq):
        for nm in g2a.ATOM14[aa]:
            names.append('CD' if aa=='I' and nm=='CD1' else nm);rid.append(r)
    xyz=np.zeros((3,len(names),3),dtype=np.float32)
    # Synthetic coherent backbone, not a chemically complete force-field geometry.
    for r in range(3):
        for nm,offset in [('N',0.),('CA',1.4),('C',2.8),('O',3.7)]:
            s=[i for i,(rr,n) in enumerate(zip(rid,names)) if rr==r and n==nm][0]
            xyz[:,s,0]=r*3.8+offset
        for i,(rr,n) in enumerate(zip(rid,names)):
            if rr==r and n not in ('N','CA','C','O'):xyz[:,i,0]=r*3.8+1.5
    arr=np.zeros((3,3,14,3),dtype=np.float32)
    for r,aa in enumerate(seq):
        for s,nm in enumerate(g2a.ATOM14[aa]):
            target='CD' if aa=='I' and nm=='CD1' else nm
            idx=next(i for i,(rr,n) in enumerate(zip(rid,names)) if rr==r and n==target)
            arr[:,r,s]=xyz[:,idx]
    raw={'coords':xyz,'atom_names':np.array(names),'resnames':np.array(res),
         'resids':np.array([1,2,3]),'res_index':np.array(rid)}
    return arr,raw,seq

class G2ATests(unittest.TestCase):
    def test_reference_slot_counts(self):
        self.assertEqual(len(g2a.ATOM14),20)
        self.assertEqual(len(g2a.ATOM14['W']),14)
        self.assertEqual(g2a.ATOM14['I'][-1],'CD1')

    def test_mapping_and_ile_rename(self):
        arr,raw,seq=synthetic();m=g2a.validate_raw_mapping(arr,raw,seq)
        self.assertEqual(m['mapping_mismatch_slots_frames'],0)
        self.assertEqual(m['incorrect_missing_slot_frames'],0)
        self.assertFalse(m['valid_mask'][1,4])  # Gly has no CB.

    def test_swapped_slots_detected(self):
        arr,raw,seq=synthetic();arr[:,2,[1,2]]=arr[:,2,[2,1]]
        m=g2a.validate_raw_mapping(arr,raw,seq)
        self.assertGreater(m['mapping_mismatch_slots_frames'],0)

    def test_missing_mask_detected(self):
        arr,raw,seq=synthetic();arr[:,1,4,0]=1
        m=g2a.validate_raw_mapping(arr,raw,seq)
        self.assertGreater(m['incorrect_missing_slot_frames'],0)

    def test_pbc_broken_peptide_flag(self):
        arr,raw,seq=synthetic();raw['coords'][:,np.where(raw['res_index']==1)[0],0]+=25
        m=g2a.validate_raw_mapping(arr,raw,seq) # mapping must fail if atom14 not updated
        self.assertGreater(m['mapping_mismatch_slots_frames'],0)
        arr[:,1,:,:][...,0] +=25 # not used: geometry uses raw
        geom=g2a.geometry_screen(m)
        self.assertGreater(geom['PEPTIDE_C_N']['violations'],0)

    def test_sequence_mismatch_rejected(self):
        arr,raw,_=synthetic()
        with self.assertRaises(ValueError):g2a.validate_raw_mapping(arr,raw,'AAA')

if __name__=='__main__':unittest.main()

class ReferenceTests(unittest.TestCase):
    def test_official_reference_mismatch_is_detected(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/'residue_constants.py'
            p.write_text("restype_name_to_atom14_names={'ALA':['N','C','CA','O','CB']}\n",encoding='utf-8')
            result=g2a.inspect_reference_constants(p)
            # Missing other amino acids fails closed instead of saying MATCH.
            self.assertEqual(result['status'],'LOAD_FAILED')

    def test_official_reference_complete_match(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/'residue_constants.py'
            mapping={next(k for k,v in g2a.THREE.items() if v==one): names+['']*(14-len(names)) for one,names in g2a.ATOM14.items()}
            p.write_text('restype_name_to_atom14_names='+repr(mapping)+'\n',encoding='utf-8')
            result=g2a.inspect_reference_constants(p)
            self.assertEqual(result['status'],'MATCH')
