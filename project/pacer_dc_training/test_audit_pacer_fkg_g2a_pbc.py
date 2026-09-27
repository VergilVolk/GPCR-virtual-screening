import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

from audit_pacer_fkg_g2a_pbc import inputs, check_metadata, valid_box, compare_raw, analyze_frame


class FakeResidues:
    resnames=np.asarray(['ALA','GLY'])
    resids=np.asarray([1,2])
    def __iter__(self):return iter([SimpleNamespace(resindex=10),SimpleNamespace(resindex=12)])

class FakeSelection:
    names=np.asarray(['N','CA','C','N','CA','C'])
    resindices=np.asarray([10,10,10,12,12,12])
    residues=FakeResidues()
    def __len__(self):return len(self.names)


class G2APBCTests(unittest.TestCase):
    def test_paths_exactly_40_and_special_apo(self):
        r=inputs(Path('S'),Path('R'))
        self.assertEqual(len(r),40)
        self.assertEqual(len({str(x['dcd']) for x in r}),8)
        self.assertEqual(len({str(x['pdb']) for x in r}),4)
        special=[x for x in r if (x['system'],x['replica'])==('apo',3)]
        self.assertEqual(len(special),5)
        self.assertTrue(all(x['dcd'].name=='trajectory_corrected.dcd' for x in special))

    def _metadata(self, f, src='project/results/pacer_dc_recovery_merged_v01/apo/replica_03/trajectory_corrected.dcd'):
        f.parent.mkdir(parents=True,exist_ok=True)
        f.write_text(json.dumps({'rows':[{'context':'apo','status':'ok','frame_count':100,'start_frame':0,'end_frame':99,'trajectory_path':src}]}),encoding='utf8')

    def test_metadata_corrected_provenance(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'t.json';self._metadata(p)
            rec={'manifest':p,'context':'apo','window':0,'system':'apo','replica':3,'raw':Path('a')}
            self.assertEqual(check_metadata(rec),(0,99))

    def test_metadata_rejects_plain_r3_apo(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'t.json';self._metadata(p,'project/results/pacer_dc_production_v01/apo/replica_03/trajectory.dcd')
            rec={'manifest':p,'context':'apo','window':0,'system':'apo','replica':3,'raw':Path('a')}
            with self.assertRaisesRegex(ValueError,'wrong trajectory provenance'):check_metadata(rec)

    def _raw(self,path,coords,mask=None):
        np.savez(path,coords=coords,atom_names=np.asarray(['N','CA','C','N','CA','C']),resnames=np.asarray(['ALA','GLY']),resids=np.asarray([1,2]),res_index=np.asarray([0,0,0,1,1,1]) if mask is None else mask)

    def test_raw_match_exact(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.npz'; coords=np.ones((100,6,3),dtype=np.float32);self._raw(p,coords)
            self.assertEqual(compare_raw(p,FakeSelection(),coords,2,0,'apo',0.05),0)

    def test_raw_rejects_coordinate_shift(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.npz';coords=np.ones((100,6,3),dtype=np.float32);self._raw(p,coords)
            moved=coords.copy();moved[7,1,0]+=0.2
            with self.assertRaisesRegex(ValueError,'coordinate difference'):compare_raw(p,FakeSelection(),moved,2,0,'apo',0.05)

    def test_raw_rejects_residue_mapping(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.npz';coords=np.ones((100,6,3),dtype=np.float32);self._raw(p,coords,mask=np.asarray([0,0,0,1,1,0]))
            with self.assertRaisesRegex(ValueError,'residue atom mapping'):compare_raw(p,FakeSelection(),coords,2,0,'apo',0.05)

    def test_box_check(self):
        self.assertTrue(valid_box([80,80,100,90,90,90]))
        self.assertFalse(valid_box([0,80,100,90,90,90]))
        self.assertFalse(valid_box(None))
        self.assertFalse(valid_box([80,80,100,0,90,90]))

    def test_pbc_split_synthetic(self):
        # First pair direct 98 A, minimum image 2 A (mock MIC to avoid dependency in unit test).
        ca=np.asarray([[1.,0,0],[99.,0,0],[102.8,0,0]])
        c=dict(valid_box_frames=0,invalid_box_frames=0,ca_pbc_split_pairs=0,ca_non_pbc_severe_pairs=0,temporal_pbc_jump_atoms=0,max_ca_direct_A=0.)
        x={k:[] for k in ('invalid_box','ca_pbc_split','ca_non_pbc_severe','temporal_pbc_jump')}
        with patch('audit_pacer_fkg_g2a_pbc.minimum_image_lengths',return_value=np.asarray([2.,3.8])):
            analyze_frame(ca,[100,100,100,90,90,90],0,None,None,c,x)
        self.assertEqual(c['ca_pbc_split_pairs'],1)
        self.assertEqual(c['ca_non_pbc_severe_pairs'],0)

if __name__=='__main__':unittest.main()
