"""Contract tests use synthetic energies only; never write prospective ranking files."""
import contextlib,io,importlib.util,json,math,sys,tempfile,types,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pacer_xr_production_v01 as c
import run_pacer_xr_full_rerun_v01 as master

class Contracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.roster=c.load_roster()
        cls.temp=tempfile.TemporaryDirectory(prefix='pacer_xr_contract_')
        cls.work=Path(cls.temp.name)

    @classmethod
    def tearDownClass(cls): cls.temp.cleanup()

    def fixture(self,n=200):
        roster=self.roster.head(n).copy(); frame=roster.copy()
        for index,channel in enumerate(c.CHANNELS): frame[channel+'_raw']=-10+np.arange(n)/n+index*.1
        return roster,frame

    def test_source_hash_and_roster(self):
        self.assertEqual(c.sha(c.SOURCE),c.SOURCE_SHA)
        self.assertEqual(len(self.roster),200)
        self.assertEqual(self.roster.canonical_smiles.nunique(),200)
        self.assertEqual(self.roster.rerun_final_rank.tolist(),list(range(1,201)))

    def test_hash_mismatch_fails_before_normalization(self):
        path=self.work/'wrong_source.csv'; path.write_text('bad')
        with patch.object(c,'SOURCE',path),self.assertRaisesRegex(ValueError,'SHA256'): c.load_roster()

    def test_schema_invalid_cases(self):
        for field,value in [('gen_id',self.roster.iloc[1].gen_id),('smiles',self.roster.iloc[1].smiles),('smiles',None),('final_rank',2),('score',float('inf'))]:
            frame=self.roster[['gen_id','smiles','score','rankpct','final_rank']].copy(); frame.loc[0,field]=value
            with self.subTest(field=field,value=value),self.assertRaises(ValueError): c.validate_roster(frame)
        with self.assertRaises(ValueError): c.validate_roster(self.roster.drop(columns=['score']))
        with self.assertRaises(ValueError): c.validate_roster(self.roster.head(3))

    def test_structural_anchors(self):
        anchors=c.verify_anchors(self.roster)
        self.assertEqual([(x['candidate_id'],x['rerun_final_rank']) for x in anchors],[('PACERGEN01755',17),('PACERGEN00094',60),('PACERGEN00123',114)])
        bad=self.roster.copy(); bad.loc[bad.candidate_id=='PACERGEN01755','canonical_smiles']='C'
        with self.assertRaises(ValueError): c.verify_anchors(bad)

    def test_pmf(self):
        pmf=c.load_pmf(); self.assertEqual(set(pmf),set(range(10)))
        self.assertAlmostEqual(pmf[9],5.219181698603361)

    def test_glide_reduction_ties(self):
        poses=[dict(glide_gscore=-7,prepared_variant_id='b',entry_index=1),dict(glide_gscore=-7,prepared_variant_id='a',entry_index=4),dict(glide_gscore=-7,prepared_variant_id='a',entry_index=2)]
        self.assertEqual(c.reduce_glide(poses)['entry_index'],2)
        for invalid in [[],[dict(glide_gscore=float('nan'),prepared_variant_id='a',entry_index=1)]]:
            with self.assertRaises(ValueError): c.reduce_glide(invalid)

    def bridge_parse(self,structures,name):
        work=self.work/name; work.mkdir(exist_ok=True)
        meta=dict(candidate_id='PACERGEN02353',canonical_smiles='CC',receptor_id='cluster_00',grid_sha256='hash',center=[1,2,3],poseviewer='unused',variants=[dict(prepared_variant_id='PACERGEN02353__s00__v01',molscrub_state_id='PACERGEN02353__s00')])
        c.write_json(work/'metadata.json',meta)
        module=types.ModuleType('schrodinger'); module.structure=types.SimpleNamespace(StructureReader=lambda path:iter(structures))
        structutils=types.ModuleType('schrodinger.structutils'); structutils.analyze=types.SimpleNamespace()
        with patch.dict(sys.modules,{'schrodinger':module,'schrodinger.structutils':structutils}):
            spec=importlib.util.spec_from_file_location('test_sch_bridge',c.P/'scripts/pacer_xr_schrodinger_v01.py'); bridge=importlib.util.module_from_spec(spec); spec.loader.exec_module(bridge)
            with patch.object(sys,'argv',['bridge','parse',str(work/'metadata.json')]): bridge.main()
        return json.loads((work/'poses.json').read_text())

    def bridge_fixture(self):
        rec=types.SimpleNamespace(property={'b_glide_receptor':True,**{'r_glide_gridbox_'+a+'cent':x for a,x in zip('xyz',[1,2,3])}})
        pose=types.SimpleNamespace(title='PACERGEN02353__s00__v01',atom=[types.SimpleNamespace(xyz=[1.,2.,3.])],property=dict(s_pacer_candidate_id='PACERGEN02353',s_pacer_canonical_smiles='CC',s_pacer_prepared_variant_id='PACERGEN02353__s00__v01',r_i_glide_gscore=-7.,r_i_glide_emodel=-999.))
        return rec,pose

    def test_glide_parser_reads_gscore_not_emodel(self):
        rec,pose=self.bridge_fixture(); out=self.bridge_parse([rec,pose],'bridge_valid')
        self.assertEqual(out[0]['glide_gscore'],-7)
        self.assertEqual(out[0]['entry_index'],2)

    def test_glide_parser_rejects_identity_geometry_scores(self):
        for key,value in [('s_pacer_candidate_id','wrong'),('s_pacer_prepared_variant_id','wrong'),('s_pacer_canonical_smiles','CCC'),('r_i_glide_gscore',float('nan'))]:
            rec,pose=self.bridge_fixture(); pose.property[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError): self.bridge_parse([rec,pose],'bridge_bad_'+key)
        rec,pose=self.bridge_fixture(); pose.atom[0].xyz[0]=float('inf')
        with self.assertRaises(ValueError): self.bridge_parse([rec,pose],'bridge_bad_xyz')
        with self.assertRaises(ValueError): self.bridge_parse([pose],'bridge_no_rec')

    def test_vina_parser_and_state_reduction(self):
        path=self.work/'pose.pdbqt'; path.write_text('REMARK VINA RESULT: -8.1 0 0\nREMARK VINA RESULT: -7.5 1 2\n')
        self.assertEqual(c.parse_vina(path),-8.1)
        states=[dict(vina_affinity=-8,state_id='b'),dict(vina_affinity=-8,state_id='a')]
        self.assertEqual(c.reduce_vina(states)['state_id'],'a')
        for text in ['', 'REMARK VINA RESULT: nan 0 0\n']:
            path.write_text(text)
            with self.assertRaises(ValueError): c.parse_vina(path)
        with self.assertRaises(ValueError): c.reduce_vina([])

    def test_both_engines_be_arithmetic(self):
        from analyze_m4_gamd_ensemble import aggregate
        for key in ['glide_gscore','vina_affinity']:
            rows=[dict(cluster_id=i,**{key:-10.+i*.25}) for i in range(10)]
            result=c.ensemble_score(rows,key); pmf=c.load_pmf(); expected=[-10+i*.25+pmf[i] for i in range(10)]
            self.assertEqual(result['BEmin'],min(expected)); self.assertAlmostEqual(result['BEavg'],sum(expected)/10)
            adapted=pd.DataFrame([dict(molecule_id='synthetic',cluster=i,vina_affinity=-10+i*.25,native_pocket_coverage=0) for i in range(10)])
            existing=aggregate(adapted).iloc[0]
            self.assertAlmostEqual(result['BEavg'],existing.BE_avg); self.assertAlmostEqual(result['BEmin'],existing.BE_min)
            for bad in [rows[:-1],rows[:-1]+[rows[0]],rows+[rows[0]]]:
                with self.assertRaises(ValueError): c.ensemble_score(bad,key)

    def test_six_channel_join(self):
        roster,frame=self.fixture(3)
        glide=frame[['candidate_id','canonical_smiles']+[x+'_raw' for x in c.CHANNELS[:3]]]
        vina=frame[['candidate_id','canonical_smiles']+[x+'_raw' for x in c.CHANNELS[3:]]]
        out=c.assemble(roster,glide.iloc[::-1],vina.iloc[::-1]); self.assertEqual(out.candidate_id.tolist(),roster.candidate_id.tolist())
        for channel in c.CHANNELS: self.assertEqual(out[channel+'_raw'].tolist(),frame[channel+'_raw'].tolist())

    def test_channel_completeness_and_identity_rejection(self):
        roster,frame=self.fixture(3)
        bads=[frame.iloc[:-1],pd.concat([frame,frame.iloc[:1]]),frame.drop(columns=['Glide_BEavg_raw'])]
        bad=frame.copy(); bad.loc[0,'canonical_smiles']='C'; bads.append(bad)
        bad=frame.copy(); bad.loc[0,'candidate_id']='UNKNOWN'; bads.append(bad)
        for value in [np.nan,np.inf]:
            bad=frame.copy(); bad.loc[0,'Vina_BEavg_raw']=value; bads.append(bad)
        for bad in bads:
            with self.assertRaises(ValueError): c.validate_scores(bad,roster)

    def fixture_reference_lock(self,name,values):
        work=self.work/name; work.mkdir(exist_ok=True); channels={}
        for channel in c.CHANNELS:
            path=work/(channel+'.csv'); pd.DataFrame({'score':values}).to_csv(path,index=False)
            channels[channel]=dict(path=str(path),column='score',sha256=c.sha(path),finite_rows=int(np.isfinite(values).sum()))
        lock=work/'lock.json'; c.write_json(lock,dict(channels=channels)); return lock,channels

    def test_reference_loader_and_tampering(self):
        lock,channels=self.fixture_reference_lock('refs_valid',np.array([-8.,-7.,np.nan]))
        arrays,_=c.load_references(lock); self.assertEqual(arrays['Glide_PDB'].tolist(),[-8.,-7.])
        Path(channels['Glide_PDB']['path']).write_text('score\n-999\n')
        with self.assertRaisesRegex(ValueError,'REFERENCE_RESTORE_BLOCKED'): c.load_references(lock)
        lock,_=self.fixture_reference_lock('refs_empty',np.array([np.nan]))
        with self.assertRaises(ValueError): c.load_references(lock)

    def test_real_reference_lock(self):
        refs,lock=c.load_references(); self.assertEqual(set(refs),set(c.CHANNELS)); self.assertEqual(lock['reference_commit'],'44798c841ee77230b1f89fe41074e41b458c7070')
        self.assertTrue(all(len(a)>100000 for a in refs.values()))

    def test_empirical_percentile_exact_ties(self):
        values=c.empirical_rank(pd.Series([-3.,-2.,-1.,0.,1.]),np.array([-2.,-1.,-1.,0.]))
        np.testing.assert_allclose(values,[1.,.875,.5,.125,0.])

    def test_fusion_mean_cascade_and_stable_order(self):
        roster,frame=self.fixture(30); refs={x:np.arange(-20.,0.) for x in c.CHANNELS}
        frame.loc[:2,[x+'_raw' for x in c.CHANNELS]]=-20.
        # Stage1 ties remain in frozen rerun order even when the input is shuffled.
        ranked,threshold=c.fuse(frame.iloc[::-1],roster,refs)
        self.assertEqual(threshold,-19.)
        self.assertEqual(ranked.iloc[:3].candidate_id.tolist(),roster.iloc[:3].candidate_id.tolist())
        self.assertTrue((ranked.iloc[:3].cascade_stage==1).all())
        np.testing.assert_allclose(ranked.PACER_XR,ranked[[x+'_rank' for x in c.CHANNELS]].mean(axis=1))
        remaining=ranked[ranked.cascade_stage==2]; self.assertTrue(remaining.PACER_XR.is_monotonic_decreasing)

    def test_natural_top20_and_anchor_insertion_do_not_change_ranking(self):
        roster,frame=self.fixture(); frame['candidate_cascade_rank']=np.arange(1,201); frame['cascade_stage']=2; frame['PACER_XR']=.5
        before=frame.copy(deep=True); natural=c.natural_top20(frame); final=c.shortlist20(frame)
        self.assertEqual(natural.candidate_id.tolist(),roster.head(20).candidate_id.tolist())
        self.assertEqual(len(final),20); self.assertEqual(final.canonical_smiles.nunique(),20)
        self.assertTrue({v[0] for v in c.ANCHORS.values()}.issubset(set(final.candidate_id)))
        self.assertEqual(set(final[final.selection_reason=='legacy_anchor_forced'].candidate_id),{'PACERGEN00094','PACERGEN00123'})
        self.assertFalse(set(roster.iloc[18:20].candidate_id)&set(final.candidate_id))
        pd.testing.assert_frame_equal(frame,before)
        self.assertEqual(final.shortlist_position.tolist(),list(range(1,21)))

    def test_selection_all_anchors_already_natural(self):
        _,frame=self.fixture(); ids={v[0] for v in c.ANCHORS.values()}
        frame=frame.assign(anchor=frame.candidate_id.isin(ids)).sort_values('anchor',ascending=False,kind='stable').drop(columns=['anchor']); frame['candidate_cascade_rank']=np.arange(1,201)
        final=c.shortlist20(frame); self.assertEqual(sum(final.selection_reason=='legacy_anchor_forced'),0)
        self.assertEqual(sum(final.selection_reason=='natural_pacer_xr_top20_legacy_anchor'),3)

    def test_selection_all_three_anchors_forced(self):
        _,frame=self.fixture(); ids={v[0] for v in c.ANCHORS.values()}
        frame=frame.assign(anchor=frame.candidate_id.isin(ids)).sort_values('anchor',ascending=True,kind='stable').drop(columns=['anchor']); frame['candidate_cascade_rank']=np.arange(1,201)
        natural=c.natural_top20(frame); final=c.shortlist20(frame)
        self.assertEqual(sum(final.selection_reason=='legacy_anchor_forced'),3)
        self.assertEqual(set(natural.tail(3).candidate_id),set(natural.candidate_id)-set(final.candidate_id))
        self.assertEqual(len(final),20)

    def test_resume_input_and_output_hashes(self):
        store=c.JobStore(self.work/'resume'); calls=[]
        def operation(work):
            calls.append(1); artifact=work/'score.json'; c.write_json(artifact,dict(score=-7)); return dict(score=-7),[artifact]
        store.execute('job',dict(input='A'),operation); store.execute('job',dict(input='A'),operation)
        self.assertEqual(len(calls),1); self.assertEqual(store.counts['resumed'],1)
        status=json.loads((store.root/'job/status.json').read_text()); Path(next(iter(status['artifacts']))).write_text('corrupt')
        store.execute('job',dict(input='A'),operation); store.execute('job',dict(input='B'),operation)
        self.assertEqual(len(calls),3); self.assertTrue((store.root/'job/attempt_0001').exists())

    def test_failure_tracking_and_retry(self):
        store=c.JobStore(self.work/'failure')
        def failure(work): raise RuntimeError('technical failure')
        with self.assertRaises(RuntimeError): store.execute('job',{'input':'A'},failure)
        status=json.loads((store.root/'job/status.json').read_text()); self.assertEqual(status['status'],'FAILED'); self.assertIn('technical failure',status['error'])
        def success(work):
            path=work/'result'; path.write_text('valid'); return {'score':-5},[path]
        store.execute('job',{'input':'A'},success)
        self.assertEqual(json.loads((store.root/'job/status.json').read_text())['attempt'],2)

    def test_stage_resume_fails_closed(self):
        self.assertFalse(master.stage_ready(self.work/'absent','glide',self.roster))

    def test_stage_resume_requires_every_job_identity(self):
        roster,frame=self.fixture(3); root=self.work/'stage_complete'; engine=root/'glide'; jobs=engine/'jobs'; jobs.mkdir(parents=True,exist_ok=True)
        scores=frame[['candidate_id','canonical_smiles']+[x+'_raw' for x in c.CHANNELS[:3]]]
        c.write_csv(engine/'glide_channels.csv',scores)
        c.write_json(engine/'channels_provenance.json',dict(complete=True,n_candidates=3,evaluations=33,source_sha256=c.SOURCE_SHA,channel_sha256=c.sha(engine/'glide_channels.csv'),codes=c.code_hashes(),assets_lock_sha256=c.sha(c.P/'config/pacer_xr_assets_v01.json')))
        artifact=root/'valid_artifact'; artifact.write_text('valid')
        for cid in roster.candidate_id:
            keys=['prepare_'+cid]+[cid+'_'+rec for rec in ['7TRS']+[f'cluster_{i:02d}' for i in range(10)]]
            for key in keys: c.write_json(jobs/key/'status.json',dict(status='PASS',fingerprint='synthetic',result={'score':-7},artifacts={str(artifact):c.sha(artifact)}))
        self.assertTrue(master.stage_ready(root,'glide',roster))
        (jobs/(roster.iloc[0].candidate_id+'_cluster_09')/'status.json').unlink()
        self.assertFalse(master.stage_ready(root,'glide',roster))

    def test_cascade_inclusive_threshold_and_existing_implementation(self):
        import apply_pacer_xr_candidates as legacy
        roster,frame=self.fixture(30)
        refs={x:np.array([-20.,-19.,-18.,-17.]) for x in c.CHANNELS}
        frame.loc[0,'Glide_BEmin_raw']=-19.
        frame.loc[1,'Glide_BEmin_raw']=-19.+1e-6
        ranked,threshold=c.fuse(frame,roster,refs)
        self.assertEqual(threshold,-19.)
        self.assertEqual(ranked.set_index('candidate_id').loc[roster.iloc[0].candidate_id,'cascade_stage'],1)
        self.assertEqual(ranked.set_index('candidate_id').loc[roster.iloc[1].candidate_id,'cascade_stage'],2)
        work=self.work/'legacy_comparison'; work.mkdir(exist_ok=True); paths={}
        for channel,values in refs.items():
            path=work/(channel+'.csv'); pd.DataFrame({'raw':values}).to_csv(path,index=False); paths[channel]=(path,'raw')
        frame.to_csv(work/'scores.csv',index=False)
        with patch.object(legacy,'FILES',paths),patch.object(legacy,'OUT',work),contextlib.redirect_stdout(io.StringIO()): legacy.apply_scores(work/'scores.csv')
        old=pd.read_csv(work/'candidate_pacer_xr_ranking.csv')
        self.assertEqual(old.candidate_id.tolist(),ranked.candidate_id.tolist())
        np.testing.assert_allclose(old.PACER_XR,ranked.PACER_XR)
        self.assertEqual(old.cascade_stage.tolist(),ranked.cascade_stage.tolist())

    def test_selection_missing_anchor_or_too_short_rejected(self):
        _,frame=self.fixture(); frame['candidate_cascade_rank']=np.arange(1,201)
        with self.assertRaises(ValueError): c.natural_top20(frame.head(19))
        reduced=frame[frame.candidate_id!='PACERGEN00123'].copy(); reduced['candidate_cascade_rank']=np.arange(1,200)
        with self.assertRaises(ValueError): c.shortlist20(reduced)

    def test_production_fusion_cannot_use_acceptance(self):
        from run_pacer_xr_rerun_v01 import run_fusion
        with self.assertRaises((FileNotFoundError,ValueError)): run_fusion()
        for name in ['pacer_xr_rerun_full200_ranking.csv','pacer_xr_natural_top20.csv','pacer_xr_final_shortlist20.csv']:
            self.assertFalse((c.OUTPUT/name).exists())

if __name__=='__main__': unittest.main()
