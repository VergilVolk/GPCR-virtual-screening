"""Execution regression tests; external docking and production launches are mocked."""
import importlib.util,json,sys,tempfile,types,unittest
from pathlib import Path
from unittest.mock import patch,Mock

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pacer_xr_glide_scheduler_v01 as s
import run_pacer_xr_glide_rerun_v01 as runner

class SchedulerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='pacer_scheduler_')
        self.work=Path(self.temp.name); self.ligand=self.work/'prepared.maegz'; self.ligand.write_bytes(b'unchanged ligand')
        self.parent=dict(candidate_id='parent',canonical_smiles='CC',ligand=str(self.ligand),variants=[])
        self.assets={'receptors':{'7TRS':{'grid':'unused','center':[1,2,3],'inner':[10]*3,'outer':[30]*3}}}
        self.rec=dict(receptor_id='7TRS',grid_sha256='grid',center=[1,2,3],inner=[10]*3,outer=[30]*3)
        self.pose=dict(candidate_id='parent',canonical_smiles='CC',prepared_variant_id='parent__s00__v01',molscrub_state_id='parent__s00',entry_index=2,glide_gscore=-7.,receptor_id='7TRS',pose_rank=1)

    def tearDown(self): self.temp.cleanup()

    def fake_bridge(self,work,mode,meta):
        s.write_json(work/'metadata.json',meta)
        if mode=='pack':
            (work/'library.maegz').write_bytes(self.ligand.read_bytes()); return {'counts':{'parent':1}}
        if mode=='parse':
            result={'parent':[self.pose]}; s.write_json(work/'batch_poses.json',result); return result
        raise AssertionError(mode)

    def command(self,work,njobs):
        jobname='unique_'+work.name
        (work/(jobname+'_pv.maegz')).write_bytes(b'valid mock output')
        cmd=dict(exit_code=0,elapsed_seconds=1,system_cpu_percent_mean=1,njobs=njobs)
        s.write_json(work/'glide.command.json',cmd); return cmd,jobname

    def context(self,command=None,bridge=None):
        return patch.multiple(s,science_contract=Mock(return_value={'science':'frozen'}),receptor_meta=Mock(return_value=self.rec),bridge=Mock(side_effect=bridge or self.fake_bridge),run_batch_command=Mock(side_effect=command or self.command))

    def test_completed_campaign_reused_without_launch_or_retry_budget(self):
        root=self.work/'batches'
        with self.context(): record=s.campaign(root,'7TRS',[self.parent],8,self.assets)
        with self.context(),patch.object(s,'remaining_campaign_limit',side_effect=AssertionError('cached result must not consume retries')):
            again=s.campaign(root,'7TRS',[self.parent],16,self.assets,legacy_root=self.work/'legacy')
            s.run_batch_command.assert_not_called(); s.bridge.assert_not_called()
        self.assertEqual(record,again)

    def test_transient_retry_clean_directory_name_and_reduced_pressure(self):
        calls=[]
        def command(work,njobs):
            calls.append((work,njobs))
            if len(calls)==1:
                (work/'glide.console.log').write_text('Job Control connection interrupted')
                raise ValueError('execution failure')
            return self.command(work,njobs)
        with self.context(command): result=s.campaign(self.work/'batches','7TRS',[self.parent],8,self.assets)
        self.assertEqual(result['attempt'],2); self.assertEqual([n for _,n in calls],[8,4])
        self.assertNotEqual(calls[0][0],calls[1][0]); self.assertTrue(calls[0][0].exists())
        for work,_ in calls:
            self.assertEqual((work/'library.maegz').read_bytes(),b'unchanged ligand')
            self.assertIn('PRECISION SP\nPOSES_PER_LIG 5\nPOSTDOCK_NPOSE 5\nFORCEFIELD OPLS_2005', (work/'campaign.in').read_text())

    def test_three_attempt_limit_includes_existing_legacy_attempt(self):
        legacy=self.work/'legacy'; (legacy/'jobs/parent_7TRS/attempt_0001').mkdir(parents=True)
        calls=[]
        def failure(work,njobs):
            calls.append(work); (work/'glide.console.log').write_text('temporary Job Control resource unavailable'); raise ValueError('transient')
        with self.context(failure):
            with self.assertRaises(ValueError): s.campaign(self.work/'batches','7TRS',[self.parent],8,self.assets,legacy_root=legacy)
            self.assertEqual(len(calls),2)
            with self.assertRaisesRegex(ValueError,'budget exhausted'): s.campaign(self.work/'batches','7TRS',[self.parent],8,self.assets,legacy_root=legacy)
            self.assertEqual(len(calls),2)

    def test_retry_budget_cannot_reset_when_batch_membership_changes(self):
        root=self.work/'batches'; old=root/'other_batch'; old.mkdir(parents=True)
        s.write_json(old/'status.json',dict(attempt=2,inputs=dict(receptor_id='7TRS',ligands={'parent':'hash'})))
        legacy=self.work/'legacy'; (legacy/'jobs/parent_7TRS/attempt_0001').mkdir(parents=True)
        with self.assertRaisesRegex(ValueError,'exhausted'): s.remaining_campaign_limit(root,legacy,'7TRS',[self.parent],{})

    def test_parser_failure_not_retried_even_with_transient_warning_in_log(self):
        def command(work,njobs):
            (work/'glide.console.log').write_text('Job Control warning'); return self.command(work,njobs)
        def bridge(work,mode,meta):
            if mode=='parse': raise ValueError('candidate identity mismatch')
            return self.fake_bridge(work,mode,meta)
        with self.context(command,bridge):
            with self.assertRaisesRegex(ValueError,'Scientific/parser'): s.campaign(self.work/'batches','7TRS',[self.parent],8,self.assets)
            self.assertEqual(s.run_batch_command.call_count,1)

    def test_modified_batch_output_invalidates_completion(self):
        with self.context(): record=s.campaign(self.work/'batches','7TRS',[self.parent],8,self.assets)
        Path(next(iter(record['artifacts']))).write_bytes(b'corrupted')
        with self.context(): self.assertFalse(s.batch_valid(record,self.assets,{'science':'frozen'},{'parent':self.parent}))

    def test_all_200_preparation_caches_reused_and_individual_corruption_rejected(self):
        roster=s.load_roster(); assets=s.load_asset_lock(); science=s.science_contract()
        results=[s.valid_preparation(row,s.OUTPUT/'glide',assets,science) for row in roster.itertuples()]
        self.assertEqual(len(results),200)
        row=next(roster.itertuples()); source=s.OUTPUT/'glide/jobs'/('prepare_'+row.candidate_id)/'status.json'
        record=json.loads(source.read_text()); status=self.work/'jobs'/('prepare_'+row.candidate_id)/'status.json'
        record['result']['variants'][0]['candidate_id']='wrong'; s.write_json(status,record)
        with self.assertRaisesRegex(ValueError,'mapping'): s.valid_preparation(row,self.work,assets,science)

    def test_scheduler_changes_preserve_scientific_contract(self):
        lock=s.science_contract()
        self.assertEqual(s.function_hash(s.P/'scripts/run_pacer_xr_glide_rerun_v01.py',['prepare_candidate']),lock['preparation_function_ast_sha256'])

    def test_resource_budget_and_bridge_change_fail_closed(self):
        config=self.work/'scheduler.json'
        for campaigns,njobs in [(11,8),(2,16),(1,17)]:
            s.write_json(config,dict(campaign_concurrency=campaigns,njobs_per_campaign=njobs,batch_bridge_sha256=s.sha(s.BRIDGE)))
            with patch.object(s,'SCHEDULER',config),self.assertRaisesRegex(ValueError,'resource budget'): s.scheduler_contract()
        s.write_json(config,dict(campaign_concurrency=2,njobs_per_campaign=8,batch_bridge_sha256='wrong'))
        with patch.object(s,'SCHEDULER',config),self.assertRaisesRegex(ValueError,'bridge changed'): s.scheduler_contract()

    def test_resume_submits_only_pending_parent_receptor(self):
        roster=s.load_roster().head(2); assets=s.load_asset_lock(); prepared={r.candidate_id:{**self.parent,'candidate_id':r.candidate_id} for r in roster.itertuples()}
        index={cid+'_'+rec:dict(state='COMPLETE') for cid in prepared for rec in assets['receptors']}
        index[roster.iloc[1].candidate_id+'_cluster_06']['state']='FAILED_RETRYABLE'
        config=dict(campaign_concurrency=2,njobs_per_campaign=8)
        with patch.multiple(s,select_scope=Mock(return_value=roster),verify_environment=Mock(),scheduler_contract=Mock(return_value=config),audit_index=Mock(return_value=(prepared,index,{})),campaign=Mock(return_value={'result':{'evaluations':{'parent':{}}}}),collect_complete=Mock(return_value=(prepared,index,{})),emit_channels=Mock(),output_root=Mock(return_value=self.work)),patch.object(s.pd,'read_csv',return_value='channels'):
            self.assertEqual(s.run_optimized(),'channels')
            self.assertEqual(s.campaign.call_count,1)
            args=s.campaign.call_args.args
            self.assertEqual(args[1],'cluster_06'); self.assertEqual([p['candidate_id'] for p in args[2]],[roster.iloc[1].candidate_id])

    def test_original_interrupted_state_classified_without_preparation_or_docking(self):
        roster=s.load_roster(); assets=s.load_asset_lock()
        baseline=json.loads((s.OPT/'baseline.json').read_text()); poses={}
        for item in baseline['evaluations']:
            if item['state']=='COMPLETE':
                result=json.loads(Path(item['status_path']).read_text())['result']
                poses[item['candidate_id']+'_'+item['receptor_id']]=[dict(glide_gscore=result['glide_gscore'],prepared_variant_id='a',entry_index=2)]
        with patch.object(s,'OPT',self.work),patch.object(s,'bridge',return_value={'completed_poses':poses}),patch.object(runner,'prepare_candidate',side_effect=AssertionError('no preparation')):
            prepared,index,results=s.audit_index(s.OUTPUT/'glide',roster,assets,True)
        from collections import Counter
        self.assertEqual(len(prepared),200); self.assertEqual(len(results),18)
        self.assertEqual(Counter(x['state'] for x in index.values()),dict(COMPLETE=18,FAILED_RETRYABLE=3,INCOMPLETE=3,NOT_STARTED=2176))

if __name__=='__main__': unittest.main()
