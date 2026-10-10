"""Two physical GPUs, one task/device, fail closed, explicit production authorization."""
import argparse
import fcntl
import os
import queue
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from build_pacer_stage4_rerun_v01 import ROOT,table


def main():
    p=argparse.ArgumentParser(); p.add_argument('--phase',choices=['prepare','production'],required=True)
    p.add_argument('--authorize-production',action='store_true'); a=p.parse_args()
    assert a.phase!='production' or a.authorize_production, 'Manual production authorization required'
    assert not os.environ.get('CUDA_VISIBLE_DEVICES'), 'Unset CUDA_VISIBLE_DEVICES to use physical indices 0/1'
    lock=(ROOT/'logs/scheduler_pacer_stage4_rerun_v01.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    # Revalidate both actual CUDA devices on every scheduler launch; failure prevents workers starting.
    subprocess.run([sys.executable,str(ROOT/'scripts/preflight_pacer_stage4_rerun_v01.py'),'--environment-only'],check=True)
    jobs=table('STAGE4_RERUN_JOB_MATRIX.csv')
    if a.phase=='production':
        import json
        for job in jobs:
            ready=ROOT/'equilibration'/job['system_id']/f"replica_{int(job['replica']):02d}"/'ready_pacer_stage4_rerun_v01.json'
            assert json.loads(ready.read_text())['passed'], 'All 36 preparation gates required'
    work=queue.Queue()
    for job in jobs: work.put(job)
    stopped=threading.Event(); errors=[]; mutex=threading.Lock()
    def worker(device):
        while not stopped.is_set():
            try: job=work.get_nowait()
            except queue.Empty: return
            command=[sys.executable,str(ROOT/'scripts/run_pacer_stage4_rerun_v01.py'),'--job',job['job_id'],'--device',str(device),'--phase',a.phase]
            if a.authorize_production: command.append('--authorize-production')
            logfile=ROOT/'logs'/f"{job['job_id']}__{a.phase}.log"
            with logfile.open('a') as log: result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode:
                stopped.set()
                with mutex: errors.append(job['job_id'])
                return
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(worker,i) for i in [0,1]]
        for future in futures: future.result()
    if errors: raise SystemExit('FAILED; no further jobs scheduled: '+', '.join(errors))


if __name__=='__main__': main()
