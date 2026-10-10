# PACER Stage4 Rerun v01 deployment

Replace USER, HOST and REMOTE_PARENT with your actual SSH destination; none is assumed or configured by this preparation.

```powershell
scp C:/projects/GPCR-virtual-screening/project/results/pacer_stage4_rerun_v01/pacer_stage4_rerun_v01_server_bundle.tar.gz USER@HOST:REMOTE_PARENT/
```

On Linux, use a new directory. Never extract into an existing campaign.

```bash
mkdir pacer_stage4_rerun_v01_deployment
cd pacer_stage4_rerun_v01_deployment
tar -xzf ../pacer_stage4_rerun_v01_server_bundle.tar.gz
cd pacer_stage4_rerun_v01
sha256sum -c SHA256SUMS.txt
conda env create -n pacer-stage4-rerun-v01 -f environment-linux.yml
conda activate pacer-stage4-rerun-v01
python scripts/preflight_pacer_stage4_rerun_v01.py --environment-only
python scripts/build_pacer_stage4_rerun_v01.py --all
python scripts/preflight_pacer_stage4_rerun_v01.py --systems
python scripts/schedule_pacer_stage4_rerun_v01.py --phase prepare
```

The inherited environment file pins CUDA NVRTC 13.2. This is an installation recipe, not a claim that the server driver supports it. Both GPU contexts must execute finite-force integration before continuing. An incompatible solver/driver result is a blocker; do not silently change the scientific stack.

Only after reviewing all 36 preparation receipts, the operator may manually start production:

```bash
python scripts/schedule_pacer_stage4_rerun_v01.py --phase production --authorize-production
```

No production command was executed locally. Scheduler uses two workers, physical GPU indices 0 and 1, with CUDA_VISIBLE_DEVICES required unset, and an exclusive scheduler lock. Each task holds a separate lock. Re-run the same command to resume. A completed task is authenticated and skipped. Resume verifies job, config, system XML, PDB topology, starting state and checkpoint hashes. Production saves immutable 100 ps segments, each with its own DCD, CSV, checkpoint and committed receipt; an interrupted segment is preserved as an orphan and replayed from the last committed checkpoint, never appended to a potentially longer DCD. Concatenate only committed segments in numeric order. No XML fallback is silently used for incompatible binary checkpoints.

Outputs: systems/SYSTEM, smoke/SYSTEM, equilibration/SYSTEM/replica_NN, production/SYSTEM/replica_NN/segments/END_STEP. logs contains GPU preflight, scheduler/task logs and readiness receipts. All paths resolve from package root. No repository checkout or external private assets are required. AM1-BCC resources come from packaged conda dependencies; force-field asset files and hashes are recorded at build time.
