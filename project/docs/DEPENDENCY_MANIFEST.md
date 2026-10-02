# PACER-M4 runtime dependency manifest v01

Dependency engineering baseline, 2026-10-02. Git HEAD: `d1f145832907c42850baf2e189c7142454f66519`.

M4 route remains `m4_2023_gpcr_loto`. The 2026 family-aug branch is a **non-decisional second opinion**. This work changes no scientific code, thresholds, weights, protocol, results or frozen assets. No training, generation, docking, MD, Stage1/2/3 or PACER-FKG computation was executed.

## Entrypoints and evidence

`requirements.txt` is the canonical Python/core entrypoint. `environment.yml` is its Conda counterpart. This is an evidence-based direct/necessary baseline, **not an all-platform, transitive or build-hash lockfile**. Existing successful environments and their selected distribution/native metadata are captured in `project/results/dependency_runtime_evidence_v01.json`; no `pip freeze` was used.

Production scope is explicit in `project/config/runtime_dependency_map_v01.json`: current generator/library/filter helpers, multi-conformation docking/preparation, 2023 bundle scorer, 2026 projection second opinion, M4-safe router, MD preparation/production/endpoints, frozen C1-BS256 and FKG apply/graph/region code, and current reporting. Recursive imports include nested function imports, internal helpers and explicit Python subprocess helpers. Training-named imported helpers are read as AST only; their fitting functions are never executed. Unreachable historical/experimental/training/audit/test entrypoints are enumerated as excluded in the audit.

Stage3D has a current `STAGE3D_FREEZE_MANIFEST.json` and report, but **no dedicated Stage3D Python implementation at HEAD**. Its route/pose/scaffold handoff is represented by the frozen assets and RDKit/pandas/NumPy dependencies; this task does not invent or reconstruct an implementation. New-molecule full encoder input uses `build_drugclip_inputs.py`, `build_drugclip_gpcr_ensemble_pockets.py` and `extract_drugclip_embeddings_cpu.py`.

Version priority: successful runtime metadata → existing environment definitions → frozen scientific source/docs → distribution Requires-Dist. No versions were selected from latest-package guesses. Installed versions are evidence of the existing runtime, not proof that a fresh Linux solve can reproduce it.

## Separate runtime stages

| Stage | Python / key versions | Canonical installation files | Existing evidence / check |
|---|---|---|---|
| Candidate chemistry, filtering, docking preparation, bundle/router, scaffold/report | 3.9.23; NumPy 1.26.4; torch 2.8.0; RDKit 2025.3.5 | `requirements.txt`, `environment.yml` | Windows `chrm4_vs`; imports/tools PASS; offline pip dry-run PASS |
| Full DrugCLIP base encoding of new inputs | 3.10.21; torch 2.0.1; RDKit 2022.09.5 | `environment-drugclip.yml` plus pinned local sources | WSL `drugclip-cpu`; PARTIAL: missing unimol checkout |
| Four-context MD, system preparation, endpoints | 3.11.16; OpenMM 8.6.1; NumPy 1.26.4 | `environment-md.yml` | WSL `pacer-dc-md`; imports/native-tool presence PASS |
| Frozen C1-BS256 extraction and FKG numerical application | 3.11.16; torch 2.6.0+cu126; NumPy 2.4.6; PyG 2.7.0 | `environment-pacer-fkg.yml`, `requirements-pacer-fkg.txt`, pinned Geom2Vec | Windows `pacer_dc_geom2vec`; imports/compiled modules PASS |

The offline PACER-200 scorer uses precomputed representations and projection weights. It does **not** need DrugCLIP, Uni-Core, LMDB, their 1.18 GB base checkpoint or torch 2.0.1. It is limited to those frozen 200 molecules; it is not a replacement for base encoding new molecules. Full DrugCLIP keeps its successful legacy runtime. NumPy 1 and NumPy 2 / PyTorch ABI stacks remain isolated. No claim is made that these environments cannot ever coexist; their merged compatibility has not been validated.

Old `project/environment*.yml` files and bootstrap scripts are preserved as provenance. They retain historical names/constraints; they are not automatically redirected to the new root definitions. Use the new root files explicitly for a fresh deployment.

## Python dependencies

Core directly/necessarily pins: NumPy 1.26.4, pandas 2.3.1, SciPy 1.13.1, scikit-learn 1.6.1, RDKit 2025.3.5, torch 2.8.0, LightGBM 4.6.0, Meeko 0.7.1, molscrub 0.1.1, Gemmi 0.7.3, matplotlib 3.9.4. Meeko's installed metadata has no Requires-Dist; Gemmi is retained because its installed source imports it. `scrubber → molscrub` is confirmed by top_level metadata. `sklearn → scikit-learn` and `yaml → PyYAML` are confirmed by distribution metadata. OpenFF namespace packages use their actual Conda names, not a fictional `pip install openff`. The core Conda definition installs Meeko through the pip entrypoint because its observed Conda build forces unused optional ProDy; NumPy/RDKit/Gemmi/SciPy are supplied explicitly.

| Runtime layer | Additional direct/necessary dependencies | Role / version source |
|---|---|---|
| MD | OpenMM 8.6.1; openmmforcefields 0.14.1; openff-toolkit-base 0.16.5; openff-interchange-base 0.3.29; openff-units 0.2.2; openff-forcefields 2026.01.0; PDBFixer 1.12; ParmEd 4.3.1; MDAnalysis 2.10.0; importlib_resources 7.1.0 | Runtime/system preparation; successful WSL native/distribution metadata. ParmEd and Interchange are necessary scientific stack support. importlib_resources is an undeclared import of openmmforcefields. |
| Full DrugCLIP | LMDB 2.3.0; biopandas 0.5.2; iopath 0.1.10; ml-collections 1.0.0; tensorboardX 2.6.2.2; tokenizers 0.23.1; tqdm 4.70.0; PyYAML 6.0.3 | Source/registration and package metadata requirements; successful WSL versions. No Jupyter/notebook stack copied. |
| DrugCLIP metadata support | setuptools 80.10.2; looseversion 1.1.2; wandb (unpinned) | setuptools is legacy import/build support. looseversion 1.1.2 is biopandas's exact metadata constraint, not the installed 1.3.0; clean repair not tested. wandb is Uni-Core's declared installation dependency, optional to CPU inference, absent in installed runtime; no tested pin exists. |
| C1-BS256 / Geom2Vec | torch-geometric 2.7.0; torch-cluster 1.6.3+pt26cu126; torch-scatter 2.1.2+pt26cu126; MDAnalysis 2.10.0; mdtraj 1.11.1.post2; einops 0.8.2; tqdm 4.70.1 | Existing Windows wheels plus Geom2Vec metadata. torch-cluster supplies radius_graph's native backend. torch-scatter is a caught optional import in equivariant source; retained as the installed backend baseline. No alternative architecture or pretraining extras added. |

`joblib`, `networkx` and related utilities occur in upstream Requires-Dist/Conda metadata and are installed transitively by scikit-learn, torch, OpenFF or MDAnalysis; they are not falsely promoted to root direct imports. Bio/biopython is not a direct import in this production closure. OpenBabel/openbabel-wheel/obabel is installed historically but is not called/imported by the selected formal code, so is not a mandatory dependency. ProDy is a caught optional Meeko import for unused covalent functionality; it is omitted from the clean core baseline. Installed ProDy 2.4.0 advertises old NumPy/Bio constraints and emits a pkg_resources warning; do not copy it into the deployment environment. No torchmd/pretraining/optimizer/Jupyter/Web extras are included merely because upstream source contains other workflows.

The per-import table (qualified import → distribution → required/optional → usage file/line → declaration → version/metadata evidence) is machine-readable in `dependency_audit_v01.json` and `runtime_dependency_map_v01.json`. Standard library and internal imports are explicitly separated.

## Conda/native and external tools

| Dependency | Purpose / version | Source / current path | Platform and Docker provision |
|---|---|---|---|
| AmberTools | Parameterization executables and native stack; 24.8 | conda-forge; `/root/miniforge3/envs/pacer-dc-md/bin/{tleap,antechamber,parmchk2}` | Linux/WSL; this project's MD environment cannot solve on win-64. Install in scientific image, activate its PATH. No simulation binary is run in this task. |
| OpenMM / OpenFF / PDBFixer / ParmEd / RDKit | MD preparation/native libraries; versions above | conda-forge; canonical `environment-md.yml` | Install with Conda; preserve compatible base Toolkit/Interchange pins. Pip-only requirements cannot reproduce AmberTools, force-field assets or native build selection. |
| Vina | Formal ensemble docking; 1.2.7 | `project/tools/vina_1.2.7.exe`; frozen SHA256 `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5` in Stage3D manifest | Windows binary now. Linux image requires Linux Vina 1.2.7 from a verified upstream build. A pip `vina` Python binding is not the invoked executable and is omitted. |
| mk_prepare_receptor | Receptor PDBQT CLI; Meeko 0.7.1 | Windows `<chrm4_vs>/Scripts/mk_prepare_receptor.exe`; provided by Meeko | Linux CLI/script must be in runtime PATH. Existing preparation code explicitly names `.exe`; this is an unresolved portability blocker, not solved by a requirements marker. |
| Git | Source checkout/build; runtime Geom2Vec source verification | Installed system executable; exact binary version not frozen | Include in source/build stage and FKG runtime where source verification calls it; no Git source fetch performed here. |
| CUDA NVRTC/native CUDA libs | OpenMM CUDA JIT; cuda-version 13.2, cuda-nvrtc 13.2.86 | Conda MD runtime metadata and historical MD environment note | Keep in MD image/runtime layer, not root pip requirements. Driver 596.08 / support ceiling 13.2 is historically recorded, not queried live in this task. |
| NVIDIA host driver / NVIDIA Container Toolkit | GPU device exposure to Linux containers | Host/container infrastructure; version not newly measured | Host driver outside image/requirements. Container Toolkit and `--gpus all` provision separately. Import-only checks do not verify CUDA platform, kernels, memory or driver compatibility. |

The MD SMIRNOFF path invokes AmberTools through its stack for charge/parameter support; no new script call is inferred where it does not exist. `tleap`, `antechamber`, `parmchk2` presence is checked as the native preparation baseline. OpenMM's glibc/native/CUDA requirements are preserved in `dependency_runtime_evidence_v01.json`; choose a Linux base satisfying those Conda build constraints (observed OpenMM build requires glibc >=2.28). Neither CUDA driver nor NVIDIA Container Toolkit belongs in requirements.txt.

## Local source dependencies

| Source import | Frozen identity | Current path / availability | Docker and installation |
|---|---|---|---|
| unimol / DrugCLIP | `bowen-gao/DrugClip@7a3a3fa33673f8668c811790f2e4681c98af44ef` | Expected `project/tools/DrugCLIP`; missing in this checkout and WSL smoke runtime | Supply verified checkout separately. Preserve the already-defined CPU device patch and its audit; this task did not apply it. Pass `--drugclip` or PYTHONPATH. No arbitrary PyPI unimol substitute. |
| unicore / Uni-Core | `dptech-corp/Uni-Core@44f6386f4dcd7137fc1e5d5e768117d635d64a26` | Expected `project/tools/Uni-Core` source missing; installed `unicore 0.0.1` in WSL imports | Existing bootstrap uses `python setup.py install --disable-cuda-ext`. Preserve CPU build choice and verify original source identity before rebuilding. Installed package version alone does not prove source identity. Optional fused CUDA modules are unavailable and CPU fallbacks import successfully. |
| geom2vec | `dinner-group/geom2vec@371d642ec1061664f16e49fcac702d07fc8d0b51` | `C:/projects/geom2vec-source`, HEAD verified; installed distribution 1.0.0 from local path. Checkout has pre-existing untracked build/egg-info. | Mount/copy exact source and install `python -m pip install --no-deps /opt/geom2vec-source` after frozen requirements. Frozen code checks Git commit and inference file equivalence. pyproject version 1.0.0 and historical setup.py 0.2.0 differ; commit/hash governs identity. Source is deliberately outside pip requirement entries. |

Local source build tools (`setuptools`, wheel, and native compiler only if rebuilding an upstream extension) belong to build-time layers. Pip build isolation may resolve additional build tools; this baseline does not claim to lock that toolchain. Use matching prebuilt PyG torch/CUDA/Python ABI wheels and `--only-binary=:all:` during a future fresh test; do not silently build substitute extensions. Torch 2.6.0+cu126 and cluster/scatter local version suffixes are intentional ABI evidence, not Windows-only artifacts.

## Scientific Runtime Assets

Model/checkpoint assets, NPZ representations, CSV benchmarks/handoffs, PDB/PDBQT receptor/pose files, MD states/trajectories, force-field files and FKG freeze/graph/region/contrast manifests are **not requirements**. No asset was downloaded, replaced, regenerated or serialized by this task. Existing expected SHA256 values are copied from their authority files; missing hashes are explicitly unrecorded, not invented. File presence below is inventory only, not model-load or prediction validation.

Full inventory with purpose, exact path, existing SHA256 and its evidence, required status, download/transfer policy, frozen status and presence: `project/results/dependency_scientific_assets_v01.json`.

| Asset | Runtime path | Existing SHA256 | Required / source policy | Frozen / present |
|---|---|---|---|---|
| PACER-200 2023 frozen representations | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_2023_frozen_representations.npz` | `c22bc605d4ed54ec7115062d364419fb781a393b2ed44fe7da081b51af042658` | yes; verified transfer only | yes / True |
| drugclip2023_m4_loto_seed20260925.projection.pt | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260925.projection.pt` | `a954c7c16fe857b5e3b02f1fbfdf7b2b1284d3559d86448b33736e16dd5d7398` | yes; verified transfer only | yes / True |
| drugclip2023_m4_loto_seed20260926.projection.pt | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260926.projection.pt` | `3cc92da9f066a36f9d65b47b7859e69ee7a456d48047a3e1bbd509283d50d4b9` | yes; verified transfer only | yes / True |
| drugclip2023_m4_loto_seed20260927.projection.pt | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260927.projection.pt` | `578ab8349fa8723f5609067e6511ebdfab8541dab8558e883bb7affeb7a8780c` | yes; verified transfer only | yes / True |
| DrugCLIP 2023 base checkpoint | `project/tools/DrugCLIP/artifacts/checkpoint_best.pt` | `dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e` | new-molecule encoding only; upstream or verified transfer; do not download now | yes / False |
| 2026 base representations | `project/results/gpcr_drugclip_screening_v01/science2026_ensemble_embeddings.npz` | `unrecorded (see authority manifest)` | second opinion only; verified transfer only | yes / False |
| 2026 candidate representations | `project/data/drugclip_muscarinic_family_aug_v01/cand200/cands_science2026.npz` | `unrecorded (see authority manifest)` | second opinion only; verified transfer only | yes / False |
| 2026 family-aug seed 20260925 | `project/results/drugclip_science2026/family_aug_v01/science2026_13target_ep80_famaug.seed20260925.projection.pt` | `unrecorded (see authority manifest)` | second opinion only; verified transfer only | yes / True |
| 2026 family-aug seed 20260926 | `project/results/drugclip_science2026/family_aug_v01/science2026_13target_ep80_famaug.seed20260926.projection.pt` | `unrecorded (see authority manifest)` | second opinion only; verified transfer only | yes / True |
| 2026 family-aug seed 20260927 | `project/results/drugclip_science2026/family_aug_v01/science2026_13target_ep80_famaug.seed20260927.projection.pt` | `unrecorded (see authority manifest)` | second opinion only; verified transfer only | yes / True |
| ViSNet frozen checkpoint | `project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth` | `b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417` | yes; pinned upstream/verified transfer; do not download now | yes / True |
| Frozen C1-BS256 specification | `project/encoder_candidate_v01/CANDIDATE_SPEC.json` | `unrecorded (see authority manifest)` | yes; repository | yes / True |
| FKG numerical freeze | `project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json` | `b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd` | yes; verified transfer only | yes / False |
| FKG normalization/RFF files | `project/results/pacer_fkg_v02_longmd_v01/calibration/*` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / directory/package inventory required |
| M4 graph | `project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| M4 regions | `project/results/pacer_dc_four_context_v01/compound110/G2_REGION_MAP_v02.json` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| M4 residue mapping | `project/results/pacer_dc_four_context_v01/compound110/G2_RESIDUE_MAPPING_v01.csv` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| FKG contrasts | `project/results/pacer_fkg_v02_longmd_v01/calibration/V02_CONTRAST_DEFINITIONS.json` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / False |
| FKG protection manifest | `project/pacer_fkg_v02/V01_PROTECTION_MANIFEST_SHA256.txt` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| M4R source ensemble | `project/data/miao2026_m4r/M4R_ensemble_clusters.pdb` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| Prepared ensemble | `project/results/m4_gamd_ensemble/receptors/*` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / directory/package inventory required |
| Common MD protein | `project/results/pacer_dc_common_protein_v01/7TRS_common_protein_pH74.pdb` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / False |
| MD start systems | `project/results/pacer_dc_restraint_release_v01/*` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / directory/package inventory required |
| Reference trajectories | `project/results/pacer_dc_long_md_v02/*` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / directory/package inventory required |
| Closure provenance | `project/results/pacer_dc_close_loop_20ns_v01/CLOSE_LOOP_EXECUTION_MANIFEST_v01.json` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / False |
| Ligand force-field asset | `openff_unconstrained-2.2.1.offxml (inside pinned openff-forcefields)` | `unrecorded (see authority manifest)` | yes; package-provided resource; exact resource presence not rechecked | yes / directory/package inventory required |
| Stage1 potency reference | `project/data/benchmarks/m4_pam_v1/potency_molecules.csv` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| Stage1 activity reference | `project/data/benchmarks/m4_pam_v1/pam_vs_inactive.csv` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / True |
| Generated library files | `project/results/generated/*` | `unrecorded (see authority manifest)` | yes; repository or exact archived transfer; no regeneration | yes / directory/package inventory required |
| Stage3D stage2_input | `project/results/pacer_candidates_v01/predock_portfolio.csv` | `0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b` | yes; verified transfer only | yes / True |
| Stage3D docking_metadata | `project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/metadata.json` | `8d1703f4f3a40b7b73cc6fb74fb7b04f63ade20ba040608e6624607c5b0746d3` | yes; verified transfer only | yes / True |
| Stage3D docking_ledger | `project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/ledger.jsonl` | `829f32c2f8c9f77cee977196b275d52f28129d395d7a28a746fc8e4268ba037a` | yes; verified transfer only | yes / True |
| Stage3D docking_framewise_scores | `project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/framewise_scores.csv` | `7eab001a2bb579e964cb6f93f21b89236f3092991bd449581bc3b74b88682d93` | yes; verified transfer only | yes / True |
| Stage3D docking_statewise_scores | `project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/statewise_scores.csv` | `c9959c05a319f9e7813400381b4179f56cfe0ee380a6c4f3d2205592eb6e7f9e` | yes; verified transfer only | yes / True |
| Stage3D router_input | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_two_model_scores.csv` | `b145b91cfca1e768761113172796be1253ea3a4aa656a5d37df6c9b9c86c4445` | yes; verified transfer only | yes / True |
| Stage3D router_output | `project/results/pacer_stage3_v01/pacer200_m4_safe_routed_ranking.csv` | `209facd6c81d70a48ae3137707e96c329f39b8960a784ec977bc17a83a1c4f9d` | yes; verified transfer only | yes / True |
| Stage3D stage3a_docking_qc | `project/results/pacer_stage3_v01/stage3a_docking_qc.json` | `935105216427530a85677b760726847185a2176157ec0ac01c7fdcdf44638f33` | yes; verified transfer only | yes / True |
| Stage3D structural_gate | `project/results/pacer_stage3_v01/pacer200_structural_gate.csv` | `b0dc774888051af3d9821df2eda0b07d9bc64e2432f1cd28b95960ba77098f81` | yes; verified transfer only | yes / True |
| Stage3D merged | `project/results/pacer_stage3_v01/pacer200_stage3d_merged.csv` | `72ef23d31cb228437228f048f5b7d6ba1e6b277465a615b6207669187e3428d6` | yes; verified transfer only | yes / True |
| Stage3D diverse_shortlist | `project/results/pacer_stage3_v01/pacer_stage3d_diverse_shortlist.csv` | `caaab138384226155fcb5e6737e4c04c5fa281122500c9530bd795007efc7974` | yes; verified transfer only | yes / True |
| Stage3D pose_manifest | `project/results/pacer_stage3_v01/md_shortlist_pose_manifest.csv` | `44d11f1b8a41b81a51bc010bd1b94cc0a539fc8fec4c47cef211d29da3f0b68d` | yes; verified transfer only | yes / True |
| PACER0100 pose | `project/results/pacer_stage3_v01/md_shortlist_poses/PACER0100_c00_s00.pdbqt` | `8b26a988d128901e7d09049cfe9c35f591137b0aa2ea6991c56686d14554bc5b` | yes; verified transfer only | yes / True |
| PACER0047 pose | `project/results/pacer_stage3_v01/md_shortlist_poses/PACER0047_c00_s00.pdbqt` | `228d2a6f0ec4afceb19bbd72920cb88699b27129de43e032199765e50bfd647b` | yes; verified transfer only | yes / True |
| PACER0046 pose | `project/results/pacer_stage3_v01/md_shortlist_poses/PACER0046_c06_s00.pdbqt` | `4b50a358ff7f6addcb5c688db8b536ffa867303c604c1210636c3c64c768e75e` | yes; verified transfer only | yes / True |
| Stage3D freeze manifest | `project/results/pacer_stage3_v01/STAGE3D_FREEZE_MANIFEST.json` | `unrecorded (see authority manifest)` | yes; repository/archive | yes / True |
| Bundle manifest | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/bundle_manifest.json` | `unrecorded (see authority manifest)` | yes; repository | yes / True |

The current checkout is not a self-contained full scientific deployment. In particular, the FKG numerical freeze/contrast files, MD start/provenance assets and base DrugCLIP checkpoint/source need exact archive provisioning where marked absent. Existing frozen MD reference/cache hierarchies need transfer/volume inventories; do not regenerate missing calibration or reference data. 2026 scorer's existing `D:\CLC` absolute paths and some source defaults are additional portability constraints. These are documented, not edited.

## Windows and Linux/Docker handling

Windows core works against the existing Vina `.exe` and Meeko launcher. Linux cannot run a copied `.exe`; the current docking script hardcodes its Windows Vina path, and receptor preparation explicitly invokes an `.exe`. A later, authorized deployment adapter/configuration change is required. This dependency-only task deliberately leaves those scientific entrypoints untouched. `rescore_cands200_famaug.py` hardcodes a `D:/CLC` tree and imports its Proj helper; it was statically scanned, never imported/executed. Frozen encoder defaults include Windows source paths; use existing CLI source options where available. No repository restructuring was done.

No Dockerfile exists in the current tracked repository; no Dockerfile, image build, upload or deployment was created. Recommended layer order:

1. Base scientific Linux image with compatible glibc and pinned Conda/Python native stage runtime.
2. Per-stage Python/Conda dependencies, cached by the corresponding environment/requirements files. Keep legacy DrugCLIP, core, MD and FKG as separate images or environment prefixes.
3. Native chemistry/MD tools and ABI-matched PyG wheels; stable Vina/AmberTools source/build identities belong in this cacheable layer.
4. Pinned local source dependencies and project source. Cache source revisions separately; preserve CPU patch audit and inference-source verification.
5. Frozen scientific assets as verified read-only mounts/object-storage artifacts; avoid baking large base checkpoints/trajectories into Git or frequently changing source layers. Small checked-in bundle manifests are metadata, not substitutes for full runtime assets.
6. Runtime entrypoint and writable output/cache volumes on native Linux filesystem. Activate the chosen environment's PATH. Configure host GPU driver/Container Toolkit separately.

Dependency structure is **PARTIAL** for Docker. Remaining blockers: clean per-platform solves/wheel availability and build locks; missing full DrugCLIP/Uni-Core source and scientific assets; existing Windows/hardcoded path assumptions; unpinned Uni-Core metadata wandb repair and untested looseversion correction; fresh Linux FKG ABI/GPU platform validation. These are not solved by an import PASS. Web/API dependency layer will be added when deployment interface is frozen.

## Validation and maintenance

Checks performed here: AST coverage, existing-stage import/presence smoke checks, Python compileall for the two new scripts, and `python -m pip install --dry-run --no-index -r requirements.txt` in the existing chrm4_vs environment. Dry-run makes no changes and succeeds because the core pins and transitive dependencies are already installed; it does **not** prove availability of every clean-install wheel. Conda/native dependencies cannot be validated by pip dry-run. No environment was installed, upgraded, downgraded or recreated.

`python project/scripts/check_runtime_dependencies.py` defaults to an overview of all runtimes. In a split environment, PARTIAL is expected. Use `--profile core|drugclip|md|fkg` with each runtime's own interpreter for actionable stage status. Imports of packages only, including `openmmforcefields.generators` and Geom2Vec ViSNet, are performed; project scripts are never imported. Native executables are only located, not invoked; their displayed versions are recorded baseline values, not live --version results. The checker does not validate GPU contexts, force-field application, checkpoint predictions or model asset contents.

Existing metadata inconsistencies: DrugCLIP `pip check` reports missing wandb and biopandas's `looseversion==1.1.2` versus installed 1.3.0. Imports pass for installed Uni-Core/biopandas; the new legacy YAML explicitly declares the metadata repairs with their validation limitations. MD `pip check` reports installed mdtraj 1.11.1 requiring NumPy ~=2.0 versus NumPy 1.26.4. Formal MD code does not import mdtraj, so the new MD definition omits it. FKG keeps its successful separate NumPy 2 + mdtraj stack. Installed core ProDy constraints are historical optional baggage and are not frozen into the new core entrypoint.

Future validation commands below are examples for a **new disposable environment**; they were not executed in this task. Never run bootstrap wrappers that download assets, retrain or execute pipelines as a dependency smoke check.

```bash
conda env create -f environment.yml
conda run -n pacer-core-frozen-v01 python project/scripts/check_runtime_dependencies.py --profile core
conda env create -f environment-md.yml  # Linux/WSL
conda run -n pacer-md-frozen-v01 python project/scripts/check_runtime_dependencies.py --profile md
conda env create -f environment-drugclip.yml
# Provision and verify pinned local DrugCLIP + CPU-built Uni-Core source separately.
conda run -n pacer-drugclip-legacy-frozen-v01 python project/scripts/check_runtime_dependencies.py --profile drugclip --source drugclip=/opt/DrugCLIP --source unicore=/opt/Uni-Core
conda env create -f environment-pacer-fkg.yml
conda run -n pacer-fkg-frozen-v01 python -m pip install --no-deps /opt/geom2vec-source
conda run -n pacer-fkg-frozen-v01 python project/scripts/check_runtime_dependencies.py --profile fkg --source geom2vec=/opt/geom2vec-source
python project/scripts/audit_runtime_dependencies.py --geom2vec-source /opt/geom2vec-source
python -m compileall -q project/scripts/check_runtime_dependencies.py project/scripts/audit_runtime_dependencies.py
```

After clean solves and native validation, retain per-platform Conda explicit/build locks and pip wheel hashes in a separate dependency-only change. Until then the direct pins, selected build evidence, source commits and frozen asset manifests are the maintainable baseline; they do not certify a clean full deployment.
