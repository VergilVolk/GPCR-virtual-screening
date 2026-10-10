# PACER Stage4 Rerun v01 preflight

Status: READY_FOR_SERVER_PREFLIGHT; NOT READY_FOR_PRODUCTION. Production not started.

PASS: frozen ranking SHA256; fixed Cascade ranks 1/2/3; cluster_01 BEmin and scoring provenance; unique source MAEGZ entries; prepared-state/parent chemical identity; formal charges; full source atom mapping; original coordinate preservation; poseviewer receptor coordinate frame; OPM rigid transform and 270-CA correspondence; 274-residue capped receptor baseline; explicit historical ACh remapping and +1 charge; native ACh exclusion; 12 contexts and 36 unique paired-seed task paths; all packaged input files; historical directory inventory unchanged.

Initial candidate/protein and candidate/ACh heavy-atom overlap screen passed at 1 A severe-clash threshold. This is not full system overlap QC. Lipid/water/ion topology, charge sums, disulfides, force-field template completeness and finite energies MUST pass the Linux builder and CUDA system smoke before simulation eligibility.

Actual constructed membrane systems: 0/12. Windows preparation Python has RDKit but lacks OpenMM/OpenFF/PDBFixer/openmmforcefields. No parameterization or successful system build is claimed. AM1-BCC and all 12 builds remain server work.

All candidates use the same cluster_01 receptor. A/P controls are physically shareable, but this package retains 12 named systems and 36 separately executed tasks. Identical controls across candidate families must not be pooled as nine independent receptor conditions.

Scientific parameters follow the supplied historical implementations: ff14SB/lipid17/TIP3P/OpenFF 2.2.1, AM1-BCC, POPC, 0.15 M, 1 nm padding/cutoff, PME, 300 K, 1 bar, HBonds, 2 fs production. Smoke 500 steps/0.5 fs; NVT/NPT each 1 ps/0.5 fs; release 500/100/10/0 each 0.5 ps/1 fs. These are inherited short equilibration stages, not convergence evidence. Trajectory output is explicitly 10 ps, following the user and final handoff correction of historical 50 ps metadata.

Pending server gates: actual RTX 5090 driver/NVRTC compatibility on BOTH devices; dependency solve; AM1-BCC; 12 builds; full topology/disulfide/clash validation; CUDA 500-step system smokes; 36 per-replica equilibration/release audits. Any failure stops the workflow. No automatic scientific parameter changes.

Concrete charge note: the winning Glide states for PACERGEN00123 and PACERGEN00462 have formal charge +2; this is a valid microstate relation to the frozen neutral canonical parent. ACh is +1. The inherited builder neutralizes the membrane baseline before inserting ligands, so final systems can have nonzero net charge. The new builder records total system charge and ion counts; it does not silently add counterions or change the inherited protocol. Any proposed change to neutralization/ion conditions must be reported for scientific review before running it.

Commands: see STAGE4_RERUN_DEPLOYMENT_GUIDE.md. Production requires an explicit manual --authorize-production flag after all server gates pass.
