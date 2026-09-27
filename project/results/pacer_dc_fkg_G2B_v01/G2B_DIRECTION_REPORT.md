# PACER-FKG G2-B common-kernel direction audit

Status: descriptive only. No gate modification, drug efficacy inference, or independent-replica significance.

Frozen input checks: 40/40; RFF dimension: 512; bootstrap: 400.

| Region | Axis | pooled R2/R3 cosine | block-bootstrap 2.5–97.5% | calibration RFF MAE |
|---|---|---:|---:|---:|
| compound110_extension | synergy_interaction | -0.171 | -0.321 to 0.024 | 0.0265 |
| compound110_extension | intrinsic_agonism | -0.021 | -0.204 to 0.187 | 0.0265 |
| compound110_extension | conditional_pam_effect | 0.122 | -0.101 to 0.288 | 0.0265 |
| cooperativity_mutagenesis | synergy_interaction | 0.011 | -0.163 to 0.255 | 0.0269 |
| cooperativity_mutagenesis | intrinsic_agonism | 0.357 | 0.204 to 0.435 | 0.0269 |
| cooperativity_mutagenesis | conditional_pam_effect | 0.442 | 0.314 to 0.528 | 0.0269 |
| distal_control | synergy_interaction | -0.604 | -0.701 to -0.213 | 0.0285 |
| distal_control | intrinsic_agonism | -0.810 | -0.810 to -0.538 | 0.0285 |
| distal_control | conditional_pam_effect | 0.134 | -0.069 to 0.326 | 0.0285 |
| intracellular_microswitches | synergy_interaction | 0.015 | -0.181 to 0.205 | 0.0258 |
| intracellular_microswitches | intrinsic_agonism | -0.308 | -0.544 to -0.014 | 0.0258 |
| intracellular_microswitches | conditional_pam_effect | 0.110 | -0.069 to 0.245 | 0.0258 |
| orthosteric_activation_core | synergy_interaction | 0.467 | 0.083 to 0.607 | 0.0245 |
| orthosteric_activation_core | intrinsic_agonism | 0.179 | -0.082 to 0.331 | 0.0245 |
| orthosteric_activation_core | conditional_pam_effect | 0.334 | 0.051 to 0.468 | 0.0245 |
| orthosteric_contact_union | synergy_interaction | 0.423 | 0.003 to 0.584 | 0.0295 |
| orthosteric_contact_union | intrinsic_agonism | 0.184 | -0.171 to 0.396 | 0.0295 |
| orthosteric_contact_union | conditional_pam_effect | 0.433 | 0.112 to 0.543 | 0.0295 |
| pam_contact_consensus | synergy_interaction | -0.135 | -0.252 to 0.078 | 0.0281 |
| pam_contact_consensus | intrinsic_agonism | -0.168 | -0.315 to 0.074 | 0.0281 |
| pam_contact_consensus | conditional_pam_effect | 0.514 | 0.397 to 0.547 | 0.0281 |
| pam_contact_union | synergy_interaction | -0.228 | -0.326 to -0.024 | 0.0266 |
| pam_contact_union | intrinsic_agonism | -0.241 | -0.385 to -0.002 | 0.0266 |
| pam_contact_union | conditional_pam_effect | 0.414 | 0.278 to 0.456 | 0.0266 |
| stable_core_control | synergy_interaction | 0.238 | 0.062 to 0.356 | 0.0272 |
| stable_core_control | intrinsic_agonism | 0.167 | -0.019 to 0.290 | 0.0272 |
| stable_core_control | conditional_pam_effect | 0.025 | -0.148 to 0.181 | 0.0272 |

## Interpretation boundary

- Random Fourier features approximate the RBF RKHS; calibration MAE is in-sample.
- R2 is calibration AND evaluated replica; R3 is held out for preprocessing and bandwidth, but not an independent compound.
- 5 contiguous windows and 20-frame blocks are serially correlated; bootstrap is descriptive, not a population CI or p-value.
- A positive cross-replica cosine is RKHS contrast direction agreement, not pharmacological sign.
- Regions overlap; do not count passing regions as independent evidence.
- Shared kernel/preprocessing differs from original per-window estimator; original gate remains unchanged.
- Two replicas alone cannot establish generalization or efficacy.
