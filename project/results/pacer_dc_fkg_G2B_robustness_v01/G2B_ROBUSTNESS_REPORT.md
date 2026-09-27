# PACER-FKG G2-B supplementary RFF/exact-kernel audit

Status: descriptive only; original G2-B gate, samples and report remain frozen.
Frozen input hashes verified: 40/40; original G2-B SHA256 `374402a7830d8c163c785cfe3899c30b5727a4375b8065fef080c98f2c5a00c5`.

| Region | Axis | full frozen 512D cosine | exact subsample cosine | 1024D max absolute error |
|---|---|---:|---:|---:|
| compound110_extension | synergy_interaction | -0.171 | -0.192 | 0.042 |
| compound110_extension | intrinsic_agonism | -0.021 | 0.037 | 0.056 |
| compound110_extension | conditional_pam_effect | 0.122 | 0.086 | 0.080 |
| cooperativity_mutagenesis | synergy_interaction | 0.011 | -0.147 | 0.051 |
| cooperativity_mutagenesis | intrinsic_agonism | 0.357 | 0.132 | 0.042 |
| cooperativity_mutagenesis | conditional_pam_effect | 0.442 | 0.357 | 0.032 |
| distal_control | synergy_interaction | -0.604 | -0.322 | 0.032 |
| distal_control | intrinsic_agonism | -0.810 | -0.629 | 0.020 |
| distal_control | conditional_pam_effect | 0.134 | 0.194 | 0.051 |
| orthosteric_activation_core | synergy_interaction | 0.467 | 0.285 | 0.035 |
| orthosteric_activation_core | intrinsic_agonism | 0.179 | 0.157 | 0.018 |
| orthosteric_activation_core | conditional_pam_effect | 0.334 | 0.053 | 0.035 |
| pam_contact_consensus | synergy_interaction | -0.135 | -0.059 | 0.064 |
| pam_contact_consensus | intrinsic_agonism | -0.168 | -0.148 | 0.068 |
| pam_contact_consensus | conditional_pam_effect | 0.514 | 0.368 | 0.056 |
| pam_contact_union | synergy_interaction | -0.228 | -0.102 | 0.053 |
| pam_contact_union | intrinsic_agonism | -0.241 | -0.152 | 0.039 |
| pam_contact_union | conditional_pam_effect | 0.414 | 0.311 | 0.048 |
| stable_core_control | synergy_interaction | 0.238 | 0.243 | 0.026 |
| stable_core_control | intrinsic_agonism | 0.167 | 0.108 | 0.046 |
| stable_core_control | conditional_pam_effect | 0.025 | -0.010 | 0.014 |

## Boundaries

- Exact kernel is empirical RBF mean-embedding cosine on the same deterministic balanced subsample, not an infinite-sample population kernel.
- Subsample estimates are NOT directly interchangeable with frozen G2-B full-500-frame-per-context estimates.
- R3 holdout refers only to frozen preprocessing/bandwidth calibration, not to molecule or chemical generalization.
- Finite 2-replica data, serially correlated frames and overlapping regions prohibit p-value or efficacy claims.
- No gate changes, no additional model fit and no random seed selected post hoc.
