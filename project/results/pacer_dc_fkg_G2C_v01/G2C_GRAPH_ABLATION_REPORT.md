# PACER-FKG G2-C signed graph ablation

Descriptive only; no gate update, efficacy or causal graph claim.

Inputs verified: 40/40; shared-node RFF dimension: 256; shuffled graphs: 12.

| Region | Axis | no graph | frozen graph | null median | true-minus-no | descriptive delta block interval |
|---|---|---:|---:|---:|---:|---:|
| compound110_extension | synergy_interaction | -0.109 | 0.035 | -0.182 | 0.144 | 0.003 to 0.209 |
| cooperativity_mutagenesis | synergy_interaction | -0.143 | -0.165 | -0.186 | -0.022 | -0.156 to 0.125 |
| distal_control | synergy_interaction | -0.690 | -0.691 | -0.669 | -0.001 | -0.031 to 0.046 |
| intracellular_microswitches | synergy_interaction | -0.124 | 0.032 | -0.142 | 0.156 | -0.053 to 0.275 |
| orthosteric_activation_core | synergy_interaction | -0.013 | -0.086 | -0.012 | -0.074 | -0.163 to 0.066 |
| orthosteric_contact_union | synergy_interaction | 0.063 | -0.018 | 0.044 | -0.080 | -0.201 to 0.076 |
| pam_contact_consensus | synergy_interaction | -0.133 | -0.260 | -0.170 | -0.128 | -0.202 to 0.027 |
| pam_contact_union | synergy_interaction | -0.134 | -0.268 | -0.156 | -0.135 | -0.233 to 0.061 |
| stable_core_control | synergy_interaction | 0.359 | 0.108 | 0.257 | -0.251 | -0.313 to -0.159 |
| compound110_extension | intrinsic_agonism | 0.321 | 0.260 | 0.398 | -0.061 | -0.191 to 0.079 |
| cooperativity_mutagenesis | intrinsic_agonism | 0.277 | 0.141 | 0.344 | -0.136 | -0.237 to -0.011 |
| distal_control | intrinsic_agonism | -0.774 | -0.665 | -0.573 | 0.109 | 0.047 to 0.279 |
| intracellular_microswitches | intrinsic_agonism | 0.053 | 0.388 | 0.141 | 0.335 | 0.167 to 0.417 |
| orthosteric_activation_core | intrinsic_agonism | 0.133 | 0.150 | 0.273 | 0.017 | -0.089 to 0.099 |
| orthosteric_contact_union | intrinsic_agonism | 0.100 | 0.162 | 0.266 | 0.062 | -0.022 to 0.123 |
| pam_contact_consensus | intrinsic_agonism | 0.118 | 0.129 | 0.154 | 0.011 | -0.118 to 0.096 |
| pam_contact_union | intrinsic_agonism | 0.083 | 0.113 | 0.139 | 0.031 | -0.103 to 0.147 |
| stable_core_control | intrinsic_agonism | 0.134 | -0.122 | 0.160 | -0.256 | -0.334 to -0.148 |
| compound110_extension | conditional_pam_effect | 0.176 | 0.290 | 0.110 | 0.114 | -0.007 to 0.221 |
| cooperativity_mutagenesis | conditional_pam_effect | -0.022 | -0.200 | -0.041 | -0.178 | -0.290 to -0.018 |
| distal_control | conditional_pam_effect | 0.034 | -0.092 | 0.002 | -0.126 | -0.171 to -0.018 |
| intracellular_microswitches | conditional_pam_effect | -0.205 | -0.192 | -0.151 | 0.013 | -0.101 to 0.106 |
| orthosteric_activation_core | conditional_pam_effect | 0.298 | 0.201 | 0.208 | -0.098 | -0.187 to 0.013 |
| orthosteric_contact_union | conditional_pam_effect | 0.300 | 0.209 | 0.191 | -0.091 | -0.181 to 0.005 |
| pam_contact_consensus | conditional_pam_effect | 0.262 | -0.056 | 0.223 | -0.318 | -0.389 to -0.191 |
| pam_contact_union | conditional_pam_effect | 0.315 | -0.028 | 0.281 | -0.343 | -0.417 to -0.129 |
| stable_core_control | conditional_pam_effect | 0.393 | 0.355 | 0.291 | -0.038 | -0.188 to 0.088 |

## Interpretation boundaries

- New shared 128D node kernel, NOT the per-region G2-B RFF estimand; across-stage cosine magnitudes must not be compared.
- Signed RFF node deltas are graph-diffused; not the nonnegative FKG sqrt-U2 scores.
- Contact edges rewired with unweighted degree preserved and backbone fixed; node weighted strengths are NOT preserved.
- Same R2-fitted channel median/MAD as frozen G2-B; shared node bandwidth fitted on R2 only.
- Null controls are descriptive conditional on one frozen graph; no p-values or randomization inference.
- Only two replica trajectories, five serially correlated windows and overlapping regions; resampling is sensitivity only.
- G2-C does not certify PAM efficacy or graph causal mechanism; G2-B remains frozen.
