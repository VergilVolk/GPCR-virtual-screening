# PACER-DC G3 Scientific Summary

Date: 2026-09-26
Status: descriptive diagnostics completed; TRAINING_GATE = CLOSED

## Scope

Compound110; four contexts: CA, A, C and apo. Two MD replicas (R2/R3), five non-overlapping 1 ns windows per replica. Analyses use the frozen G2 ACh pocket (12 residues), L1 representations and saved receptor coordinates.

## Analysis definitions

dPAM = z_CA - z_A; dAGO = z_C - z_apo; dINT = dPAM - dAGO. These are operational representation-space differences, not pharmacological validation.

## Main findings

1. L1 ACh pocket dINT cross-replica cosine is positive in all five matched windows (mean 0.548038). W3 shows strong dPAM/dAGO bias alignment (cosine 0.915128) and cancellation (ratio 0.322876). This does not establish convergence.

2. Across all ten window pairs, eight ACh pocket dINT shift cosines are positive (mean 0.433684). The W2-to-W4 shift cosine is 0.981211, but its magnitude differs substantially between R2 and R3.

3. In the CA context, ACh pocket mean RMSD increases from 0.526261 to 0.693423 A in R2, and from 0.587286 to 0.609036 A in R3. RMSD is calculated after whole-receptor CA alignment using each trajectory own W2 first frame.

4. The A203-W413 CA distance increases by 0.7094 A in R2 but decreases by 0.6785 A in R3 between W2 and W4. This residue pair was selected after examining 66 pairs.

5. Four-context geometric interaction based on 66 CA distances has a W2-to-W4 cross-replica cosine of 0.295906, versus 0.981211 for L1 dINT shifts. These cosines describe different feature spaces.

6. Across ten window pairs, the descriptive Pearson correlation between representation and geometry consistency is 0.518315. Leave-one-window results range from -0.0459 to 0.8160, demonstrating sensitivity to window selection.

## Limitations

Only one compound, two MD replicas and five temporally related windows per replica were studied. Window pairs are not independent observations. Structural and representation-space associations are descriptive, not causal. No statistical significance, convergence, PAM activity or pharmacological mechanism has been established. The contribution of insufficient MD sampling remains unresolved.

## Independent validation plan

Use independently sampled longer MD trajectories to examine cumulative-time stability, replica reproducibility, conformational-state occupancy and predefined ACh pocket structural features. Preserve trajectory autocorrelation and avoid selecting validation features based on favorable new results. Keep TRAINING_GATE = CLOSED until separate training authorization and validation requirements are satisfied.
