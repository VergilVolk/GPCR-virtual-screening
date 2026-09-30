# PACER-DC 20 ns Retrospective Closure — Final Report v01

## 1. Scope

This experiment tested whether the frozen PACER-DC / PACER-FKG analysis stack could recover a mechanistically interpretable signal for LY2119620, a known M4 positive allosteric modulator, without outcome-driven retraining, recalibration, threshold selection, or modification of the historical reference analysis.

The closure used three paired replicas and a matched 20 ns observation window:

- 400 frames per trajectory;
- 50 ps frame spacing;
- 20-frame / 1 ns non-overlapping blocks;
- 20 blocks per trajectory;
- paired random seeds R1 = 27101, R2 = 38201, R3 = 49301.

Historical apo, probe-only, and compound110 trajectories were truncated to their first 20 ns. LY2119620 was simulated directly for 20 ns in candidate-only and candidate+probe contexts.

The three frozen factorial contrasts were:

- ΔPAM = candidate+probe − probe;
- ΔAGO = candidate − apo;
- ΔINT = candidate+probe − probe − candidate + apo.

No block or frame was treated as an independent biological replicate. The independent replication unit remained the simulation seed.

## 2. Track B — frozen C1-BS256 + PACER-FKG-v02

Track B completed successfully and is the primary executed retrospective closure.

All six LY trajectories passed the frozen C1-BS256 extraction contract. The resulting 400 × 270 × 256 representations were divided into twenty 1 ns blocks and mapped through the previously frozen PACER-FKG-v02 STATE_MOTION and SIGNED_DRIFT branches.

No Phase-2 recalibration was performed. Frozen normalization, bandwidths, RFF coordinates, graph diffusion and region definitions were retained.

Historical PACER-FKG-v02 freeze:

`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`

Phase-1 closure freeze receipt:

`58b0c76459683018b1c850dfe12fd1b4b616a7fe79727aa818d2aba0d26070ad`

Phase-2A frozen-apply receipt:

`e42d3445b3374a3837a6c12d6b1c5e2e97b54824fd2747d08902394faad6744f`

Matched 20 ns graph/region result:

`0ad0c2a0be008b38ccf8f5b0551aa6dab5552e459699cb4ffe43d5287e542c0f`

## 3. Main scientific result

The strongest positive closure result is the interaction contrast ΔINT in the preregistered `compound110_extension` region.

For STATE_MOTION, LY2119620 showed an R1/R3 direction cosine of approximately 0.63. Its ΔINT vector magnitude was of the same order as the historical compound110 matched-20ns reference, and LY-versus-compound110 direction agreement was strong in R1 and R3.

The independent SIGNED_DRIFT branch provided concordant evidence. LY-versus-compound110 direction cosines were positive in all three replicas, while LY/compound110 vector-magnitude ratios remained close to unity.

Therefore the most defensible positive conclusion is:

**The frozen PACER-FKG-v02 representation detects a reproducible candidate×probe interaction-associated conformational response for LY2119620, most clearly in the same preregistered `compound110_extension` region that was the strongest historical interaction-associated signal for compound110.**

This is evidence that the method is sensitive to a relevant non-additive receptor response outside the compound110 trajectory set.

## 4. What did not validate

Direct ΔPAM direction was not replica-stable over matched 20 ns windows.

This instability was observed not only for LY2119620 but also in the matched 20 ns compound110 reference. Consequently, the current closure does not justify describing ΔPAM direction as a reliable 20 ns PAM classification endpoint.

ΔAGO also produced measurable candidate-alone conformational responses. These should not be equated with pharmacological agonism. In LY2119620, several of the more reproducible candidate-alone responses were localized to PAM-contact or allosteric regions, while intracellular microswitch direction was considerably less stable.

Accordingly, the ΔAGO signal is interpreted as candidate-induced receptor conformational perturbation rather than evidence of intrinsic agonist efficacy.

## 5. Track B conclusion

Track B is classified as:

**PARTIAL POSITIVE retrospective closure.**

What validated:

- frozen apply-only execution on a new known PAM;
- candidate×probe interaction sensitivity;
- strongest signal localized to the historical interaction-associated `compound110_extension` region;
- comparable LY and compound110 ΔINT magnitude;
- concordant interaction signal across STATE_MOTION and SIGNED_DRIFT branches.

What did not validate:

- a stable direct ΔPAM directional endpoint;
- an automatic PAM-efficacy rule;
- a PAM/non-PAM classifier;
- a new quantitative qualification threshold.

No new threshold was created from LY2119620.

## 6. Track A — historical 128D common-kernel

The historical 128D common-kernel analysis remains frozen and available as a compound110 reference.

Its frozen long-MD analysis had already identified `dINT / compound110_extension` as the strongest interaction-associated signal, with all three replicas positive relative to the stable-core control.

However, the synchronized closure repository did not contain the legacy `GEOM2VEC_LONG_MD_MANIFEST_v01.json` required to authenticate the exact upstream 128D representation used by this historical track. Only the expected manifest SHA was retained in the frozen contract.

Because the strict closure requires apply-only use of the identical upstream representation, a new LY2119620 128D extraction was not reconstructed after observing the Track B result.

Track A therefore remains a historical reference rather than a newly executed LY validation track.

## 7. Track C — PACER-MCV

PACER-MCV v0.1 provides a preregistered, encoder-free physical collective-variable panel.

The feature configuration was frozen before this closure, including receptor distances, side-chain distances, regional compactness variables, region-centroid distances, and the three factorial axes.

However, the closure execution manifest records only:

`FROZEN_CONFIG`

No pre-LY frozen R2-derived scaling artifact was materialized.

The historical MCV implementation estimates separate robust scales from R2 before application to another replica. Reconstructing that calibration after observing the LY result would violate the strict no-post-outcome-refit policy.

Track C was therefore deliberately not executed in this retrospective closure.

This omission is a provenance/calibration constraint, not a negative biological result.

## 8. Overall interpretation

The strict retrospective closure produced one completed independent computational validation track: frozen PACER-FKG-v02.

That track provides positive evidence that PACER-FKG detects an interaction-associated conformational response for a known PAM in a new MD dataset, with the strongest signal appearing in the same predefined receptor region implicated by the historical compound110 analysis.

The closure does not establish a validated PAM classifier or efficacy predictor.

The strongest supported capability is narrower:

**PACER-FKG-v02 acts as an interaction-sensitive mechanistic representation capable of detecting reproducible non-additive receptor responses across independently generated trajectories.**

The main unresolved limitation is replica instability of the direct ΔPAM directional endpoint at the 20 ns timescale.

## 9. Next preregistered test

The next experiment is the inactive hard-negative Stage B using CM00734.

The LY2119620 result is frozen before CM00734 is examined.

No encoder weights, preprocessing parameters, bandwidths, RFF coordinates, graph parameters, region definitions, thresholds, or qualification rules should be changed as a consequence of either the LY2119620 or CM00734 outcomes.

The principal Stage B question is whether the interaction-associated response observed for LY2119620 is absent, substantially attenuated, or qualitatively different for the experimental inactive control.

Until that test is complete, all conclusions remain descriptive mechanistic validation rather than functional classification.