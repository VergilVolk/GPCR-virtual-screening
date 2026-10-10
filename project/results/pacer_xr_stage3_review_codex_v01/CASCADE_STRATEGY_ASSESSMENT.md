# Cascade strategy assessment

## Tested sensitivity

Frozen Cascade was exactly reproduced before the supplied read-only ablations. The official Glide-BEmin top-1% threshold is -9.12206074744592; this review does not change it.

|Scenario|Spearman vs Cascade|Kendall|Top10 shared|Top20 shared|Mean absolute movement|
|---|--:|--:|---:|---:|--:|
|cascade_v01|1.000|1.000|10/10|20/20|0.0|
|pure_six_mean|0.971|0.957|3/10|16/20|4.0|
|glide_only|0.777|0.583|2/10|11/20|29.5|
|vina_only|0.788|0.608|2/10|11/20|26.9|
|pdb_only|0.786|0.618|2/10|13/20|25.8|
|ensemble_only|0.896|0.732|5/10|11/20|19.3|
|bemin_only|0.846|0.664|4/10|9/20|23.6|
|beavg_only|0.881|0.707|3/10|11/20|20.3|
|loo_Glide_PDB|0.923|0.778|4/10|14/20|16.3|
|loo_Glide_BEmin|0.927|0.815|2/10|16/20|12.9|
|loo_Glide_BEavg|0.961|0.900|3/10|16/20|7.2|
|loo_Vina_PDB|0.940|0.818|3/10|12/20|13.1|
|loo_Vina_BEmin|0.959|0.866|3/10|16/20|9.5|
|loo_Vina_BEavg|0.963|0.890|3/10|16/20|8.0|
|cascade_gate_0.5pct|0.986|0.978|5/10|18/20|2.1|
|cascade_gate_2pct|0.976|0.948|9/10|11/20|4.9|
|cascade_gate_5pct|0.898|0.820|9/10|11/20|15.4|

## Assessment

- Pure six-channel mean remains globally similar to Cascade (Spearman 0.971, Kendall 0.957; Top20 overlap 16/20) but shares only 3/10 Top10 candidates. The operational priority end is gate-sensitive.
- Removing the gate changes PACERGEN02015 4->122, PACERGEN00462 3->104, PACERGEN00123 2->100, and PACERGEN02355 7->66. These movements show that their Cascade positions are principally gate membership rather than six-channel-mean support.
- Gate sensitivity is material: the 0.5% gate shares 5/10 Top10 entries, while the 5% gate has mean absolute movement 15.4. The frozen gate remains unchanged.
- Full200 matched engine correlations are modest: G-PDB/V-PDB 0.213, G-BEmin/V-BEmin 0.257, and G-BEavg/V-BEavg 0.462. Within-engine ensemble correlations are higher: Glide BEmin/BEavg 0.736 and Vina BEmin/BEavg 0.896.
- Cascade is traceable and defensible as a pre-existing hypothesis prioritizing exceptional Glide-BEmin ensemble performance. It is not structurally condition-robust enough to claim superiority without independent pose and functional validation. The balanced pure-XR Top3 are useful comparison cases, not a mandate to replace the frozen protocol.
