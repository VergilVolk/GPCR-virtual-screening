# -*- coding: utf-8 -*-
"""Fail-fast integrity checks for the PACER-M4 v1 release."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import pandas as pd
from rdkit import Chem
P=Path(__file__).resolve().parents[1];R=P/"results";checks={}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 pot=pd.read_csv(P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv");checks["potency_n"]=len(pot);assert len(pot)==430
 assert (pot.source_component.value_counts()>=8).sum()==12;checks["evaluable_series"]=12
 pre=pd.read_csv(R/"pacer_candidates_v01"/"predock_portfolio.csv");dock=pd.read_csv(R/"pacer_candidates_v01"/"7trs_docking"/"features.csv");final=pd.read_csv(R/"pacer_candidates_v01"/"final"/"final_candidate_hypotheses.csv")
 assert len(pre)==200 and pre.candidate_id.nunique()==200;assert len(dock)==200 and dock.error.fillna("").eq("").all();assert len(final)==24 and final.candidate_id.nunique()==24
 known={Chem.MolToSmiles(Chem.MolFromSmiles(x),isomericSmiles=True) for x in pot.canonical_smiles};assert not (set(final.canonical_smiles)&known)
 assert final.PAINS.eq(0).all() and final.druglike.eq(1).all();assert final.canonical_smiles.map(lambda x:Chem.MolFromSmiles(x) is not None).all();assert final.potency_reference_lower90.notna().all() and final.potency_reference_upper90.notna().all();assert final.murcko_scaffold.nunique()==20
 assert set(final.portfolio_role)=={"local_exploitation","local_diversification","exploratory_hypothesis"};assert final.pocket_residue_coverage.between(0,1).all()
 sdf=list(Chem.SDMolSupplier(str(R/"pacer_candidates_v01"/"final"/"final_candidate_hypotheses.sdf"),removeHs=False));assert len(sdf)==24 and all(x is not None for x in sdf)
 checks.update({"predock_n":200,"docking_failures":0,"final_n":24,"unique_murcko_scaffolds":20,"known_exact_overlap":0,"PAINS":0,"roles":final.portfolio_role.value_counts().to_dict()})

 # Thompson-Miao exact-score reproduction and PACER-XR frozen benchmark.
 exact=pd.read_csv(R/"m4_official_exact_reproduction"/"m4_table4_exact_reproduction.csv");assert len(exact)==6
 xr=pd.read_csv(R/"pacer_xr_v01"/"m4_blind_test_metrics.csv");audit=json.loads((R/"pacer_xr_v01"/"audit.json").read_text())
 assert audit["m4_labels_used_for_training_or_selection"] is False and audit["m4_intersection_n"]==114321
 assert audit["PACER_XR_AUC"]>audit["best_single_AUC"] and audit["delta_AUC_CI95"][0]>0
 cascade=audit["PACER_XR_Cascade1pct"]
 assert cascade["ROC_AUC"]>audit["PACER_XR_AUC"] and cascade["delta_AUC_CI95"][0]>0
 assert cascade["EF_0.5pct"]>20 and cascade["EF_1pct"]>13
 semantic=json.loads((R/"pacer_xr_functional_semantics_v01"/"audit.json").read_text())
 semantic_metrics=pd.read_csv(R/"pacer_xr_functional_semantics_v01"/"apparent_metrics_nonpromotable.csv")
 assert semantic["analysis_role"]=="semantic_stress_test_only"
 assert semantic["mapped_unique_molecules"]==299 and semantic["A_tier_confirmed_negative"]==0
 assert semantic["B_tier_database_text_negative"]==33 and semantic["negative_also_ASD_active"]==33
 assert not semantic_metrics.promotable_functional_metric.astype(bool).any()
 transfer=pd.read_csv(R/"pacer_xr_candidate_transfer_v01"/"candidate_six_channel_raw_scores.csv")
 raw_cols=[x for x in transfer.columns if x.endswith("_raw")]
 assert len(transfer)==24 and len(raw_cols)==6 and transfer[raw_cols].isna().all().all()
 deltar=json.loads((R/"pacer_deltar_v01"/"audit.json").read_text())
 meta=json.loads((R/"pacer_meta_metric_v01"/"audit.json").read_text())
 residue=json.loads((R/"m4_residue_potency_association_v01"/"audit.json").read_text())
 mechanism=json.loads((R/"pacer_candidate_7trq_mechanism_v01"/"audit.json").read_text())
 assert deltar["preregistered_threshold"]["delta_ci_lower_gt_zero"] is False
 assert meta["preregistered_pass"] is False and meta["PACER_Hybrid_macro_Spearman"]>0.263
 assert residue["significant_consistent_features"]==3
 assert residue["physchem_adjusted_significant_features"]==0
 assert mechanism["n_candidates"]==24 and mechanism["docking_failures"]==0
 assert mechanism["all_three_geometry_checks"]==0 and mechanism["use_as_ranking_score"] is False
 assert (R/"pacer_biology_v01"/"PACER_M4_BIOLOGY_SUMMARY.png").exists()
 match=pd.read_csv(R/"m4_gamd_public_recluster_v01"/"representative_hungarian_match.csv")
 cluster=json.loads((R/"m4_gamd_public_recluster_v01"/"audit.json").read_text())
 assert len(match)==10 and match.pocket_RMSD_A.lt(2).all() and cluster["public_frames"]==3000
 batch=pd.read_csv(R/"pacer_assay_closed_loop_v01"/"round1_assay_batch.csv")
 anchors=set(batch.loc[batch.assay_role=="calibration_anchor","candidate_id"])
 assert len(batch)==12 and anchors=={"PACER0024","PACER0153","PACER0039"}
 template=pd.read_csv(R/"pacer_assay_closed_loop_v01"/"round1_results_template.csv")
 required_mechanism={"assay_probe","assay_readout","functional_pKB","log_alpha_beta","log_tauB","operational_model_qc"}
 assert required_mechanism.issubset(template.columns)
 # Independent Suven external set and post-hoc multidomain development.
 suven=pd.read_csv(P/"data"/"benchmarks"/"m4_pam_v1"/"external_suven_2025"/"external_patent_potency.csv")
 suven_audit=json.loads((R/"pacer_external_suven_v01"/"audit.json").read_text())
 assert len(suven)==45 and suven.smiles.map(lambda x:Chem.MolFromSmiles(x) is not None).all()
 assert suven.exact_training_overlap.sum()==0 and suven_audit["n_exact_training_overlap"]==0
 assert suven_audit["external_labels_used_for_training_or_model_selection"] is False
 assert suven_audit["preregistered_external_support"] is False
 assert suven_audit["primary_absolute_spearman"]["ci95"][0] < 0
 assert suven_audit["absolute_vs_similarity_knn"]["ci95"][0] > 0
 multidomain=json.loads((R/"pacer_multidomain_sar_v01"/"audit.json").read_text())
 assert multidomain["all_hyperparameters_selected_inside_outer_training"] is True
 assert multidomain["independent_external_claim"] is False
 sequential=pd.read_csv(R/"pacer_multidomain_vu_sequential_v01"/"metrics.csv")
 pooled=sequential.loc[sequential.method=="Pooled-Ridge"].iloc[0]
 assert pooled.ceiling_Spearman>0.94 and pooled.exact_permutation_p_two_sided<0.01
 stereo=json.loads((R/"pacer_generation4_stereo_v01"/"audit.json").read_text())
 assert stereo["n_primary_stereoisomers"]==4 and stereo["primary_method_GaMD_BEavg_spearman"]==0.0
 monash=json.loads((R/"pacer_monash_allostery_v01"/"audit.json").read_text())
 monash_data=json.loads((P/"data"/"benchmarks"/"m4_pam_v1"/"external_monash_ly2033298"/"audit.json").read_text())
 assert monash_data["n_molecules"]==16 and monash_data["n_complete_functional_triplet"]==13
 assert monash_data["n_exact_history_overlap"]==1 and monash["mechanism_support"] is False
 mech=json.loads((R/"pacer_mechanistic_multidomain_v01"/"audit.json").read_text())
 assert mech["preregistered_support"] is False and mech["paired_bootstrap_delta"]["ci95"][0]<0
 robust=pd.read_csv(R/"pacer_assay_anchor_robustness_v01"/"random_anchor_summary.csv")
 ab=robust[(robust.endpoint=="cAMP_log_alpha_beta")&(robust.method=="ThreeShot-Offset")].iloc[0]
 tau=robust[(robust.endpoint=="cAMP_log_tauB")&(robust.method=="ThreeShot-Offset")].iloc[0]
 assert ab.q025_Spearman>0.5 and tau.q025_Spearman>0
 # Acadia 2025 untouched campaign: active-moiety deduplication, external failure,
 # anchor sensitivity, and multi-objective functional biology.
 acadia_dir=P/"data"/"benchmarks"/"m4_pam_v1"/"external_acadia_2025"
 acadia=pd.read_csv(acadia_dir/"external_molecule_benchmark.csv")
 acadia_build=json.loads((acadia_dir/"external_benchmark_audit.json").read_text())
 acadia_audit=json.loads((R/"pacer_external_acadia_2025_v01"/"audit.json").read_text())
 acadia_robust=json.loads((R/"pacer_external_acadia_anchor_robustness_v01"/"audit.json").read_text())
 acadia_biology=json.loads((R/"acadia_pam_agonism_tradeoff_v01"/"summary.json").read_text())
 assert len(acadia)==45 and acadia.eligible_exact_potency.sum()==43
 assert acadia.canonical_smiles.nunique()==45 and acadia_build["n_duplicate_active_moiety_groups"]==2
 assert acadia.exact_training_overlap.sum()==0
 assert acadia_audit["n_unique_active_moieties"]==45 and acadia_audit["n_exact_potency"]==43
 assert acadia_audit["primary_external_support"] is False
 assert acadia_audit["primary_vs_similarity_knn_spearman"]["estimate"]<0
 assert acadia_robust["n_random_label_blind_triples"]==1000
 assert acadia_robust["primary_fraction_rank_better_than_knn"]<0.65
 assert acadia_biology["correlations"]["pEC50_vs_intrinsic_agonist_RE"]["ci95"][0]>0
 assert acadia_biology["correlations"]["PAM_RE_vs_intrinsic_agonist_RE"]["ci95"][0]>0
 # Function-Space v3 development failure and fifth-campaign ordinal stress.
 adaptive=json.loads((R/"pacer_adaptive_anchors_v01"/"audit.json").read_text())
 us_data=pd.read_csv(P/"data"/"benchmarks"/"m4_pam_v1"/"external_us20260055116"/"named_structure_functional_subset.csv")
 us_unique=pd.read_csv(R/"pacer_external_us20260055116_v01"/"external_unique_active_moieties.csv")
 us_audit=json.loads((R/"pacer_external_us20260055116_v01"/"audit.json").read_text())
 us_div=json.loads((R/"us20260055116_functional_divergence_v01"/"summary.json").read_text())
 context=json.loads((R/"pacer_context_kernel_v01"/"audit.json").read_text())
 function_panel=json.loads((R/"pacer_candidate_function_space_v01"/"audit.json").read_text())
 assert adaptive["preregistered_development_pass"] is False
 assert len(us_data)==26 and len(us_unique)==25 and us_unique.canonical_smiles.nunique()==25
 assert not us_unique.exact_train_overlap.astype(bool).any()
 assert us_audit["labels_visible_during_source_qualification"] is True
 assert us_audit["model_or_feature_changed_after_label_view"] is False
 assert us_audit["primary_support_against_every_baseline"] is False
 assert us_audit["n_query_after_7_anchors"]==18
 assert us_audit["cross_readout_audit"]["rat_m4_perk"]["exact_class_agreement"]<0.20
 assert us_div["human_pERK_vs_rat_pERK"]["human_stronger_bins"]==21
 assert us_div["human_pERK_vs_rat_pERK"]["two_sided_sign_test_p_non_ties"]<0.001
 assert us_div["human_pERK_vs_human_GTPgammaS"]["human_stronger_bins"]==6
 assert us_div["n_local_pairs_with_ge2_bin_endpoint_cliff"]>=20
 assert context["all_hyperparameters_selected_inside_outer_training"] is True
 assert context["development_pass"] is False
 assert context["primary_macro_ordinal_concordance"]["PACER-ContextKRR"] < context["primary_macro_ordinal_concordance"]["PooledNoContext-KRR"]
 assert function_panel["n_candidates"]==8
 assert function_panel["transferred_priors_are_pam_labels"] is False
 assert function_panel["n_candidates_with_local_functional_prior"]==0
 assert (R/"pacer_function_space_v3_figure"/"PACER_FUNCTION_SPACE_V3_EVIDENCE.png").exists()
 checks.update({"official_methods_reproduced":6,"pacer_xr_auc":audit["PACER_XR_AUC"],
                "pacer_xr_delta_ci95":audit["delta_AUC_CI95"],"public_gamd_frames":3000,
                "pacer_xr_cascade_auc":cascade["ROC_AUC"],
                "pacer_xr_cascade_delta_ci95":cascade["delta_AUC_CI95"],
                "semantic_overlap_n":semantic["mapped_unique_molecules"],
                "semantic_A_tier_negatives":semantic["A_tier_confirmed_negative"],
                "semantic_metric_promotable":False,
                "candidate_xr_contract_n":len(transfer),"candidate_xr_channels":len(raw_cols),
                "candidate_xr_fabricated_scores":0,
                "pacer_deltar_preregistered_pass":False,
                "pacer_hybrid_macro_spearman":meta["PACER_Hybrid_macro_Spearman"],
                "pacer_hybrid_preregistered_pass":False,
                "residue_primary_signals":residue["significant_consistent_features"],
                "residue_physchem_adjusted_signals":residue["physchem_adjusted_significant_features"],
                "candidate_all_three_geometry_checks":mechanism["all_three_geometry_checks"],
                "public_cluster_matches_below_2A":int(match.pocket_RMSD_A.lt(2).sum()),
                "suven_external_n":len(suven),
                "suven_external_preregistered_support":False,
                "suven_absolute_spearman":suven_audit["primary_absolute_spearman"]["estimate"],
                "vu_sequential_pooled_spearman":float(pooled.ceiling_Spearman),
                "generation4_gamd_spearman":stereo["primary_method_GaMD_BEavg_spearman"],
                "monash_mechanism_support":False,
                "mechanistic_cross_publication_support":False,
                "cAMP_alpha_beta_three_anchor_median_spearman":float(ab.median_Spearman),
                "cAMP_tauB_three_anchor_median_spearman":float(tau.median_Spearman),
                "mechanism_template_fields_present":True,
                "acadia_unique_active_moieties":len(acadia),
                "acadia_exact_potency":int(acadia.eligible_exact_potency.sum()),
                "acadia_primary_external_support":False,
                "acadia_anchor_robustness_rank_win_fraction":acadia_robust["primary_fraction_rank_better_than_knn"],
                "acadia_potency_agonism_rho":acadia_biology["correlations"]["pEC50_vs_intrinsic_agonist_RE"]["rho"],
                "adaptive_anchor_preregistered_pass":False,
                "us20260055116_unique_named_molecules":len(us_unique),
                "us20260055116_primary_support":False,
                "us20260055116_human_rat_exact_class_agreement":us_audit["cross_readout_audit"]["rat_m4_perk"]["exact_class_agreement"],
                "us20260055116_human_stronger_than_rat_non_ties":"21/22",
                "us20260055116_local_ge2bin_cliffs":us_div["n_local_pairs_with_ge2_bin_endpoint_cliff"],
                "context_kernel_development_pass":False,
                "function_space_candidates":function_panel["n_candidates"],
                "function_space_confirmed_pams":0,
                "function_space_evidence_figure":True,
                "round1_candidates":int((batch.is_control==0).sum()),"round1_controls":int((batch.is_control==1).sum())})
 files=[P/"PACER_M4_FINAL_REPORT.md",P/"docs"/"PACER_M4_RESEARCH_CHARTER.md",
        P/"docs"/"THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md",P/"docs"/"PACER_M4_LEAD_PAM_HYPOTHESES.md",
        P/"docs"/"PACER_XR_CASCADE_METHOD.md",
        R/"pacer_xr_functional_semantics_v01"/"VALIDATION_REPORT.md",
        R/"pacer_candidates_v01"/"final"/"final_candidate_hypotheses.csv"]
 checks["sha256"]={str(x.relative_to(P)):sha(x) for x in files};out=R/"pacer_release_v1_integrity.json";out.write_text(json.dumps(checks,indent=2),encoding="utf-8");print(json.dumps(checks,indent=2))
if __name__=="__main__":main()
