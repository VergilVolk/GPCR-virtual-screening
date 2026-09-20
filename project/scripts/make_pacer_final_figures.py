# -*- coding: utf-8 -*-
"""Create publication-ready summary figures from frozen PACER-M4 outputs."""
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];O=P/"results"/"pacer_final_figures";O.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":180})
def save(fig,name):fig.tight_layout();fig.savefig(O/f"{name}.png",bbox_inches="tight");fig.savefig(O/f"{name}.svg",bbox_inches="tight");plt.close(fig)
def main():
 # Frozen series-LOSO comparisons.
 names=["RF","ChemBERTa","Chemprop","LightGBM v1","LightGBM v2","ACh-IFP","Nested fusion"]
 macro=[.170,.148,.147,.199,.207,.191,.151];worst=[np.nan,-.352,-.206,-.130,-.138,-.085,-.148]
 fig,ax=plt.subplots(figsize=(8,4));x=np.arange(len(names));ax.bar(x-.18,macro,.36,label="Macro Spearman",color="#3568b8");ax.bar(x+.18,np.nan_to_num(worst),.36,label="Worst series",color="#e47b40");ax.axhline(0,color="black",lw=.8);ax.set_xticks(x,names,rotation=25,ha="right");ax.set_ylabel("Source-LOSO Spearman");ax.legend(frameon=False);save(fig,"figure1_loso_baselines")
 # Evidence funnel.
 stages=["Generated\nunique","Eligible\nchemistry/domain","ACh-state\ndocked","Final\nhypotheses"];vals=[3271,1772,200,24]
 fig,ax=plt.subplots(figsize=(7,4));bars=ax.bar(stages,vals,color=["#9aa5b1","#6c91bf","#4c78a8","#2d5f8b"]);ax.set_yscale("log");ax.set_ylabel("Molecules (log scale)");ax.bar_label(bars,labels=[str(v) for v in vals],padding=3);save(fig,"figure2_evidence_funnel")
 # Final portfolio: bounded novelty vs orthogonal structure gate.
 d=pd.read_csv(P/"results"/"pacer_candidates_v01"/"final"/"final_candidate_hypotheses.csv");colors={"local_exploitation":"#2f6f9f","local_diversification":"#43a047","exploratory_hypothesis":"#e07b39"}
 fig,ax=plt.subplots(figsize=(7,5))
 for role,q in d.groupby("portfolio_role"):
  ax.scatter(q.max_tanimoto,q.pocket_residue_coverage,s=40+120*q.QED,c=colors[role],label=role.replace("_"," "),alpha=.85,edgecolor="white",linewidth=.5)
 ax.axvline(.55,color="gray",ls="--",lw=1);ax.axvline(.65,color="gray",ls=":",lw=1);ax.set_xlabel("Max ECFP4 similarity to known potency set");ax.set_ylabel("ACh-state pocket residue coverage");ax.legend(frameon=False,fontsize=8);save(fig,"figure3_final_portfolio")
 # Generator reality check.
 g=json.loads((P/"results"/"generated"/"generator_comparison.json").read_text());order=["Fragment","LSTM","GPT","Diffusion"];count=[g[x]["n"] for x in order];sim=[g[x].get("mean_Tanimoto_vs_known",0) for x in order]
 fig,ax=plt.subplots(figsize=(7,4));b=ax.bar(order,count,color="#707f8d");ax.set_yscale("symlog",linthresh=1);ax.set_ylabel("Valid outputs (symlog)");ax.bar_label(b,labels=[str(x) for x in count],padding=3);a=ax.twinx();a.plot(order,sim,"o-",color="#d8573c",label="Mean max similarity");a.set_ylabel("Mean max similarity to known set",color="#d8573c");a.set_ylim(0,.5);save(fig,"figure4_generator_audit")
if __name__=="__main__":main()
