# -*- coding: utf-8 -*-
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
P=Path(__file__).resolve().parents[1];O=P/"results"/"pacer_m4r_route_audit";O.mkdir(parents=True,exist_ok=True)
def main():
 full=json.loads((P/"results"/"pacer_m4r_enrichment_v01"/"metrics.json").read_text());ext=json.loads((P/"results"/"pacer_m4r_external_pu_v01"/"metrics.json").read_text())["methods"]
 fig,ax=plt.subplots(1,2,figsize=(10,4.2));names=["Vina","Descriptors","ECFP","LambdaRank"];keys=["PDB_Vina","DescriptorClassifier","ECFP_Classifier","ECFP_LambdaRank"];ef=[full[x]["aggregate"]["macro_EF0.5%"] for x in keys];ax[0].bar(names,ef,color=["#777","#d98841","#4776a8","#4d9b67"]);ax[0].axhline(1,color="black",ls="--",lw=1,label="random");ax[0].set_ylabel("Scaffold-held-out EF(0.5%)");ax[0].set_title("Artificial-decoy benchmark");ax[0].tick_params(axis="x",rotation=25);ax[0].legend(frameon=False)
 names2=["Similarity","Descriptor PU","ECFP PU","DUD-E supervised"];keys2=["PositiveSimilarity","Descriptor_PU","ECFP_PU","DUD_E_Supervised"];pooled=[ext[x]["pooled_ROC_AUC"] for x in keys2];macro=[ext[x]["macro_group_AUC"] for x in keys2];x=np.arange(4);ax[1].bar(x-.18,pooled,.36,label="Pooled AUC",color="#9aa5b1");ax[1].bar(x+.18,macro,.36,label="Macro series AUC",color="#4c78a8");ax[1].axhline(.5,color="black",ls="--",lw=1);ax[1].set_xticks(x,names2,rotation=25,ha="right");ax[1].set_ylim(0,1);ax[1].set_title("External experimental inactive test");ax[1].legend(frameon=False);fig.tight_layout();fig.savefig(O/"decoy_shortcut_audit.png",dpi=220,bbox_inches="tight");fig.savefig(O/"decoy_shortcut_audit.svg",bbox_inches="tight")
if __name__=="__main__":main()
