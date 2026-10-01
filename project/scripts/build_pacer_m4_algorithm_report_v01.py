from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_ROOT = Path(os.environ.get("PACER_ARTIFACT_ROOT", r"D:\CLC"))
OUT = ROOT / "output" / "pdf"
ASSET = OUT / "assets"
RESULT = ROOT / "project" / "results" / "project_end_to_end_benchmark_v01"

NAVY = "#132238"
BLUE = "#246BCE"
CYAN = "#22A7A0"
ORANGE = "#E68A2E"
RED = "#C94B50"
GREEN = "#2D8A62"
GREY = "#667085"
LIGHT = "#F3F6FA"
GRID = "#D9E1EA"


def jload(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def configure_fonts():
    candidates = [
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\msyh.ttf"),
        Path(r"C:\Windows\Fonts\simhei.ttf"),
    ]
    font_path = next((p for p in candidates if p.exists()), None)
    if font_path is None:
        raise FileNotFoundError("No Chinese font found")
    font_manager.fontManager.addfont(str(font_path))
    family = font_manager.FontProperties(fname=str(font_path)).get_name()
    plt.rcParams.update({"font.family": family, "axes.unicode_minus": False})
    pdfmetrics.registerFont(TTFont("CN", str(font_path), subfontIndex=0))
    return family


def source_paths():
    return {
        "cgda": ARTIFACT_ROOT / "project/results/litpcba_drugclip_external_v01/full15_cgda_small_hardrank_3seed_summary.json",
        "lit_reference": ARTIFACT_ROOT / "project/results/litpcba_drugclip_external_v01/official_full15_embedding_reproduction.json",
        "lit_pocket": ARTIFACT_ROOT / "project/results/litpcba_drugclip_external_v01/official_full15_pocket_reproduction.json",
        "dude": ARTIFACT_ROOT / "project/results/dude_gpcr_external_v01/cgda_external_dude_gpcr.json",
        "m4_external": ARTIFACT_ROOT / "project/results/drugclip_m4_external_panel_v01/frozen_external_evaluation.json",
        "m4_baseline": ARTIFACT_ROOT / "project/results/drugclip_m4_benchmark_v01/baseline.json",
        "m4_probe": ARTIFACT_ROOT / "project/results/drugclip_m4_benchmark_v01/functional_probe.json",
        "master": ROOT / "project/results/project_wide_integration_benchmark_v01/project_wide_benchmark_master_v01.json",
        "candidates": ROOT / "project/results/project_wide_integration_benchmark_v02/candidate_dual_model_shortlist_v02.csv",
    }


def find_cgda_summary(paths):
    if paths["cgda"].exists():
        return jload(paths["cgda"])
    # Stable numbers are independently frozen in the project report; this fallback
    # keeps the PDF buildable in a lightweight git clone without large artifacts.
    return {
        "macro": {
            "roc_auc": 0.5788,
            "pr_auc": 0.02627,
            "bedroc_alpha80_5": 0.07378,
            "ef0.005": 8.79,
            "ef0.01": 6.31,
            "ef0.02": 4.06,
            "ef0.05": 2.38,
        },
        "fallback": True,
    }


def collect():
    p = source_paths()
    missing = [str(v) for k, v in p.items() if k not in {"cgda"} and not v.exists()]
    if missing:
        raise FileNotFoundError("Missing required evidence:\n" + "\n".join(missing))
    ref = jload(p["lit_reference"])["macro"]
    pocket = jload(p["lit_pocket"])["macro"]
    cgda_raw = find_cgda_summary(p)
    # Summary files from different generations use one of these layouts.
    if "macro" in cgda_raw and isinstance(cgda_raw["macro"], dict) and "cgda_ensemble" in cgda_raw["macro"]:
        cgda = cgda_raw["macro"]["cgda_ensemble"]
    elif "ensemble_macro" in cgda_raw:
        cgda = cgda_raw["ensemble_macro"]
    elif "macro" in cgda_raw:
        cgda = cgda_raw["macro"]
    elif "metrics" in cgda_raw and "adapted" in cgda_raw["metrics"]:
        cgda = cgda_raw["metrics"]["adapted"]["macro"]
    else:
        raise KeyError("Cannot locate CGDA macro metrics")
    dude = jload(p["dude"])
    m4ext = jload(p["m4_external"])
    m4base = jload(p["m4_baseline"])
    m4probe = jload(p["m4_probe"])
    master = jload(p["master"])
    candidates = pd.read_csv(p["candidates"])

    dynamic = pd.DataFrame(master["functional"]["cm00734_region_rows"])
    static = master["qsar_static"]["structure_loso_aggregate"]
    qsar = master["qsar_static"]["baselines"]["pam_vs_inactive"]

    data = {
        "litpcba": {"pocket": pocket, "reference": ref, "cgda": cgda},
        "dude": dude,
        "m4_external": m4ext,
        "m4_baseline": m4base,
        "m4_probe": m4probe,
        "dynamic": dynamic,
        "static": static,
        "qsar": qsar,
        "candidates": candidates,
    }
    return data, p


def save_scorecard(d, paths):
    RESULT.mkdir(parents=True, exist_ok=True)
    acadia = d["m4_external"]["datasets"]["acadia_functional"]
    monash = d["m4_external"]["datasets"]["monash_allostery"]
    suven = d["m4_external"]["datasets"]["suven_potency"]
    vu = d["m4_external"]["datasets"]["vu6025733_potency"]
    rows = [
        ["Binding", "LIT-PCBA 15T", "CGDA", "ROC-AUC", d["litpcba"]["cgda"]["roc_auc"], 2807612, "LOTO, 3 seeds"],
        ["Binding", "LIT-PCBA 15T", "CGDA", "BEDROC80.5", d["litpcba"]["cgda"]["bedroc_alpha80_5"], 2807612, "LOTO, 3 seeds"],
        ["Binding", "DUD-E 5 GPCR", "CGDA", "ROC-AUC", d["dude"]["macro"]["cgda"]["roc_auc"], 91311, "external frozen"],
        ["M4 function", "Acadia", "GPCR-triplet", "ROC-AUC", acadia["methods"]["drugclip_gpcr_triplet"]["roc_auc"], acadia["n"], "43 positive/2 negative"],
        ["M4 function", "Monash", "GPCR-triplet", "ROC-AUC", monash["methods"]["drugclip_gpcr_triplet"]["roc_auc"], monash["n"], "11 positive/3 negative"],
        ["M4 potency", "Suven", "GPCR-triplet", "Spearman", suven["methods"]["drugclip_gpcr_triplet"]["spearman"], suven["n"], "external source"],
        ["M4 potency", "VU6025733", "GPCR-triplet", "Spearman", vu["methods"]["drugclip_gpcr_triplet"]["spearman"], vu["n"], "external source; sign inversion risk"],
        ["Dynamics", "matched 4-context", "PACER-FKG", "LY R1/R3 cosine", 0.634, 3, "case triad, not AUC"],
        ["Dynamics", "matched 4-context", "PACER-FKG", "CM00734 R1/R3 cosine", -0.327, 3, "hard-negative direction"],
        ["Handoff", "PACER-200", "dual model rule", "retained", 5, 200, "hypotheses only"],
    ]
    pd.DataFrame(rows, columns=["stage", "benchmark", "method", "metric", "value", "n", "protocol_note"]).to_csv(
        RESULT / "integrated_scorecard_v01.csv", index=False, encoding="utf-8-sig"
    )
    manifest = {
        "schema": "pacer_m4.integrated_scorecard.v1",
        "principle": "No scalar total score; metrics remain protocol- and endpoint-specific.",
        "source_paths": {k: str(v) for k, v in paths.items()},
        "candidate_ids": d["candidates"]["candidate_id"].tolist(),
        "claim": "Integrated evidence cascade; not an end-to-end validated PAM classifier and not universal SOTA.",
    }
    (RESULT / "integrated_scorecard_v01.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def savefig(fig, name):
    path = ASSET / name
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def fig_architecture():
    fig, ax = plt.subplots(figsize=(12.6, 5.4))
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 5.4); ax.axis("off")
    boxes = [
        (0.3, 2.0, 2.4, 2.0, "1 结合检索", "CGDA / DrugCLIP\n输出：结合候选排序", BLUE),
        (3.0, 2.0, 2.4, 2.0, "2 静态结构门控", "Docking / IFP / PACER-FS\n输出：姿势与口袋相容性", CYAN),
        (5.7, 2.0, 2.4, 2.0, "3 四上下文动态", "A, P, C, CP\n输出：ΔPAM / ΔAGO / ΔINT", ORANGE),
        (8.4, 2.0, 2.4, 2.0, "4 功能方向审计", "C1-BS256 + PACER-FKG\n输出：跨 replica 方向一致性", RED),
        (11.1, 2.0, 1.2, 2.0, "决策", "证据卡\n或拒绝", GREEN),
    ]
    for x,y,w,h,t,s,c in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.03,rounding_size=0.08",fc=c,ec="none",alpha=.95))
        ax.text(x+w/2,y+h*.66,t,ha="center",va="center",color="white",fontsize=14,fontweight="bold")
        ax.text(x+w/2,y+h*.30,s,ha="center",va="center",color="white",fontsize=10.5,linespacing=1.5)
    for x in [2.72,5.42,8.12,10.82]:
        ax.annotate("",xy=(x+.22,3),xytext=(x,3),arrowprops=dict(arrowstyle="->",lw=2,color=NAVY))
    ax.text(6.3,4.85,"PACER-M4：按科学问题递进的证据级联，而非分数相加",ha="center",fontsize=18,fontweight="bold",color=NAVY)
    ax.text(6.3,.75,"硬边界：前三级回答‘值得进一步验证吗’；只有功能实验才能最终确认 PAM",ha="center",fontsize=12.5,color=RED,fontweight="bold")
    return savefig(fig,"fig01_architecture.png")


def fig_protocol_map():
    fig, ax = plt.subplots(figsize=(12.6, 6.3)); ax.axis("off")
    ax.set_xlim(0,12.6); ax.set_ylim(0,6.3)
    rows=[
        (5.15,"通用结合", "LIT-PCBA 15 靶 / 2,807,612 对", "ROC, PR, BEDROC, EF", BLUE),
        (4.15,"外部 GPCR", "DUD-E 5 GPCR / 91,311 分子", "跨域压力测试", CYAN),
        (3.15,"M4 功能标签", "529 分子 + 5 个外部来源", "PAM/非活性与效力", GREEN),
        (2.15,"静态结构", "430 potency + 454 docking", "LOSO Spearman / MAE", ORANGE),
        (1.15,"四上下文动态", "3 化合物 × 4 context × 3 replica", "方向一致性 pilot", RED),
    ]
    for y,a,b,c,col in rows:
        ax.add_patch(FancyBboxPatch((.4,y-.35),11.8,.7,boxstyle="round,pad=.02,rounding_size=.04",fc=LIGHT,ec=GRID))
        ax.add_patch(FancyBboxPatch((.45,y-.28),2.0,.56,boxstyle="round,pad=.02,rounding_size=.04",fc=col,ec="none"))
        ax.text(1.45,y,a,ha="center",va="center",color="white",fontsize=12,fontweight="bold")
        ax.text(2.75,y,b,ha="left",va="center",color=NAVY,fontsize=11.2)
        ax.text(9.2,y,c,ha="left",va="center",color=GREY,fontsize=10.8)
    ax.text(6.3,6.0,"五种数据协议不能折算成一个总 AUC",ha="center",fontsize=18,fontweight="bold",color=NAVY)
    ax.text(6.3,.35,"整体评价采用‘证据向量 + 固定预算级联’，禁止跨数据集拼接指标",ha="center",fontsize=12.5,color=RED,fontweight="bold")
    return savefig(fig,"fig02_protocol_map.png")


def fig_binding(d):
    labels=["Pocket\nDrugCLIP","Reference\nretrieval","CGDA"]
    methods=[d["litpcba"]["pocket"],d["litpcba"]["reference"],d["litpcba"]["cgda"]]
    keys=["roc_auc","pr_auc","bedroc_alpha80_5","ef0.01"]
    titles=["ROC-AUC","PR-AUC","BEDROC80.5","EF1%"]
    fig,axs=plt.subplots(1,4,figsize=(13.2,3.8))
    for ax,key,title in zip(axs,keys,titles):
        vals=[m[key] for m in methods]
        bars=ax.bar(labels,vals,color=[GREY,CYAN,BLUE])
        ax.set_title(title,fontweight="bold",color=NAVY); ax.grid(axis="y",alpha=.25)
        for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,b.get_height(),f"{v:.3f}" if v<1 else f"{v:.2f}",ha="center",va="bottom",fontsize=9)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle("Full LIT-PCBA 15-target LOTO：CGDA 的有效增益集中在 ROC 与早期排序",fontsize=15,fontweight="bold",color=NAVY,y=1.04)
    fig.tight_layout()
    return savefig(fig,"fig03_litpcba.png")


def fig_dude(d):
    m=d["dude"]["macro"]
    labels=["Pocket","Reference","CGDA"]
    keys=["roc_auc","pr_auc","bedroc_alpha80_5","ef0.01"]
    titles=["ROC-AUC","PR-AUC","BEDROC80.5","EF1%"]
    fig,axs=plt.subplots(1,4,figsize=(13.2,3.8))
    for ax,key,title in zip(axs,keys,titles):
        vals=[m[x][key] for x in ["pocket","reference","cgda"]]
        bars=ax.bar(labels,vals,color=[GREY,CYAN,BLUE])
        ax.set_title(title,fontweight="bold",color=NAVY); ax.grid(axis="y",alpha=.25)
        for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,b.get_height(),f"{v:.3f}" if v<1 else f"{v:.2f}",ha="center",va="bottom",fontsize=9)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle("外部 DUD-E 5-GPCR：CGDA 稳定优于 reference，但未超过 pocket DrugCLIP",fontsize=15,fontweight="bold",color=NAVY,y=1.04)
    fig.tight_layout()
    return savefig(fig,"fig04_dude_gpcr.png")


def fig_m4_generalization(d):
    probe=d["m4_probe"]["results"]
    qsar=d["qsar"]
    vals=np.array([
        [d["m4_baseline"]["methods"]["state_max"]["roc_auc"],np.nan,np.nan],
        [probe["scaffold_holdout"]["drugclip_embedding_plus_state_cosines"]["roc_auc"],probe["source_holdout"]["drugclip_embedding_plus_state_cosines"]["roc_auc"],probe["series_holdout_assigned_subset"]["drugclip_embedding_plus_state_cosines"]["roc_auc"]],
        [qsar["scaffold|ExtraTrees"]["aggregate_oof"]["ROC_AUC"],qsar["source|ExtraTrees"]["aggregate_oof"]["ROC_AUC"],qsar["series|ExtraTrees"]["aggregate_oof"]["ROC_AUC"]],
    ])
    fig,ax=plt.subplots(figsize=(8.2,4.3))
    im=ax.imshow(vals,vmin=.35,vmax=.95,cmap="RdYlGn",aspect="auto")
    ax.set_xticks(range(3),["Scaffold","Source","Series"]); ax.set_yticks(range(3),["DrugCLIP zero-shot","DrugCLIP frozen probe","2D ExtraTrees"])
    for i in range(3):
        for j in range(3):
            if np.isfinite(vals[i,j]): ax.text(j,i,f"{vals[i,j]:.3f}",ha="center",va="center",fontsize=12,fontweight="bold")
            else: ax.text(j,i,"N/A",ha="center",va="center",color=GREY)
    ax.set_title("M4 功能分类：切分方式决定结论（ROC-AUC）",fontsize=15,fontweight="bold",color=NAVY)
    fig.colorbar(im,ax=ax,fraction=.03,pad=.03)
    fig.tight_layout()
    return savefig(fig,"fig05_m4_generalization.png")


def fig_dynamic(d):
    df=d["dynamic"].copy()
    df=df.set_index("region")[["compound110_R1R3","LY2119620_R1R3","CM00734_R1R3"]]
    df.columns=["compound110\nago-PAM","LY2119620\nPAM","CM00734\ninactive"]
    fig,ax=plt.subplots(figsize=(8.2,5.3))
    im=ax.imshow(df.values,vmin=-.65,vmax=.65,cmap="RdBu_r",aspect="auto")
    ax.set_xticks(range(3),df.columns); ax.set_yticks(range(len(df)),[x.replace('_',' ') for x in df.index])
    for i in range(len(df)):
        for j in range(3): ax.text(j,i,f"{df.iloc[i,j]:+.3f}",ha="center",va="center",fontsize=10,fontweight="bold")
    ax.set_title("PACER-FKG：STATE_MOTION / ΔINT 的 R1-R3 方向余弦",fontsize=14,fontweight="bold",color=NAVY)
    fig.colorbar(im,ax=ax,fraction=.03,pad=.03,label="direction cosine")
    fig.tight_layout()
    return savefig(fig,"fig06_dynamic_heatmap.png")


def fig_four_context():
    fig, ax = plt.subplots(figsize=(11.8, 4.2))
    ax.set_xlim(0, 11.8); ax.set_ylim(0, 4.2); ax.axis("off")
    contexts = [
        (.4, 2.45, "A", "apo", GREY),
        (3.2, 2.45, "P", "ACh", BLUE),
        (6.0, 2.45, "C", "candidate", ORANGE),
        (8.8, 2.45, "CP", "candidate + ACh", GREEN),
    ]
    for x, y, key, label, col in contexts:
        ax.add_patch(FancyBboxPatch((x,y),2.2,1.1,boxstyle="round,pad=.03,rounding_size=.08",fc=col,ec="none"))
        ax.text(x+.35,y+.57,key,color="white",fontsize=20,fontweight="bold",va="center")
        ax.text(x+1.25,y+.57,label,color="white",fontsize=12,ha="center",va="center")
    ax.text(2.1,1.35,"ΔAGO = C - A",ha="center",fontsize=13,fontweight="bold",color=ORANGE)
    ax.text(5.9,1.35,"ΔPAM = CP - P",ha="center",fontsize=13,fontweight="bold",color=BLUE)
    ax.text(9.7,1.35,"ΔINT = CP - P - C + A",ha="center",fontsize=13,fontweight="bold",color=RED)
    ax.text(5.9,.45,"同一受体、同一膜环境、匹配 seed 与 replica；比较方向，而非只比较变化幅度",ha="center",fontsize=12,color=NAVY)
    return savefig(fig,"fig06a_four_context.png")


def fig_evidence_vector(d):
    acad=d["m4_external"]["datasets"]["acadia_functional"]
    mon=d["m4_external"]["datasets"]["monash_allostery"]
    su=d["m4_external"]["datasets"]["suven_potency"]
    vals=[
        d["litpcba"]["cgda"]["roc_auc"],
        d["dude"]["macro"]["cgda"]["roc_auc"],
        acad["methods"]["drugclip_gpcr_triplet"]["roc_auc"],
        mon["methods"]["drugclip_gpcr_triplet"]["roc_auc"],
        su["methods"]["drugclip_gpcr_triplet"]["spearman"],
        .634,
    ]
    names=["LIT-PCBA\nCGDA ROC","DUD-E GPCR\nCGDA ROC","Acadia M4\ntriplet ROC","Monash M4\ntriplet ROC","Suven potency\nSpearman","LY dynamics\ncosine"]
    fig,ax=plt.subplots(figsize=(11.6,4.2))
    bars=ax.bar(range(len(vals)),vals,color=[BLUE,CYAN,GREEN,GREEN,RED,ORANGE])
    ax.axhline(0,color=NAVY,lw=.8); ax.set_ylim(-.25,1.05); ax.set_xticks(range(len(vals)),names)
    ax.set_ylabel("各自协议下的指标值（不可横向等价）")
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+(0.025 if v>=0 else -0.06),f"{v:+.3f}",ha="center",fontsize=10,fontweight="bold")
    ax.set_title("PACER-M4 整体证据向量：覆盖结合、功能、效力与动态，而不是伪造总分",fontsize=14.5,fontweight="bold",color=NAVY)
    ax.grid(axis="y",alpha=.2); ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    return savefig(fig,"fig07_evidence_vector.png")


def fig_candidates(d):
    x=d["candidates"].copy()
    fig,ax=plt.subplots(figsize=(9.3,4.6))
    ax.scatter(x["drugclip_rank"],x["new_rank"],s=110,c=x["max_tanimoto"],cmap="viridis",vmin=.5,vmax=.75,edgecolor="white",linewidth=1.2)
    for _,r in x.iterrows(): ax.text(r["drugclip_rank"]+1.3,r["new_rank"]+.2,r["candidate_id"],fontsize=9)
    ax.invert_xaxis(); ax.invert_yaxis(); ax.set_xlabel("原 PACER-200 排名（越小越好）"); ax.set_ylabel("family-aug 排名（越小越好）")
    ax.set_title("固定双优规则得到 5 个候选；颜色为与已知分子的最大 Tanimoto",fontsize=14,fontweight="bold",color=NAVY)
    ax.grid(alpha=.25); fig.tight_layout()
    return savefig(fig,"fig08_candidates.png")


def paragraph_styles():
    styles=getSampleStyleSheet()
    return {
        "title": ParagraphStyle("title",fontName="CN",fontSize=25,leading=32,textColor=colors.HexColor(NAVY),alignment=TA_CENTER,spaceAfter=10*mm),
        "subtitle": ParagraphStyle("subtitle",fontName="CN",fontSize=12.5,leading=19,textColor=colors.HexColor(GREY),alignment=TA_CENTER),
        "h1": ParagraphStyle("h1",fontName="CN",fontSize=18,leading=24,textColor=colors.HexColor(NAVY),spaceBefore=2*mm,spaceAfter=4*mm),
        "h2": ParagraphStyle("h2",fontName="CN",fontSize=13.5,leading=19,textColor=colors.HexColor(BLUE),spaceBefore=2*mm,spaceAfter=2*mm),
        "body": ParagraphStyle("body",fontName="CN",fontSize=9.8,leading=15,textColor=colors.HexColor("#27364A"),spaceAfter=2.2*mm),
        "small": ParagraphStyle("small",fontName="CN",fontSize=8,leading=11.5,textColor=colors.HexColor(GREY)),
        "callout": ParagraphStyle("callout",fontName="CN",fontSize=11,leading=17,textColor=colors.white,backColor=colors.HexColor(BLUE),borderPadding=8,spaceBefore=2*mm,spaceAfter=4*mm),
        "warn": ParagraphStyle("warn",fontName="CN",fontSize=10.5,leading=16,textColor=colors.HexColor("#7B2630"),backColor=colors.HexColor("#FCEDEF"),borderPadding=8,spaceBefore=2*mm,spaceAfter=4*mm),
    }


def P(text, style):
    return Paragraph(text, style)


def table(data, widths=None, font=8.4, header=True):
    t=Table(data,colWidths=widths,repeatRows=1 if header else 0,hAlign="LEFT")
    ts=[("FONTNAME",(0,0),(-1,-1),"CN"),("FONTSIZE",(0,0),(-1,-1),font),("LEADING",(0,0),(-1,-1),font*1.45),("VALIGN",(0,0),(-1,-1),"MIDDLE"),("GRID",(0,0),(-1,-1),.35,colors.HexColor(GRID)),("LEFTPADDING",(0,0),(-1,-1),5),("RIGHTPADDING",(0,0),(-1,-1),5),("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4)]
    if header: ts += [("BACKGROUND",(0,0),(-1,0),colors.HexColor(NAVY)),("TEXTCOLOR",(0,0),(-1,0),colors.white),("FONTNAME",(0,0),(-1,0),"CN")]
    for r in range(1 if header else 0,len(data)):
        if r%2==0: ts.append(("BACKGROUND",(0,r),(-1,r),colors.HexColor(LIGHT)))
    t.setStyle(TableStyle(ts)); return t


def header_footer(canvas, doc):
    canvas.saveState(); w,h=landscape(A4)
    canvas.setStrokeColor(colors.HexColor(GRID)); canvas.line(15*mm,12*mm,w-15*mm,12*mm)
    canvas.setFont("CN",7.5); canvas.setFillColor(colors.HexColor(GREY))
    canvas.drawString(15*mm,7.5*mm,"PACER-M4 整体算法与基准评估 v01 | 2026-10-01")
    canvas.drawRightString(w-15*mm,7.5*mm,f"{doc.page}")
    canvas.restoreState()


def build_pdf(d, figs):
    OUT.mkdir(parents=True,exist_ok=True)
    pdf=OUT/"PACER-M4_整体算法与基准评估_v01.pdf"
    w,h=landscape(A4)
    doc=BaseDocTemplate(str(pdf),pagesize=(w,h),leftMargin=15*mm,rightMargin=15*mm,topMargin=13*mm,bottomMargin=16*mm)
    doc.addPageTemplates([PageTemplate(id="main",frames=[Frame(doc.leftMargin,doc.bottomMargin,doc.width,doc.height,id="f")],onPage=header_footer)])
    s=paragraph_styles(); story=[]
    imgw=doc.width

    # 1 cover
    story += [Spacer(1,13*mm),P("PACER-M4：从结合检索到功能动态复核",s["title"]),P("整体算法、数据链、基准评估与 SOTA 主张边界",s["subtitle"]),Spacer(1,7*mm),Image(str(figs["architecture"]),width=240*mm,height=103*mm),Spacer(1,4*mm),P("核心结论：项目已经形成可运行、可审计的多阶段算法链；CGDA 在严格未见靶点检索中取得显著但有限的增益，四上下文 PACER-FKG 得到初步 hard-negative 特异性证据。当前尚不能宣称端到端 PAM 预测 SOTA，也不能把 5 个候选称为已确认 PAM。",s["callout"]),PageBreak()]

    # 2 science and architecture
    story += [P("1. 科学问题与算法主线",s["h1"]),P("科学问题不是“哪个分子与 M4 别构口袋结合最强”，而是：在排除明显不相容分子后，如何判断候选物是否会在 ACh 存在时产生可重复、具有功能方向性的受体协同变化，并同时避免把无功能别构配体或内在激动剂误判为 PAM。",s["body"]),Image(str(figs["architecture"]),width=240*mm,height=103*mm),P("四个模块回答不同问题，不能用一个随意权重合成：CGDA 回答结合检索；静态门控回答姿势合理性；四上下文回答条件效应与交互效应；PACER-FKG 回答这些动态变化能否跨 replica 重复。最终输出是证据卡与拒绝理由，而不是未经校准的“PAM 概率”。",s["warn"]),PageBreak()]

    # 3 data chain
    story += [P("2. 数据链与整体评价原则",s["h1"]),Image(str(figs["protocol"]),width=235*mm,height=117*mm),Spacer(1,2*mm),P("为什么不能给一个总 AUC：LIT-PCBA 评价跨靶点结合检索；M4 文献集评价功能标签或效力；四上下文 MD 目前只有 3 个化合物。样本宇宙、终点和统计单位不同。把它们标准化后相加会掩盖样本量、标签含义与泄漏风险。因而本报告给出一组整体证据向量，并在样本真正一致的子任务内做公平比较。",s["body"]),PageBreak()]

    # 4 CGDA
    story += [P("3. DrugCLIP 主线：CGDA，而不是只看 family-aug",s["h1"]),P("CGDA（Context-Gated DrugCLIP Adapter）冻结 DrugCLIP 编码器，用靶点口袋与共晶参考配体构造上下文，门控选择两个 rank-2 低秩残差专家。训练使用靶点内 active/decoy、跨靶点检索、困难负样本排序和表示保持损失；留出靶点的标签不进入训练。",s["body"]),Image(str(figs["binding"]),width=240*mm,height=69*mm),Spacer(1,2*mm),table([
        ["Full LIT-PCBA 15T", "ROC-AUC", "PR-AUC", "BEDROC80.5", "EF1%"],
        ["Pocket DrugCLIP", "0.5672", "0.02379", "0.06109", "5.36"],
        ["Reference retrieval", "0.5693", "0.02594", "0.07103", "6.06"],
        ["CGDA hard-rank ensemble", "0.5788", "0.02627", "0.07378", "6.31"],
    ],widths=[60*mm,30*mm,30*mm,35*mm,30*mm],font=8.7),Spacer(1,2*mm),P("相对 reference retrieval，target-bootstrap 的 ROC 增益 +0.00954（95% CI +0.00134 至 +0.01840），BEDROC 增益 +0.00274（+0.00045 至 +0.00504）。这是严格可主张的结果；PR 与多数 EF 截断点不稳定。family-aug 的 13T ROC 0.641 是另一个协议下的域适配结果，不能与 Full LIT-PCBA 直接比较，也不能替代 CGDA。",s["warn"]),PageBreak()]

    # 5 external and M4
    acad=d["m4_external"]["datasets"]["acadia_functional"]; mon=d["m4_external"]["datasets"]["monash_allostery"]
    story += [P("4. 外部迁移与 M4 专化：有效，但并未解决效力",s["h1"]),Image(str(figs["dude"]),width=240*mm,height=69*mm),Spacer(1,2*mm),P("在冻结后的 DUD-E 5-GPCR 外部压力测试中，CGDA 相对 reference retrieval 的 ROC、PR 和 BEDROC 均有正向 target-bootstrap 区间；但原始 pocket DrugCLIP 的 ROC/BEDROC 更高。因此可声称“CGDA 改善 reference-assisted retrieval 的跨 GPCR 迁移”，不能声称“CGDA 普遍优于 DrugCLIP”。",s["body"]),table([
        ["M4 外部来源", "n（正/负）", "Official DrugCLIP", "GPCR-triplet", "结论"],
        ["Acadia 功能阈值", "45（43/2）", "0.837", "0.930", "数值改善；差值 CI 跨零，负样本太少"],
        ["Monash PAM/Inactive", "14（11/3）", "0.576", "0.788", "方向改善；差值 CI 跨零"],
        ["Suven potency", "45", "ρ=-0.066", "ρ=-0.106", "不能排序效力"],
        ["VU6025733 potency", "18", "ρ=-0.643", "ρ=-0.637", "强系列/方向混杂，不是可用效力模型"],
    ],widths=[45*mm,28*mm,35*mm,33*mm,96*mm],font=8.4),Spacer(1,2*mm),P("因此 triplet 的真正价值是改善部分 M4 功能分类与毒蕈碱亚型排序表征；它没有把静态 DrugCLIP 变成 PAM 效力预测器。",s["warn"]),PageBreak()]

    # 6 generalization and static
    story += [P("5. 为什么不能只报一个漂亮 AUC",s["h1"]),Image(str(figs["m4gen"]),width=155*mm,height=82*mm),Spacer(1,2*mm),P("同一 529 分子，DrugCLIP 冻结探针在 scaffold holdout 看似达到 0.856，但 source holdout 降至 0.400；2D ExtraTrees 从 0.919 降至 0.629。说明模型主要受化学系列和来源差异影响。零样本 pocket score 约 0.491，更直接表明“结合相似度”不能替代 PAM 功能。",s["body"]),table([
        ["M4 potency LOSO（12 series）", "macro Spearman", "worst series", "解释"],
        ["Vina only", "0.043", "-0.226", "静态亲和打分几乎无效"],
        ["2D LightGBM", "0.199", "-0.130", "系列内结构信号仍占主导"],
        ["Structure RF", "0.170", "-0.030", "结构特征较稳但增益有限"],
        ["2D + Structure LightGBM", "0.208", "-0.152", "宏均值小幅提升，最差系列恶化"],
    ],widths=[58*mm,40*mm,40*mm,100*mm],font=8.7),Spacer(1,2*mm),P("结论不是“MD 一定更好”，而是静态模型无法承担功能效力终点；动态模块必须用已知 PAM 与 hard negative 做独立机制复核。",s["callout"]),PageBreak()]

    # 7 four context method
    story += [P("6. 四上下文 PACER-DC / PACER-FKG",s["h1"]),P("四个匹配体系：A=apo，P=ACh probe only，C=candidate without probe，CP=candidate+ACh。由同一冻结 C1-BS256 表征生成三种差分：",s["body"]),Image(str(figs["four_context"]),width=210*mm,height=75*mm),Spacer(1,2*mm),table([
        ["差分", "公式", "生物学问题"],
        ["ΔPAM", "CP - P", "在已有 ACh 条件下，加入候选物改变了什么"],
        ["ΔAGO", "C - A", "候选物在无 ACh 时是否自行推动受体状态"],
        ["ΔINT", "CP - P - C + A", "候选物与 ACh 的非加性交互变化"],
    ],widths=[35*mm,55*mm,145*mm],font=9.2),Spacer(1,4*mm),P("编码器并非 OneProt 最终头。当前冻结主线从预训练 Geom2Vec/ViSNet 的 C1 中间 message-passing state 读取原子 invariant 特征（scalar 64 + vector norm 64），按残基 backbone/sidechain 均值拼成 BS256。PACER-FKG 在预定义区域、共同核空间和匹配 replica 中比较差分方向。",s["body"]),P("最重要的统计纠偏：核距离的大小没有方向，且不同区域若使用不同核尺度就不能相减。当前主结果采用共同坐标系中的带符号方向余弦；时间块只作相关子样本，不能冒充独立生物学重复。",s["warn"]),PageBreak()]

    # 8 dynamic conclusion
    story += [P("7. 四上下文结果：partial specificity support",s["h1"]),Image(str(figs["dynamic"]),width=155*mm,height=100*mm),Spacer(1,2*mm),table([
        ["对象", "药理身份", "主区域 ΔINT R1/R3 cosine", "可解释结论"],
        ["compound110", "ago-PAM / allosteric agonist", "+0.584", "局部方向可重复；其他区域并不一致"],
        ["LY2119620", "已知 PAM", "+0.634", "多个机制区为正，支持已知 PAM 的方向模式"],
        ["CM00734", "实验 inactive hard negative", "-0.327", "7 个报告区域均为负，提供冻结 hard-negative 支持"],
    ],widths=[38*mm,55*mm,52*mm,92*mm],font=8.6),Spacer(1,2*mm),P("这是目前最有价值的动态结论：区分力来自跨 replica 的方向一致性，而不是 ΔINT 幅度。CM00734 的幅度仍可达到 LY 的 0.53-0.80，单看‘变化有多大’会误判。由于只有 2 个功能分子和 1 个 hard negative，不能计算可信 AUC、阈值或效力回归，也不能宣称 SOTA classifier。",s["warn"]),PageBreak()]

    # 9 integrated vector
    story += [P("8. 整体算法的一系列指标",s["h1"]),Image(str(figs["vector"]),width=235*mm,height=85*mm),Spacer(1,2*mm),P("这张图不是让柱高彼此竞争，而是展示全流程每个科学终点已有多少证据。当前链条在“通用结合检索”上有大规模严格验证，在“M4 功能迁移”上有小样本外部提示，在“功能效力排序”上仍失败，在“四上下文动态”上只有三化合物机制闭环。",s["body"]),table([
        ["阶段", "当前最强证据", "当前缺口", "比赛中应如何表述"],
        ["结合检索", "CGDA：2.8M LIT-PCBA + 外部 5-GPCR", "未全面胜过 pocket DrugCLIP/所有 SOTA", "严格未见靶点增益"],
        ["M4 专化", "triplet：Acadia/Monash 数值改善", "负样本少，效力不改善", "外部压力测试，不称确认"],
        ["静态结构", "姿势、IFP、LOSO 结构特征", "Vina 与效力弱相关", "作为门控与解释，不作功能分数"],
        ["动态功能", "LY/CM 冻结 hard-negative 对照", "n=3，无通用阈值", "部分特异性机制支持"],
    ],widths=[35*mm,73*mm,68*mm,62*mm],font=8.3),PageBreak()]

    # 10 candidates
    cand=d["candidates"]
    story += [P("9. 候选物交接：有可审计 shortlist，没有已确认 PAM",s["h1"]),Image(str(figs["candidates"]),width=185*mm,height=91*mm),Spacer(1,2*mm),table([
        ["候选", "旧 rank", "family-aug rank", "DrugCLIP p", "famaug score", "max Tanimoto"],
        *[[r.candidate_id,f"{int(r.drugclip_rank)}",f"{int(r.new_rank)}",f"{r.drugclip_m4_probability:.3f}",f"{r.famaug_m4_score:.3f}",f"{r.max_tanimoto:.3f}"] for r in cand.itertuples()],
    ],widths=[38*mm,30*mm,40*mm,38*mm,38*mm,38*mm],font=8.5),Spacer(1,2*mm),P("固定规则为旧 PACER-200 rank≤50 且 family-aug rank≤25。5 个候选均通过 PAINS/drug-like 检查，但属于相近 fragment chemotype；尚未完成各自的四上下文 MD，更没有湿实验。因此这里只能称为优先验证候选。",s["warn"]),PageBreak()]

    # 11 baseline/SOTA
    story += [P("10. 应与哪些算法对标，以及怎样证明 SOTA",s["h1"]),table([
        ["任务", "必须 baseline", "同协议指标", "当前状态"],
        ["大规模结合检索", "DrugCLIP pocket/reference；ECFP；random", "ROC, PR, BEDROC, EF0.5/1/2/5", "CGDA 已完成同协议 LOTO"],
        ["结构重打分", "Vina；GNINA；RTMScore；EquiScore；DeepRLI", "同一 pose/同一靶点/同一候选库", "Vina 已跑；其余只能列 published reference，不能宣称胜出"],
        ["M4 功能分类", "ECFP/RF；DrugCLIP frozen；GPCR-triplet", "source/series holdout + 外部来源", "已完成，但样本不均衡且外部负样本少"],
        ["M4 效力排序", "2D RF/LGBM；static structure fusion", "LOSO Spearman, MAE, worst-series", "整体仍弱；不主张解决"],
        ["MD 表征", "PCA/tICA/VAMP；OneProt；Geom2Vec final", "同轨迹、同窗口、replica holdout", "C1-BS256 方向稳定性更好；仅特定体系"],
        ["PAM 功能动态", "state-only；无图；单差分；置乱图", "预注册区域、方向、hard negative", "三化合物 pilot，尚无领域标准 benchmark"],
    ],widths=[43*mm,73*mm,70*mm,56*mm],font=8.2),Spacer(1,4*mm),P("可接受的 SOTA 证明必须同时满足：同一公开数据、同一输入信息、同一切分、同一指标、预先冻结模型选择、对目标或样本单位 bootstrap，并优于最强可运行 baseline。当前项目只有局部的“best among tested baselines”；不具备“端到端 PAM SOTA”证据。更有力量的竞赛表述是：我们识别了静态结合与功能 PAM 之间的评估缺口，并给出一个带 hard-negative 和拒绝机制的动态复核框架。",s["callout"]),PageBreak()]

    # 12 conclusion refs
    story += [P("11. 最终结论与提交前优先级",s["h1"]),table([
        ["优先级", "必须完成", "理由"],
        ["P0", "冻结 CGDA、family-aug、PACER-FKG 的用途边界，统一候选 ID 与证据卡", "避免把不同协议和不同分子库拼接"],
        ["P0", "对 5 个候选至少选择 1-2 个完成四上下文 MD；同时保留已知 PAM 与 inactive 对照", "把候选链真正接到动态模块"],
        ["P1", "补 GNINA/RTMScore 中至少一个同 pose 重打分 baseline，或明确 NOT RUN", "回应‘docking 太古早’的质疑"],
        ["P1", "扩充跨 chemotype PAM 与 hard negative，预注册 ΔINT 主区域", "把 n=3 pilot 推向可估计性能"],
        ["P2", "湿实验验证 ACh 条件下曲线位移、Emax 与 intrinsic agonism", "唯一能确认 PAM 身份的终点"],
    ],widths=[24*mm,140*mm,78*mm],font=8.6),Spacer(1,5*mm),P("一句话定位",s["h2"]),P("PACER-M4 是一个把通用结合检索、静态结构门控和 ACh 条件化四上下文动态复核连接起来的证据级联系统；其创新点是显式区分“可能结合”与“可能产生功能性协同”，并以冻结 hard negative 和拒绝机制约束结论。",s["callout"]),P("主要参考",s["h2"]),P("DrugCLIP, NeurIPS 2023, https://papers.nips.cc/paper_files/paper/2023/hash/8bd31288ad8e9a31d519fdeede7ee47d-Abstract-Conference.html<br/>LIT-PCBA, JCIM 2020, https://doi.org/10.1021/acs.jcim.0c00155<br/>GNINA 1.3, J. Cheminformatics 2025, https://doi.org/10.1186/s13321-025-00973-x<br/>RTMScore, J. Med. Chem. 2022, https://doi.org/10.1021/acs.jmedchem.2c00991<br/>EquiScore, Nature Machine Intelligence 2024, https://doi.org/10.1038/s42256-024-00849-z<br/>VAMPnets, Nature Communications 2018, https://doi.org/10.1038/s41467-017-02388-1<br/>Geom2Vec, arXiv:2409.19838, https://arxiv.org/abs/2409.19838",s["small"])]

    doc.build(story)
    return pdf


def main():
    configure_fonts(); ASSET.mkdir(parents=True,exist_ok=True)
    d,paths=collect(); save_scorecard(d,paths)
    figs={
        "architecture":fig_architecture(),
        "protocol":fig_protocol_map(),
        "binding":fig_binding(d),
        "dude":fig_dude(d),
        "m4gen":fig_m4_generalization(d),
        "dynamic":fig_dynamic(d),
        "four_context":fig_four_context(),
        "vector":fig_evidence_vector(d),
        "candidates":fig_candidates(d),
    }
    pdf=build_pdf(d,figs)
    print(pdf)


if __name__ == "__main__":
    main()
