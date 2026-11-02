"""09_figures.py -- publication figures (Okabe-Ito palette, Arial, individual PNG+vector PDF panels)."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT,FIG=os.path.join(BASE,"results"),os.path.join(BASE,"data"),os.path.join(BASE,"figures")
os.makedirs(FIG,exist_ok=True)
OKABE={"blue":"#0072B2","orange":"#E69F00","sky":"#56B4E9","green":"#009E73",
       "yellow":"#F0E442","vermil":"#D55E00","purple":"#CC79A7","black":"#000000"}
rcParams.update({"font.family":"Arial","font.size":9,"axes.linewidth":0.9,
                 "axes.labelsize":9,"xtick.labelsize":8,"ytick.labelsize":8,
                 "axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})
def save(fig,name,w,h):
    fig.set_size_inches(w/2.54,h/2.54)
    fig.savefig(os.path.join(FIG,name+".png"),dpi=300,bbox_inches="tight",facecolor="white")
    fig.savefig(os.path.join(FIG,name+".pdf"),bbox_inches="tight")
    plt.close(fig); print("saved",name)

deg=pd.read_csv(os.path.join(RES,"DEG_GSE36791.csv"))
cand=pd.read_csv(os.path.join(RES,"candidate_RBC_SAH_genes.csv"))
hub=pd.read_csv(os.path.join(RES,"hub_genes_annotated.csv"))

# ---------- Fig1 volcano ----------
fig,ax=plt.subplots()
x=deg.logFC; y=-np.log10(deg["P.Value"].clip(lower=1e-300))
ns=(deg["adj.P.Val"]>=0.05); up=(deg["adj.P.Val"]<0.05)&(deg.logFC>=1); dn=(deg["adj.P.Val"]<0.05)&(deg.logFC<=-1)
ax.scatter(x[ns],y[ns],s=6,c="#B9BEC4",linewidths=0)
ax.scatter(x[up],y[up],s=18,c=OKABE["vermil"],linewidths=0)
ax.scatter(x[dn],y[dn],s=18,c=OKABE["blue"],linewidths=0)
lab=["HP","MARC1","NFE2","CST7","IL2RB","SLC40A1","FTH1","SLC25A37","GYPB","MMP9","S100A12","ANXA3","EOMES","FCER1A"]
_texts=[]
for g in lab:
    r=deg[deg.gene==g]
    if len(r):
        r=r.iloc[0]
        _texts.append(ax.annotate(g,(r.logFC,-np.log10(max(r["P.Value"],1e-300))),fontsize=6.5,color=OKABE["black"]))
try:
    from adjustText import adjust_text
    adjust_text(_texts, ax=ax, expand=(1.9,2.6), expand_text=(1.4,2.2), maxiter=1200,
                min_arrow_len=4, force_text=(1.2,1.8), force_points=(0.6,0.6),
                arrowprops=dict(arrowstyle="-",lw=0.4,color="#888888"))
except Exception as e:
    print("adjustText skipped:", e)
ax.axhline(-np.log10(0.05),ls="--",lw=0.7,c="grey"); ax.axvline(0,ls=":",lw=0.7,c="grey")
ax.axvline(1,color="grey",ls=":",lw=0.7); ax.axvline(-1,color="grey",ls=":",lw=0.7)
ax.set_xlabel("log2FC (SAH vs Control)"); ax.set_ylabel(r"$-$log$_{10}$ P")
ax.set_title(f"GSE36791 acute SAH blood  (up={up.sum()}, down={dn.sum()})")
save(fig,"Fig1_volcano_GSE36791",12,9)

# ---------- Fig2 overlap DEGs vs RBC signature ----------
fig,ax=plt.subplots()
sig=set(open(os.path.join(RES,"RBC_signature_genes.txt")).read().split())
nA=len(deg[deg["adj.P.Val"]<0.05]); nB=len(sig); nAB=len(cand)
rect=lambda xs,ys,w_,h_,c: ax.add_patch(plt.Rectangle((xs,ys),w_,h_,fc="none",ec=c,lw=2))
ax.bar([0],[nA],color=OKABE["blue"],width=0.5); ax.bar([1],[nB],color=OKABE["green"],width=0.5)
ax.bar([2],[nAB],color=OKABE["vermil"],width=0.5)
for i,v in enumerate([nA,nB,nAB]): ax.text(i,v+15,str(v),ha="center",fontsize=9,weight="bold")
ax.set_xticks([0,1,2]); ax.set_xticklabels(["SAH DEGs\n(FDR<0.05)","RBC/HRR\nsignature","Candidates\n(intersection)"],fontsize=8)
ax.set_ylabel("# genes"); ax.set_title("Cross-analysis: SAH transcriptome $\\cap$ red-cell axis")
save(fig,"Fig2_overlap",8,6.5)

# ---------- Fig3 heatmap candidates ----------
def lnorm(tag):
    lg=np.load(os.path.join(DAT,f"lnorm_{tag}.npy")); g=np.array(open(os.path.join(DAT,f"genes_{tag}.txt")).read().split())
    ph=pd.read_csv(os.path.join(DAT,f"expr_linear_{tag}.tsv"),sep="\t",index_col=0).columns
    pd_=pd.read_csv(os.path.join(DAT,f"pheno_aligned_{tag}.csv")).set_index("gsm").loc[list(ph)]
    return lg,g,pd_["group"].values
genes=list(dict.fromkeys(list(hub.gene)+list(cand.sort_values("P.Value").gene)))
lg,g,grp=lnorm("GSE36791"); gi={x:i for i,x in enumerate(g)}
rows=[gi[x] for x in genes if x in gi]; M=lg[rows]
Z=(M-M.mean(1,keepdims=True))/(M.std(1,keepdims=True)+1e-9)
order=np.argsort((grp=="SAH"))
fig,ax=plt.subplots()
im=ax.imshow(Z[:,order],aspect="auto",cmap="RdBu_r",vmin=-2.5,vmax=2.5)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([g[i] for i in rows],fontsize=6)
ax.set_xticks([]); ax.set_title("Candidate + hub genes, SAH (right) vs Control (left)\nz-scored, GSE36791")
cb=fig.colorbar(im,ax=ax,fraction=0.03,pad=0.02); cb.set_label("row z",fontsize=7)
save(fig,"Fig3_heatmap_candidates",10,8)

# ---------- Fig4 enrichment (ORA + GSEA focus) ----------
fig,axes=plt.subplots(1,2,gridspec_kw={"wspace":1.15})
ax=axes[0]
ora=pd.read_csv(os.path.join(RES,"ORA_candidates_GO_Biological_Process_2023.csv")).head(8).iloc[::-1]
ax.barh(range(len(ora)),-np.log10(ora["Adjusted P-value"]),color=OKABE["green"])
ax.set_yticks(range(len(ora))); ax.set_yticklabels([t.split(" (")[0][:38] for t in ora.Term],fontsize=6)
ax.set_xlabel(r"$-$log$_{10}$ adj.P (candidates)"); ax.set_title("GO BP — candidate genes",fontsize=8)
ax=axes[1]
g=pd.read_csv(os.path.join(RES,"GSEA_GSE36791_all.csv"))
focus=["H:heme Metabolism","H:Reactive Oxygen Species Pathway","H:TNF-alpha Signaling via NF-kB",
       "H:Hypoxia","K:HIF-1 signaling pathway","K:Ferroptosis","K:Complement and coagulation cascades",
       "R:Neutrophil Degranulation","H:IL-6/JAK/STAT3 Signaling"]
sub=g[g.Term.isin(focus)].sort_values("NES")
cols=[OKABE["vermil"] if v>0 else OKABE["blue"] for v in sub.NES]
ax.barh(range(len(sub)),sub.NES,color=cols)
ax.set_yticks(range(len(sub))); ax.set_yticklabels([t.split(":")[1][:34] for t in sub.Term],fontsize=6)
ax.axvline(0,lw=0.8,c="k"); ax.set_xlabel("Normalised enrichment (SAH)"); ax.set_title("GSEA — pathways of interest",fontsize=8)
save(fig,"Fig4_enrichment",14,7)

# ---------- Fig5 ML hub ranking ----------
fig,ax=plt.subplots()
rank=pd.read_csv(os.path.join(RES,"ml_ranking_all.csv"),index_col=0)
top=rank.head(12).iloc[::-1]
b=0.35; y=np.arange(len(top))
ax.barh(y+b/2,top["lasso"],b*0.9,color=OKABE["blue"],label="LASSO")
ax.barh(y-b/2,top["rf_gini"],b*0.9,color=OKABE["orange"],label="Random forest")
ax.barh(y,top["svm_rfe"]*1.001,b*0.5,color=OKABE["green"],label="SVM-RFE")
ax.set_yticks(y); ax.set_yticklabels(top.index,fontsize=7)
ax.set_xlabel("percentile importance (0=most important)"); ax.invert_xaxis()
ax.legend(fontsize=7,loc="lower right"); ax.set_title("ML feature ranking (top 12)")
save(fig,"Fig5_ml_ranking",9,7)

# ---------- Fig6 program score + ROC + external ----------
fig=plt.figure(); gs=fig.add_gridspec(1,3,wspace=0.5)
ax=fig.add_subplot(gs[0,0])
sc=pd.read_csv(os.path.join(RES,"program_activity_scores.csv"))
for i,(coh,cc) in enumerate([("GSE36791",OKABE["blue"]),("GSE73378",OKABE["orange"])]):
    d=sc[sc.cohort==coh]
    a=d[d.group=="SAH"].score; b_=d[d.group=="Control"].score
    pos=np.arange(2)+i*3.0
    bp=ax.boxplot([a,b_],positions=pos,widths=0.7,patch_artist=True,showfliers=False)
    for p,c in zip(bp["boxes"],[cc,"#9AA0A6"]): p.set_facecolor(c); p.set_alpha(0.7); p.set_linewidth(0.6)
    for m in bp["medians"]: m.set_color("k")
    ax.scatter(pos[0]+np.random.RandomState(1).normal(0,.08,len(a)),a,s=6,c="k",alpha=.4,zorder=5)
    ax.scatter(pos[1]+np.random.RandomState(2).normal(0,.08,len(b_)),b_,s=6,c="k",alpha=.4,zorder=5)
ax.axvline(2.0,color="grey",lw=0.7,ls="--")
ax.text(0.5,0.98,"acute SAH",ha="center",fontsize=6.5,transform=ax.get_xaxis_transform())
ax.text(3.5,0.98,"remote aSAH",ha="center",fontsize=6.5,transform=ax.get_xaxis_transform())
ax.set_xticks([0,1,3,4]); ax.set_xticklabels(["SAH","Ctrl","SAH","Ctrl"],fontsize=7)
ax.set_ylabel("RBC/HRR program score"); ax.set_title("bulk program: acute p=0.13 (ns)\nremote p=0.88 (ns) — see gene/pathway evidence",fontsize=7)
ax.set_xlim(-0.8,4.8)
ax=fig.add_subplot(gs[0,1])
r=pd.read_csv(os.path.join(RES,"roc_hub_discovery.csv"))
ax.plot(r.fpr,r.tpr,color=OKABE["vermil"],lw=1.5,label="discovery CV AUC=0.90")
ax.plot([0,1],[0,1],ls="--",lw=0.8,c="grey")
ax.set_xlabel("1-Specificity"); ax.set_ylabel("Sensitivity"); ax.set_title("Hub score ROC (discovery)",fontsize=8); ax.legend(fontsize=6)
ax=fig.add_subplot(gs[0,2])
s=pd.read_csv(os.path.join(RES,"score_validation_summary.csv"))
xx=np.arange(len(s)); w=0.38
ax.bar(xx-w/2,s.AUC_discovery_CV20x5,w,color=OKABE["blue"],label="discovery CV")
ax.bar(xx+w/2,s.AUC_external_GSE73378,w,color=OKABE["orange"],label="external (remote aSAH)")
ax.axhline(0.5,ls="--",lw=0.8,c="grey"); ax.axhline(0.75,ls=":",lw=0.8,c=OKABE["green"])
ax.set_xticks(xx); ax.set_xticklabels([m.replace("M_","") for m in s.model],fontsize=6)
ax.set_ylabel("AUC"); ax.set_ylim(0,1); ax.set_title("Internal vs external AUC",fontsize=8); ax.legend(fontsize=6)
save(fig,"Fig6_score_and_validation",16,6)

print("DONE 09_figures")
