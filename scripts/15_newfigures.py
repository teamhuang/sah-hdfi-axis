"""15_newfigures.py -- WGCNA, single-cell summary, XGBoost comparison, axis schematic."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,FIG=os.path.join(BASE,"results"),os.path.join(BASE,"figures")
OKABE={"blue":"#0072B2","orange":"#E69F00","sky":"#56B4E9","green":"#009E73","vermil":"#D55E00","purple":"#CC79A7"}
rcParams.update({"font.family":"Arial","font.size":9,"axes.linewidth":0.9,
                 "axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})
def save(fig,name,w,h):
    fig.set_size_inches(w/2.54,h/2.54)
    fig.savefig(os.path.join(FIG,name+".png"),dpi=300,bbox_inches="tight",facecolor="white")
    fig.savefig(os.path.join(FIG,name+".pdf"),bbox_inches="tight"); plt.close(fig); print("saved",name)

# ---------- Fig7 WGCNA ----------
mt=pd.read_csv(os.path.join(RES,"wgcna_module_trait.csv"),index_col=0)
net=pd.read_csv(os.path.join(RES,"wgcna_gene_network.csv"))
fig,axes=plt.subplots(1,2,gridspec_kw={"width_ratios":[1,1.1]})
ax=axes[0]
mt=mt.sort_values("cor_SAH")
cols=[OKABE["vermil"] if p<0.05 else "#B9BEC4" for p in mt.p_SAH]
ax.barh(range(len(mt)),mt.cor_SAH,color=cols)
ax.set_yticks(range(len(mt))); ax.set_yticklabels(mt.index,fontsize=6)
ax.set_xlabel("module–SAH correlation"); ax.axvline(0,lw=.8,c="k")
ax.set_title("WGCNA modules vs SAH\n(coloured: FDR-signif. eigengene)",fontsize=8)
ax=axes[1]
top=mt.index[0]; tmod=int(top.replace("ME",""))
sub=net[net.module==tmod]
ax.scatter(sub.GS,sub.MM,s=6,c="#B9BEC4",linewidths=0)
axis=["HP","HMOX1","FTH1","SLC40A1","NFE2","MARC1","IL2RB","BCL11A","GPX4","ACSL4","SLC25A37","CST7"]
for g in axis:
    r=net[net.gene==g]
    if len(r):
        r=r.iloc[0]; ax.scatter(r.GS,r.MM,s=30,c=OKABE["blue"],zorder=5)
        ax.annotate(g,(r.GS,r.MM),fontsize=6.5,xytext=(2,2),textcoords="offset points")
ax.set_xlabel("gene significance (SAH)"); ax.set_ylabel("module membership (MM)")
ax.set_title(f"MM vs GS, {top} (r={mt.loc[top,'cor_SAH']:.2f})",fontsize=8)
save(fig,"Fig7_wgcna",14,7)

# ---------- Fig8 single-cell summary ----------
R=pd.read_csv(os.path.join(RES,"scrna_heme_score_by_celltype.csv"))
me=pd.read_csv(os.path.join(RES,"scrna_marker_meanexpr.csv"),index_col=0)
fig,axes=plt.subplots(1,2,gridspec_kw={"width_ratios":[1,1.2]})
ax=axes[0]
R=R.sort_values("delta")
cols=[OKABE["vermil"] if p<0.001 else OKABE["orange"] if p<0.05 else "#B9BEC4" for p in R.MWU_p]
ax.barh(range(len(R)),R.delta,color=cols)
ax.set_yticks(range(len(R))); ax.set_yticklabels(R.cell_type,fontsize=7)
ax.set_xlabel("heme/iron score Δ (SAH − Sham)"); ax.axvline(0,lw=.8,c="k")
ax.set_title("SAH raises heme/iron program\nacross brain cell types (MWU)",fontsize=8)
ax=axes[1]
mk=[g for g in ["Hp","Cd163","Mrc1","Hmox1","Fth1","Slc40a1","Nfe2"] if g in me.columns]
M=me[mk]
Mz=(M-M.mean())/ (M.std()+1e-9)   # per-gene z across cell types to reveal specificity
im=ax.imshow(Mz.values,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
ax.set_xticks(range(len(mk))); ax.set_xticklabels(mk,rotation=45,ha="right",fontsize=7)
ax.set_yticks(range(len(M.index))); ax.set_yticklabels(M.index,fontsize=6.5)
ax.set_title("axis markers: cell-type z-score\n(mean expr standardised per gene)",fontsize=8)
fig.colorbar(im,ax=ax,fraction=0.03)
save(fig,"Fig8_scrna_summary",14,7)

# ---------- Fig9 XGBoost comparison ----------
X=pd.read_csv(os.path.join(RES,"xgb_model_comparison.csv"))
fig,ax=plt.subplots()
xx=np.arange(len(X)); w=0.38
ax.bar(xx-w/2,X.CV_AUC_mean,w,color=OKABE["blue"],label="discovery CV AUC",yerr=X.CV_AUC_sd,capsize=3)
ax.bar(xx+w/2,X.external_AUC_GSE73378,w,color=OKABE["orange"],label="external AUC (remote aSAH)")
ax.axhline(0.5,ls="--",lw=.8,c="grey"); ax.axhline(0.75,ls=":",lw=.8,c=OKABE["green"])
ax.set_xticks(xx); ax.set_xticklabels([s.replace("_"," ") for s in X.feature_set],fontsize=7)
ax.set_ylabel("AUC"); ax.set_ylim(0,1); ax.legend(fontsize=7)
ax.set_title("Gradient-boosting: heme panel generalises best externally",fontsize=8)
save(fig,"Fig9_xgb",9,6)

# ---------- Fig10 axis schematic ----------
fig,ax=plt.subplots(); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,6.4)
def box(x,y,w,h,txt,c,fs=8.5):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.10",fc=c,ec="none",alpha=0.92))
    ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=fs,color="white",weight="bold")
def arrow(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=13,lw=1.6,color="#444"))
# row 1 (left->right)
y1=4.7; w=1.72; h=1.05; xs=[0.15,2.15,4.15,6.15,8.15]
box(xs[0],y1,w,h,"SAH\nRBC lysis",OKABE["vermil"])
box(xs[1],y1,w,h,"free Hb\nup",OKABE["orange"])
box(xs[2],y1,w,h,"HP up\n(Hb scavenger)",OKABE["green"])
box(xs[3],y1,w,h,"CD163+\nmacrophage",OKABE["blue"])
box(xs[4],y1,w,h,"HMOX1\nheme->Fe",OKABE["purple"])
for i in range(4): arrow(xs[i]+w,y1+h/2,xs[i+1],y1+h/2)
# row 2 (left->right), connected by elbow from HMOX1
y2=2.5
box(xs[1],y2,w,h,"iron / ROS\nFTH1·SLC40A1",OKABE["orange"],fs=7.5)
box(xs[2],y2,w,h,"ferroptosis\nGPX4-low·ACSL4",OKABE["vermil"],fs=7.5)
box(xs[3],y2,w,h,"neuronal\ninjury",OKABE["purple"],fs=8)
arrow(xs[4]+w/2,y1, xs[1]+w/2, y2+h)          # HMOX1 -> iron (elbow across)
arrow(xs[1]+w,y2+h/2,xs[2],y2+h/2)            # iron -> ferroptosis
arrow(xs[2]+w,y2+h/2,xs[3],y2+h/2)            # ferroptosis -> neuronal
box(0.15,0.4,9.7,1.1,"HRR = Hb / RDW  —  peripheral index of the oxygen-carrying (Hb)\nand erythrocyte-turnover (RDW) arms of this axis",OKABE["sky"],fs=8.5)
ax.text(5,6.1,"HP-CD163-HMOX1-ferroptosis axis after SAH (integrative evidence)",ha="center",fontsize=11,weight="bold")
save(fig,"Fig10_axis_schematic",14,6.4)
print("DONE 15_newfigures")
