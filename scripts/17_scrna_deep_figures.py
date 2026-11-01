"""17_scrna_deep_figures.py -- Fig12: cell-type-specific axis DE, neuron ferroptosis panel, regulon activity."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import rcParams
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,FIG=os.path.join(BASE,"results"),os.path.join(BASE,"figures")
OKABE={"blue":"#0072B2","orange":"#E69F00","green":"#009E73","vermil":"#D55E00","purple":"#CC79A7"}
rcParams.update({"font.family":"Arial","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})
def save(fig,name,w,h):
    fig.set_size_inches(w/2.54,h/2.54)
    fig.savefig(os.path.join(FIG,name+".png"),dpi=300,bbox_inches="tight",facecolor="white")
    fig.savefig(os.path.join(FIG,name+".pdf"),bbox_inches="tight"); plt.close(fig); print("saved",name)

DE=pd.read_csv(os.path.join(RES,"scrna_celltype_DE_panels.csv"))
REG=pd.read_csv(os.path.join(RES,"scrna_regulon_proxy.csv"))

fig=plt.figure(figsize=(13,4.2)); gs=fig.add_gridspec(1,3,wspace=0.55)
# A: log2FC heatmap by cell type
ax=fig.add_subplot(gs[0,0])
piv=DE.pivot_table(index="gene",columns="cell_type",values="log2FC")
order=["Neutrophil","Microglia/Mac","Neuron"]
piv=piv.reindex(columns=[c for c in order if c in piv.columns])
im=ax.imshow(piv.values,cmap="RdBu_r",vmin=-3,vmax=3,aspect="auto")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels([c.replace("/","/\n") for c in piv.columns],fontsize=7)
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index,fontsize=7)
ax.set_title("SAH−Sham log2FC\n(cell-type-specific)",fontsize=8)
fig.colorbar(im,ax=ax,fraction=0.04)
# B: neuron ferroptosis panel
ax=fig.add_subplot(gs[0,1])
neu=DE[DE.cell_type=="Neuron"].sort_values("log2FC")
cols=[OKABE["vermil"] if v>0 else OKABE["blue"] for v in neu.log2FC]
ax.barh(range(len(neu)),neu.log2FC,color=cols)
ax.set_yticks(range(len(neu))); ax.set_yticklabels(neu.gene,fontsize=8)
ax.axvline(0,lw=.8,c="k"); ax.set_xlabel("log2FC (SAH−Sham)")
ax.set_title("Neuronal ferroptosis program\n(Gpx4 down; Acsl4/Slc7a11/Ptgs2/Hmox1 up)",fontsize=8)
# C: regulon activity SAH vs Sham
ax=fig.add_subplot(gs[0,2])
x=np.arange(len(REG)); w=0.38
ax.bar(x-w/2,REG.Sham_activity,w,color="#9AA0A6",label="Sham")
ax.bar(x+w/2,REG.SAH_activity,w,color=OKABE["green"],label="SAH")
for i,r in REG.iterrows():
    ax.text(i+w/2,r.SAH_activity+0.02,f"p={r.SAH_vs_Sham_p:.1e}",ha="center",fontsize=5.5)
ax.set_xticks(x); ax.set_xticklabels([s.replace("_"," ") for s in REG.regulon],fontsize=6.5,rotation=20)
ax.set_ylabel("regulon activity score"); ax.legend(fontsize=7)
ax.set_title("SCENIC-lite curated regulons\n(NRF2 up in SAH)",fontsize=8)
save(fig,"Fig12_scrna_deep",13,4.6)
print("DONE 17")
