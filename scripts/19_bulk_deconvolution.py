"""19_bulk_deconvolution.py -- human-blood (GSE36791) marker-based immune deconvolution +
heme/iron and ferroptosis panel scores: SAH vs control, and ferroptosis~heme correlation.
Addresses the composition-confounding limitation quantitatively."""
import os, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
from matplotlib import rcParams
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT=os.path.join(BASE,"results"),os.path.join(BASE,"data"); FIG=os.path.join(BASE,"figures")
OKABE={"blue":"#0072B2","orange":"#E69F00","green":"#009E73","vermil":"#D55E00","purple":"#CC79A7","grey":"#9AA0A6"}
rcParams.update({"font.family":"Arial","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})

lg=np.load(os.path.join(DAT,"lnorm_GSE36791.npy")); g=np.array(open(os.path.join(DAT,"genes_GSE36791.txt")).read().split())
gi={x:i for i,x in enumerate(g)}
ph=pd.read_csv(os.path.join(DAT,"expr_linear_GSE36791.tsv"),sep="\t",index_col=0).columns.tolist()
pheno=pd.read_csv(os.path.join(DAT,"pheno_aligned_GSE36791.csv")).set_index("gsm").loc[ph]
grp=pheno["group"].values; is_sah=grp=="SAH"
Z=(lg-lg.mean(1,keepdims=True))/(lg.std(1,keepdims=True)+1e-9)

sets={
 "Neutrophil":["FCGR3B","CSF3R","ELANE","MPO","CXCR2","S100A8","S100A9","CAMP","LCN2","ITGAM","PGLYRP1"],
 "Monocyte":["CD14","FCN1","VCAN","S100A12","CD68","CSF1R","ANXA3","SERPINA1"],
 "T cell":["CD3D","CD3E","IL7R","TRAC","CD2","CD27","CD6","CCR7","IL2RB"],
 "B cell":["CD79A","CD79B","MS4A1","BANK1","PAX5"],
 "NK":["GNLY","NKG7","KLRD1","PRF1","EOMES","CD247"],
 "Heme_iron":["HP","HMOX1","FTH1","FTL","SLC40A1","SLC25A37","TFRC","AHSP","FECH","NFE2"],
 "Ferroptosis":["GPX4","ACSL4","SLC7A11","HMOX1","FTH1","FTL","TFRC","PTGS2","ALOX15","SLC40A1"],
}
scores={}
for name,gl in sets.items():
    idx=[gi[x] for x in gl if x in gi]
    scores[name]=Z[idx].mean(axis=0)
    print(f"{name}: {len(idx)}/{len(gl)} markers | SAH mean={Z[idx].mean(axis=0)[is_sah].mean():.3f} Ctrl={Z[idx].mean(axis=0)[~is_sah].mean():.3f}")
S=pd.DataFrame(scores); S["group"]=grp
rows=[]
for name in sets:
    a=S.loc[S.group=="SAH",name]; b=S.loc[S.group=="Control",name]
    t,p=stats.ttest_ind(a,b); u,pu=stats.mannwhitneyu(a,b)
    d=(a.mean()-b.mean())/np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
    rows.append({"signature":name,"SAH_mean":a.mean(),"Ctrl_mean":b.mean(),"delta":a.mean()-b.mean(),
                 "cohens_d":d,"t_p":p,"MWU_p":pu})
DE=pd.DataFrame(rows).sort_values("MWU_p"); DE.to_csv(os.path.join(RES,"bulk_deconvolution_scores.csv"),index=False)
print("\nSAH vs Control signature scores:\n",DE.round(4).to_string(index=False))

# ferroptosis ~ heme correlation across samples
rho,p=stats.spearmanr(S.Ferroptosis,S.Heme_iron)
rho2,p2=stats.spearmanr(S.Ferroptosis,S.Neutrophil)
print(f"\nFerroptosis ~ Heme_iron (Spearman rho={rho:.2f}, p={p:.2g}); ~ Neutrophil (rho={rho2:.2f}, p={p2:.2g})")
pd.DataFrame({"pair":["Ferroptosis~Heme","Ferroptosis~Neutrophil"],"rho":[rho,rho2],"p":[p,p2]}).to_csv(
    os.path.join(RES,"bulk_ferroptosis_heme_correlation.csv"),index=False)

# Fig14
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
ax=axes[0]
order=["Neutrophil","Monocyte","Heme_iron","Ferroptosis","T cell","NK","B cell"]
dd=DE.set_index("signature").loc[order].sort_values("delta")
cols=[OKABE["vermil"] if v>0 else OKABE["blue"] for v in dd.delta]
ax.barh(range(len(dd)),dd.delta,color=cols)
for i,(nm,r) in enumerate(dd.iterrows()):
    ax.text(r.delta+(0.02 if r.delta>=0 else -0.02),i,f"p={r.MWU_p:.1e}",va="center",
            ha="left" if r.delta>=0 else "right",fontsize=6)
ax.set_yticks(range(len(dd))); ax.set_yticklabels(dd.index,fontsize=8); ax.axvline(0,lw=.8,c="k")
ax.set_xlabel("Δ score (SAH − Control)"); ax.set_title("Blood immune/heme/ferroptosis signatures",fontsize=9)
ax=axes[1]
ax.scatter(S.loc[S.group=="Control","Heme_iron"],S.loc[S.group=="Control","Ferroptosis"],
           s=28,c=OKABE["blue"],label="Control")
ax.scatter(S.loc[S.group=="SAH","Heme_iron"],S.loc[S.group=="SAH","Ferroptosis"],
           s=28,c=OKABE["vermil"],label="SAH")
ax.set_xlabel("heme/iron score"); ax.set_ylabel("ferroptosis score")
ax.set_title(f"Ferroptosis ~ heme (rho={rho:.2f}, p={p:.1e})",fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig14_bulk_deconvolution.png"),dpi=300,facecolor="white")
fig.savefig(os.path.join(FIG,"Fig14_bulk_deconvolution.pdf")); plt.close(fig)
print("saved Fig14_bulk_deconvolution")
print("DONE 19")
