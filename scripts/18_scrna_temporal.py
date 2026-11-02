"""18_scrna_temporal.py -- temporal dynamics (Sham vs SAH-24h vs SAH-72h) of the heme program and
neuronal ferroptosis markers, plus cell-type composition shifts (neutrophil infiltration)."""
import os, numpy as np, pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
from matplotlib import rcParams
sc.settings.verbosity=0
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,FIG=os.path.join(BASE,"results"),os.path.join(BASE,"figures")
SCRATCH="H:/徐伟明论文/sah_rbc_bioinfo/_processed"
OKABE={"blue":"#0072B2","orange":"#E69F00","green":"#009E73","vermil":"#D55E00","purple":"#CC79A7","grey":"#9AA0A6"}
rcParams.update({"font.family":"Arial","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})

a=sc.read_h5ad(os.path.join(SCRATCH,"sah_scrna_annot.h5ad"))
# neuronal / general ferroptosis execution score (pro-ferroptosis genes)
ferro=[g for g in ["Acsl4","Ptgs2","Slc7a11","Hmox1","Tfrc","Slc40a1","Alox15"] if g in a.var_names]
sc.tl.score_genes(a,ferro,score_name="ferro_score")
a.obs["cond"]=a.obs.condition
conds=["Sham","SAH_24h","SAH_72h"]

# --- heme + ferroptosis score by cell type x condition ---
rows=[]
for ct in ["Neutrophil","Microglia/Mac","Neuron"]:
    d=a[a.obs.cell_type==ct]
    for c in conds:
        for score in ["heme_score","ferro_score"]:
            v=d.obs.loc[d.obs.cond==c,score]
            rows.append({"cell_type":ct,"condition":c,"score":score,"mean":v.mean(),"n":len(v)})
TS=pd.DataFrame(rows); TS.to_csv(os.path.join(RES,"scrna_temporal_scores.csv"),index=False)
# Kruskal test per cell type/score across conditions
tests=[]
for ct in ["Neutrophil","Microglia/Mac","Neuron"]:
    for score in ["heme_score","ferro_score"]:
        d=a[a.obs.cell_type==ct]
        groups=[d.obs.loc[d.obs.cond==c,score].values for c in conds]
        h,p=stats.kruskal(*groups)
        tests.append({"cell_type":ct,"score":score,"kruskal_H":h,"p":p})
TT=pd.DataFrame(tests); TT.to_csv(os.path.join(RES,"scrna_temporal_tests.csv"),index=False)
print("temporal tests (Sham vs 24h vs 72h):\n",TT.round(4).to_string(index=False))
print("\nheme/ferro means:\n",TS.pivot_table(index=["cell_type","score"],columns="condition",values="mean").round(3).to_string())

# --- composition by condition (fraction per cell type) ---
comp=pd.crosstab(a.obs.cell_type,a.obs.cond,normalize="columns")*100
comp.to_csv(os.path.join(RES,"scrna_temporal_composition.csv"))
print("\ncell-type composition (% per condition):\n",comp.round(1).to_string())

# Fig13
fig,axes=plt.subplots(1,3,figsize=(14,4.2),gridspec_kw={"width_ratios":[1.3,1.1,1]})
ax=axes[0]
piv=TS[TS.score=="heme_score"].pivot(index="cell_type",columns="condition",values="mean")
x=np.arange(len(piv)); w=0.26
for i,c in enumerate(conds):
    ax.bar(x+(i-1)*w,piv[c],w,color={"Sham":OKABE["grey"],"SAH_24h":OKABE["orange"],"SAH_72h":OKABE["vermil"]}[c],label=c)
ax.set_xticks(x); ax.set_xticklabels(piv.index,fontsize=8); ax.set_ylabel("heme/iron score")
ax.set_title("Heme program: temporal",fontsize=9); ax.legend(fontsize=6.5)
ax=axes[1]
piv=TS[TS.score=="ferro_score"].pivot(index="cell_type",columns="condition",values="mean")
for i,c in enumerate(conds):
    ax.bar(x+(i-1)*w,piv[c],w,color={"Sham":OKABE["grey"],"SAH_24h":OKABE["orange"],"SAH_72h":OKABE["vermil"]}[c],label=c)
ax.set_xticks(x); ax.set_xticklabels(piv.index,fontsize=8); ax.set_ylabel("ferroptosis score")
ax.set_title("Ferroptosis exec: temporal",fontsize=9)
ax=axes[2]
c2=comp.loc[["Neutrophil","Microglia/Mac","T cell","Astrocyte","Neuron"]]
c2.plot(kind="bar",ax=ax,color=[OKABE["vermil"],OKABE["green"],OKABE["blue"],OKABE["orange"],OKABE["purple"]],width=0.7)
ax.set_ylabel("% of cells"); ax.set_title("Composition shift",fontsize=9); ax.tick_params(labelsize=7)
handles,labels=ax.get_legend_handles_labels(); ax.legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig13_scrna_temporal.png"),dpi=300,facecolor="white")
fig.savefig(os.path.join(FIG,"Fig13_scrna_temporal.pdf")); plt.close(fig)
print("saved Fig13_scrna_temporal")
print("DONE 18")
