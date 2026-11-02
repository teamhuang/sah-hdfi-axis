"""16_scrna_deep.py -- deep single-cell: (1) cell-type-specific SAH vs Sham DE for the axis panels,
(2) microglia/macrophage pseudotime trajectory (homeostatic->heme-processing), (3) SCENIC-lite
curated-regulon activity for NFE2L2(NRF2)/NFE2/BACH1 (a labelled proxy; true pySCENIC needs motif DB)."""
import os, numpy as np, pandas as pd
import scanpy as sc, anndata as ad
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
OKABE={"blue":"#0072B2","orange":"#E69F00","green":"#009E73","vermil":"#D55E00","purple":"#CC79A7","sky":"#56B4E9"}
rcParams.update({"font.family":"Arial","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})

adata=sc.read_h5ad(os.path.join(SCRATCH,"sah_scrna_annot.h5ad"))
print("loaded",adata.shape,"| X is log-norm:",float(adata.X.max())<20)

panels={
 "Neutrophil":["Hp","Hmox1","Fth1","Ftl1","S100a8","S100a9","Lcn2","Mmp9"],
 "Microglia/Mac":["Cd163","Hmox1","Mrc1","Slc40a1","Fth1","Acs14","Gpx4"],
 "Neuron":["Slc7a11","Gpx4","Acsl4","Ptgs2","Hmox1"],
}
rows=[]
for ct,gl in panels.items():
    d=adata[adata.obs.cell_type==ct].copy()
    if d.n_obs<50: continue
    d.obs['grp']=pd.Categorical(pd.Series(d.obs.is_sah.values).map({True:'SAH',False:'Sham'}).values)
    sc.tl.rank_genes_groups(d, groupby="grp", groups=['SAH'], reference='Sham', method="wilcoxon")
    rg=d.uns["rank_genes_groups"]
    pmap=dict(zip(np.asarray(rg['names']['SAH']), np.asarray(rg['pvals_adj']['SAH'])))
    sah_mask=(d.obs.grp=='SAH').values; sham_mask=(d.obs.grp=='Sham').values
    for g in gl:
        if g not in d.var_names: 
            rows.append({"cell_type":ct,"gene":g,"SAH_mean":np.nan,"Sham_mean":np.nan,"log2FC":np.nan,"p_adj":np.nan}); continue
        col=np.asarray(d[:,g].X.todense()).ravel()
        ms=col[sah_mask].mean(); mc=col[sham_mask].mean()
        rows.append({"cell_type":ct,"gene":g,"SAH_mean":ms,"Sham_mean":mc,
                     "log2FC":np.log2((ms+1e-3)/(mc+1e-3)),"p_adj":pmap.get(g,np.nan)})
DE=pd.DataFrame(rows); DE.to_csv(os.path.join(RES,"scrna_celltype_DE_panels.csv"),index=False)
print("\ncell-type-specific axis-gene DE (SAH vs Sham):\n",DE.round(4).to_string(index=False))

# ---------- pseudotime on microglia/macrophage ----------
mg=adata[adata.obs.cell_type=="Microglia/Mac"].copy()
sc.pp.highly_variable_genes(mg,n_top_genes=2000,flavor='seurat')
mgs=mg[:,mg.var.highly_variable].copy()
sc.pp.pca(mgs,n_comps=15); sc.pp.neighbors(mgs,n_neighbors=15)
sc.tl.diffmap(mgs)
homeo=mgs.obs_names[mgs[:,"P2ry12"].X.toarray().ravel().argmax()] if "P2ry12" in mgs.var_names else mgs.obs_names[0]
mgs.uns['iroot']=list(mgs.obs_names).index(homeo)
sc.tl.dpt(mgs)
mgs.obs['dpt']=mgs.obs['dpt_pseudotime']
traj=[]
for g in ["Cd163","Hmox1","Fth1","Ftl1","Slc40a1","Hp","P2ry12","Csf1r"]:
    if g in mgs.var_names:
        e=np.asarray(mgs[:,g].X.todense()).ravel()
        rho,p=stats.spearmanr(mgs.obs.dpt.values, e)
        traj.append({"gene":g,"spearman_vs_pseudotime":rho,"p":p})
T=pd.DataFrame(traj); T.to_csv(os.path.join(RES,"scrna_pseudotime_correlation.csv"),index=False)
print("\nmicroglia->heme-processing pseudotime correlations:\n",T.round(4).to_string(index=False))

fig,axes=plt.subplots(1,3,figsize=(12,3.6))
dpt=mgs.obs.dpt.values
axes[0].hist(dpt,bins=30,color=OKABE["blue"]); axes[0].set_xlabel("DPT pseudotime"); axes[0].set_ylabel("#cells")
for ax,g,c in zip(axes[1:],["Cd163","Hmox1"],[OKABE["vermil"],OKABE["green"]]):
    if g in mgs.var_names:
        e=np.asarray(mgs[:,g].X.todense()).ravel()
        o=np.argsort(dpt); ax.scatter(dpt[o],e[o],s=4,c=c,alpha=.5,linewidths=0)
        ax.set_xlabel("pseudotime"); ax.set_ylabel(f"{g} (log-norm)"); ax.set_title(f"{g} vs pseudotime",fontsize=8)
fig.suptitle("Microglia/macrophage trajectory (homeostatic -> heme-processing)",fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig11_scrna_trajectory.png"),dpi=300,facecolor="white")
fig.savefig(os.path.join(FIG,"Fig11_scrna_trajectory.pdf")); plt.close(fig)

# ---------- SCENIC-lite curated regulon activity ----------
regulons={
 "NFE2L2_NRF2":["Hmox1","Nqo1","Gclm","Gclc","Slc7a11","Txnrd1","Ftl1","Fth1","Gsr","Prdx1","G6pd","Pgd","Me1"],
 "NFE2_erythroid":["Hbb-bs","Hba-a1","Alas2","Slc4a1","Ahsp","Gypa","Gypb","Blvra","Blvrb"],
 "BACH1":["Hmox1","Ftl1","Fth1","Bilz1","Slc40a1","Mt1","Mt2"],
}
tf_of={"NFE2L2_NRF2":"Nfe2l2","NFE2_erythroid":"Nfe2","BACH1":"Bach1"}
act=[]
for name,gl in regulons.items():
    gl=[g for g in gl if g in adata.var_names]
    sc.tl.score_genes(adata,gl,score_name=f"reg_{name}")
    tf=tf_of[name]; a=adata.obs
    if tf in adata.var_names:
        tfe=np.asarray(adata[:,tf].X.todense()).ravel()
        rho,p=stats.spearmanr(tfe, a[f"reg_{name}"].values)
    else: rho,p=np.nan,np.nan
    sah=a.loc[a.is_sah,f"reg_{name}"].mean(); sham=a.loc[~a.is_sah,f"reg_{name}"].mean()
    u,pp=stats.mannwhitneyu(a.loc[a.is_sah,f"reg_{name}"],a.loc[~a.is_sah,f"reg_{name}"])
    act.append({"regulon":name,"n_targets":len(gl),"TF":tf,
                "TF_vs_regulon_spearman":rho,"p":p,"SAH_activity":sah,"Sham_activity":sham,"SAH_vs_Sham_p":pp})
REG=pd.DataFrame(act); REG.to_csv(os.path.join(RES,"scrna_regulon_proxy.csv"),index=False)
print("\nSCENIC-lite curated regulon activity (proxy, NOT pySCENIC):\n",REG.round(4).to_string(index=False))
print("DONE 16_scrna_deep")
