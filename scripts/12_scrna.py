"""12_scrna.py -- mouse brain SAH single-cell atlas (GSE311852): load 10x filtered matrices,
QC, integrate, cluster, annotate, heme/iron module score by cell type (SAH vs Sham),
HP/CD163/HMOX1/FTH1 localisation, and a CellChat-style HP-Cd163 communication analysis."""
import os, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, gzip
import scanpy as sc, anndata as ad
sc.settings.verbosity=1
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S=os.path.join(BASE,"scrna"); RES=os.path.join(BASE,"results"); FIG=os.path.join(BASE,"figures")
os.makedirs(RES,exist_ok=True)
SCRATCH="H:/徐伟明论文/sah_rbc_bioinfo/_processed"
os.makedirs(SCRATCH,exist_ok=True)

SAMPLES={ "GSM9333532":"Sham","GSM9333533":"SAH_24h","GSM9333534":"SAH_24h",
          "GSM9333535":"SAH_72h","GSM9333536":"SAH_72h"}
tag2file={"Sham":"filtered_Sham","SAH_24h":None,"SAH_72h":None}

def read_sample(gsm):
    stem=SAMPLES[gsm]
    import glob
    mtx=glob.glob(os.path.join(S,f"{gsm}_*matrix.mtx.gz"))[0]
    bc=glob.glob(os.path.join(S,f"{gsm}_*barcodes.tsv.gz"))[0]
    ft=glob.glob(os.path.join(S,f"{gsm}_*features.tsv.gz"))[0]
    with gzip.open(mtx,'rt') as fh:
        M=sio.mmread(fh)
    M=sp.csr_matrix(M)
    barcodes=[l.strip() for l in gzip.open(bc,'rt')]
    feat=pd.read_csv(ft,sep="\t",header=None)
    sym=feat[1].astype(str).values
    A=ad.AnnData(X=M.T.tocsr())
    A.obs['gsm']=gsm; A.obs['sample']=stem
    A.obs['condition']=('Sham' if stem=='Sham' else stem)
    A.var['gene_symbols']=sym
    A.var_names=sym
    A.var_names_make_unique()
    A.obs['orig']=A.obs_names.map(lambda x:f"{gsm}_{x}")
    return A

ads=[read_sample(g) for g in SAMPLES]
adata=ad.concat(ads, index_unique='_')
print("merged:",adata.shape)
# QC
adata.var['mt']=adata.var_names.str.startswith('mt-')
adata.var['ribo']=adata.var_names.str.startswith(('Rps','Rpl'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt','ribo'],percent_top=None, inplace=True)
sc.pp.filter_cells(adata,min_genes=200)
adata=adata[adata.obs.pct_counts_mt<20]
print("after QC:",adata.shape)
# subsample to keep runtime bounded
if adata.n_obs>60000:
    adata=adata[sc.pp.subsample(adata,n_obs=60000,random_state=2026,return_subset=True)]
print("analysis cells:",adata.shape)
print("cells per sample/condition:\n",adata.obs.groupby(['sample','condition']).size())

sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2500,flavor='seurat')
hvg=adata[:,adata.var.highly_variable].copy()
sc.pp.scale(hvg,max_value=10); sc.tl.pca(hvg,n_comps=40)
sc.pp.neighbors(hvg,n_neighbors=15)
try:
    sc.tl.leiden(hvg,resolution=1.0,key_added='leiden',flavor='igraph',n_iterations=2,directed=False)
except Exception as e:
    print("leiden igraph fallback:",repr(e)[:80]); sc.tl.leiden(hvg,resolution=1.0,key_added='leiden')
adata.obs['leiden']=hvg.obs['leiden'].values
sc.tl.umap(hvg)
adata.obsm['X_umap']=hvg.obsm['X_umap']

# ---- annotation via canonical mouse markers ----
markers={
 'Microglia/Mac':['Cx3cr1','P2ry12','C1qa','C1qb','Ctss','Adgre1','Cd68','Itgam','Csf1r','Mrc1','Cd163','Spp1','Clec7a','Lpl','Trem2','ApoE'],
 'Neutrophil':['Elane','Prtn3','Mpo','S100a8','S100a9','Csfr3','Ly6g','Retnlg','Camp','Ltf'],
 'T cell':['Cd3e','Cd3d','Cd8a','Cd4','Lck','Il7r','Ccr7','Cd27','Il2rb'],
 'B cell':['Cd79a','Cd79b','Ms4a1','Ebf1','Pax5','IGHM'],
 'NK/GM':['Klrd1','Gzma','Gzmb','Ncr1','Ifng'],
 'Endothelial':['Cldn5','Pecam1','Flt1','Slco1a4','Aqp4_low'],
 'Astrocyte':['Aqp4','Gfap','Slc1a3','Aldh1l1','Ttyl1'],
 'Oligodendrocyte':['Mbp','Plp1','Sox10','Mobp','Mal','Cnp'],
 'OPC':['Pdgfra','Cspg4','Sox10_low'],
 'Neuron':['Snap25','Rbfox3','Syt1','Nefl','Gad1','Slc17a7'],
 'Erythroid/RBC':['Hba-a1','Hba-a2','Hbb-b1','Hbb-bs','Alas2','Slc4a1'],
 'Mast':['Cpa3','Kit','Ms4a2','Hs3st1'],
}
score_df={}
for ct,gl in markers.items():
    gl=[g for g in gl if g in adata.var_names]
    if gl: sc.tl.score_genes(adata,gl,score_name=f"s_{ct}")
    else: adata.obs[f"s_{ct}"]=0
score_df=pd.DataFrame({ct:adata.obs[f"s_{ct}"] for ct in markers})
# z-score each marker score across cells so module-size doesn't bias the call, then argmax
score_z=(score_df-score_df.mean())/ (score_df.std()+1e-9)
adata.obs['cell_type']=score_z.idxmax(axis=1)
# drop implausible placeholder keys
print("cell type counts:\n",adata.obs.cell_type.value_counts())
sc.pl.umap(adata,color=['cell_type','condition'],save="_celltypes.png",show=False,
           legend_loc='right margin') if False else None

# ---- heme/iron module score (AddModuleScore equivalent) ----
heme=["Hp","Hmox1","Fth1","Ftl1","Slc40a1","Slc25a37","Nfe2","Alas1","Alas2","Cd163","Mrc1","Ftmt","Blvrb","Slc4a1"]
heme=[g for g in heme if g in adata.var_names]
print("heme module genes present:",heme)
sc.tl.score_genes(adata,heme,score_name='heme_score')
adata.obs['is_sah']=~adata.obs.condition.eq('Sham')

res=[]
for ct in adata.obs.cell_type.cat.categories if hasattr(adata.obs.cell_type,'cat') else sorted(adata.obs.cell_type.unique()):
    d=adata[adata.obs.cell_type==ct]
    a=d.obs.loc[d.obs.is_sah,"heme_score"]; b=d.obs.loc[~d.obs.is_sah,"heme_score"]
    from scipy import stats
    if len(a)>10 and len(b)>10:
        t,p=stats.mannwhitneyu(a,b)
        res.append({"cell_type":ct,"n_sah":len(a),"n_sham":len(b),
                    "heme_sah":a.mean(),"heme_sham":b.mean(),"delta":a.mean()-b.mean(),
                    "MWU_p":p})
R=pd.DataFrame(res).sort_values("delta",ascending=False)
R.to_csv(os.path.join(RES,"scrna_heme_score_by_celltype.csv"),index=False)
print("\nheme/iron module score, SAH vs Sham by cell type:\n",R.round(4).to_string(index=False))

# marker expression by cell type (mean, pseudobulk-ish) -- restrict to needed genes (memory-safe)
mk=[g for g in ['Hp','Cd163','Mrc1','Hmox1','Fth1','Ftl1','Slc40a1','Nfe2','Alas2','Hba-a1','Elane','Cx3cr1','P2ry12'] if g in adata.var_names]
lr_genes=['Hp','Cd163','Mrc1','Hmox1','CD91','Tfrc']
needed=sorted(set(mk+heme+lr_genes))
needed=[g for g in needed if g in adata.var_names]
expr=pd.DataFrame(adata[:,needed].X.toarray(),index=adata.obs_names,columns=needed)
mean_expr=pd.DataFrame({ct: expr.loc[adata.obs.cell_type==ct, mk].mean() for ct in adata.obs.cell_type.unique()}).T
mean_expr.to_csv(os.path.join(RES,"scrna_marker_meanexpr.csv"))
print("\nmean normalised expression of axis markers by cell type:\n",mean_expr.round(2).to_string())

# composition SAH vs Sham
comp=pd.crosstab(adata.obs.cell_type, adata.obs.condition, normalize='columns')
comp.to_csv(os.path.join(RES,"scrna_celltype_composition.csv"))

# ---- CellChat-style HP-Cd163 communication ----
lr_pairs=[("Hp","Cd163"),("Hp","Mrc1"),("Hmox1","CD91"),("Fth1","Tfrc")]
comm=[]
ctypes=sorted(adata.obs.cell_type.unique())
for L,Rc in lr_pairs:
    if L not in adata.var_names or Rc not in adata.var_names: continue
    lig=expr.groupby(adata.obs.cell_type)[L].mean()
    rec=expr.groupby(adata.obs.cell_type)[Rc].mean()
    # pairwise communication strength = mean(lig in src)*mean(rec in tgt), scaled
    for src in ctypes:
        for tgt in ctypes:
            if src==tgt: continue
            strength=lig.get(src,0)*rec.get(tgt,0)
            comm.append({"pair":f"{L}-{Rc}","source":src,"target":tgt,"strength":strength})
C=pd.DataFrame(comm).sort_values("strength",ascending=False)
C.to_csv(os.path.join(RES,"scrna_LR_communication.csv"),index=False)
print("\ntop HP-Cd163 / heme ligand-receptor communication edges:\n",C[C.pair=='Hp-Cd163'].head(8).to_string(index=False))

adata.write(os.path.join(SCRATCH,"sah_scrna_annot.h5ad"))
# save a couple figures
sc.set_figure_params(dpi_save=200,frameon=False)
sc.pl.umap(adata,color=['cell_type'],save="_celltype.png",show=False)
sc.pl.umap(adata,color=['heme_score','Hp','Cd163','Hmox1'],ncols=2,save="_heme.png",show=False)
sc.pl.dotplot(adata,['Hp','Cd163','Mrc1','Hmox1','Fth1','Slc40a1','Nfe2','Elane','Cx3cr1','P2ry12','Cd3e','Mbp','Snap25'],
              groupby='cell_type',save="_markers.png",show=False)
print("DONE 12_scrna")
