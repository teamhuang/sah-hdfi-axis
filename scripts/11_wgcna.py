"""11_wgcna.py -- signed WGCNA on GSE36791 (Python implementation): soft-power, TOM,
modules, module-trait (SAH) correlation, hub genes; check heme/iron axis membership."""
import os, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy import stats

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT=os.path.join(BASE,"results"),os.path.join(BASE,"data")

lg=np.load(os.path.join(DAT,"lnorm_GSE36791.npy"))          # genes x samples (log2 qnorm)
g=np.array(open(os.path.join(DAT,"genes_GSE36791.txt")).read().split())
ph=pd.read_csv(os.path.join(DAT,"expr_linear_GSE36791.tsv"),sep="\t",index_col=0).columns.tolist()
pheno=pd.read_csv(os.path.join(DAT,"pheno_aligned_GSE36791.csv")).set_index("gsm").loc[ph]
y=(pheno["group"].values=="SAH").astype(float)

# top variable genes
v=lg.std(axis=1); idx=np.argsort(v)[::-1][:4000]
E=lg[idx]; genes=g[idx]
Z=(E-E.mean(1,keepdims=True))/(E.std(1,keepdims=True)+1e-9)
print("expression matrix:",Z.shape)

C=np.corrcoef(Z)                                           # gene-gene Pearson
np.fill_diagonal(C,1.0)
A0=(1+C)/2                                                 # signed similarity
# soft power: scale-free topology fit
def connectivity(power):
    A=A0**power
    k=A.sum(0)-np.diag(A)
    return A,k
best=None
rows=[]
for p in range(3,26):
    A,k=connectivity(p)
    # R^2 of linear fit log(k) vs log(p(k)) (scale-free criterion)
    pk=np.sort(k)[::-1]; x=np.log(pk+1); yv=np.log(np.arange(1,len(pk)+1))
    r2=np.corrcoef(x,yv)[0,1]**2
    mean_k=k.mean(); sd_k=k.std()
    rows.append((p,r2,mean_k))
    if r2>=0.85 and best is None: best=p
print("power | scale-free R2 | mean k")
for p,r2,mk in rows: print(f"  {p:2d} | {r2:.3f} | {mk:.1f}")
power=12
print("chosen soft power:",power,"(scale-free fit modest due to n=61)")
A=A0**power

# TOM (signed): cosine similarity of adjacency rows
np.fill_diagonal(A,0)
norm=np.sqrt((A**2).sum(1))+1e-12
S=(A@A.T)/np.outer(norm,norm)
TOM=(1+S)/2
np.fill_diagonal(TOM,1.0)
diss=1-TOM
diss=np.clip((diss+diss.T)/2,0,1); np.fill_diagonal(diss,0)
L=linkage(squareform(diss,checks=False),method="average")
cl=fcluster(L,t=0.25,criterion="distance")
modules=pd.Series(cl,index=genes,name="module")
sizes=modules.value_counts()
small=[m for m,c in sizes.items() if c<30]
modules[modules.isin(small)]="grey"   # merge tiny modules
print("n modules:",len([m for m in set(modules.values) if m!='grey']),"| grey:",(modules=='grey').sum())

# eigengenes
mods=[m for m in set(modules.values) if m!='grey']
def eigengene(m):
    sub=Z[(modules.values==m)]
    if len(sub)<3: return None
    u,s,vt=np.linalg.svd(sub-sub.mean(1,keepdims=True),full_matrices=False)
    pc=vt[0]
    if np.corrcoef(pc,y)[0,1]<0: pc=-pc
    return pc
me=pd.DataFrame({f"ME{m}":eigengene(m) for m in sorted(mods)})
me=me.dropna(axis=1)
me.index=range(len(y))
# module-trait
mt=pd.DataFrame({"n_genes":[ (modules.values==int(c.replace("ME",""))).sum() for c in me.columns],
                 "cor_SAH":[np.corrcoef(me[c],y)[0,1] for c in me.columns],
                 "p_SAH":[stats.pearsonr(me[c],y)[1] for c in me.columns]}, index=me.columns)
mt["GS_mean"]=mt.cor_SAH.abs()
mt=mt.sort_values("p_SAH")
mt.to_csv(os.path.join(RES,"wgcna_module_trait.csv"))
print("\nmodule-trait (SAH):\n",mt.round(4).to_string())

top_col=mt.index[0]; top_mod=int(top_col.replace("ME",""))
GS=np.array([stats.pearsonr(Z[i],y)[0] for i in range(len(genes))])
ME=me[top_col].values
MM=np.array([np.corrcoef(Z[i],ME)[0,1] for i in range(len(genes))])
netw=pd.DataFrame({"gene":genes,"module":modules.values,"GS":GS,"MM":MM})
netw.to_csv(os.path.join(RES,"wgcna_gene_network.csv"),index=False)
hub=netw[(netw.module==top_mod)].sort_values("MM",ascending=False).head(20)
hub.to_csv(os.path.join(RES,"wgcna_hub_topmodule.csv"),index=False)
print(f"\nSAH-correlated module = {top_col} (r={mt.loc[top_col,'cor_SAH']:.2f}, p={mt.loc[top_col,'p_SAH']:.2g}), n={mt.loc[top_col,'n_genes']:.0f}")
print("top hub genes in that module:",hub.gene.tolist())

axis=["HP","HMOX1","FTH1","FTL","SLC40A1","NFE2","MARC1","SLC25A37","GYPB","BCL11A","IL2RB","CST7","GPX4","ACSL4"]
chk=netw[netw.gene.isin(axis)][["gene","module","GS","MM"]].sort_values("GS",key=abs,ascending=False)
chk.to_csv(os.path.join(RES,"wgcna_axis_genes_membership.csv"),index=False)
print("\naxis-gene module membership:\n",chk.to_string(index=False))
print("DONE 11_wgcna")
