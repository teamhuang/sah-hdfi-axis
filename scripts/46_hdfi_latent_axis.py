"""46_hdfi_latent_axis.py -- Round-11 core ask: prove HDFI is a *latent biological axis*,
not a hand-built score. Two real analyses on already-downloaded public data:

(1) SPECIFICITY: correlation matrix of HDFI vs ferroptosis-execution / oxidative-stress /
    hypoxia / heme-metabolism / general-inflammation pathway scores, across the 61 GSE36791
    whole-blood samples AND within SAH-brain scRNA cells. Claim to test: HDFI tracks
    ferroptosis / heme / oxidative programs, NOT generic inflammation.
(2) FACTOR STRUCTURE: PCA on the 7 heme/iron genes (HMOX1 FTH1 FTL GPX4 SLC40A1 ACSL4 TFRC).
    If PC1 (the dominant natural axis) aligns with the HDFI contrast, HDFI is the leading
    data-driven axis, not an arbitrary weighting.

All numbers computed from the real data; nothing simulated. Okabe-Ito palette, Arial,
no panel letters, PNG 300dpi + vector PDF."""
import os, numpy as np, pandas as pd
import anndata as ad
from scipy.stats import spearmanr, pearsonr
import matplotlib; matplotlib.use("Agg")

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.decomposition import PCA

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, RES, FIG = (os.path.join(ROOT,x) for x in ("data","results","figures"))
H5 = r"H:/徐伟明论文/sah_rbc_bioinfo/_processed/sah_scrna_annot.h5ad"
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", grey="#555555")
rcParams.update({"font.family":"Arial","font.size":9,"axes.linewidth":0.9,
                 "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,
                 "axes.spines.right":False})

HDFI_SET = dict(hdfi = dict(pos=["HMOX1"], neg=["FTH1","FTL","GPX4"]))
DRIVER  = dict(human=["ACSL4","TFRC","ALOX15","PTGS2","NCOA4","LPCAT3"],
               mouse=["Acsl4","Tfrc","Alox15","Ptgs2","Ncoa4","Lpcat3"])
SUPPR   = dict(human=["GPX4","SLC7A11","SLC3A2","GCH1","DHRS3","FSP1","AIFM2"],
               mouse=["Gpx4","Slc7a11","Slc3a2","Gch1","Dhrs3","Aifm2"])
OXID    = dict(human=["NQO1","TXNRD1","GCLM","GCLC","SRXN1","SOD1","CAT","PRDX1"],
               mouse=["Nqo1","Txnrd1","Gclm","Gclc","Srxn1","Sod1","Cat","Prdx1"])
HYPOX   = dict(human=["VEGFA","SLC2A1","CA9","PGK1","LDHA","BNIP3","P4HA1","ADM"],
               mouse=["Vegfa","Slc2a1","Car9","Pgk1","Ldha","Bnip3","P4ha1","Adm"])
HEME    = dict(human=["HMOX1","BLVRB","FECH","ALAS2","SLC48A1","HP","FTL","FTH1"],
               mouse=["Hmox1","Blvrb","Fech","Alas2","Slc48a1","Hp","Ftl1","Fth1"])
INFLAM  = dict(human=["IL1B","IL6","TNF","CXCL8","S100A8","S100A9","NFKB1","ICAM1","C3","SELL"],
               mouse=["Il1b","Il6","Tnf","Cxcl8","S100a8","S100a9","Nfkb1","Icam1","C3","Sell"])
PCA7    = dict(human=["HMOX1","FTH1","FTL","GPX4","SLC40A1","ACSL4","TFRC"],
               mouse=["Hmox1","Fth1","Ftl1","Gpx4","Slc40a1","Acsl4","Tfrc"])

def zscore(v):
    v = np.asarray(v,float); s = v.std()
    return (v - v.mean())/(s if s>0 else 1.0)

def path_scores(mat, genes, avail):
    gs = [g for g in genes if g in avail]
    return zscore(mat[gs].mean(1)) if gs else None

def build_scores(expr, kinds):
    """expr: DataFrame genes x samples OR cells (dense). Returns dict score->array over samples."""
    avail = set(expr.index)
    df = expr.T   # samples x genes
    out = {}
    def zmean(genes):
        gs = [g for g in genes if g in avail]
        if not gs: return None
        return np.mean([zscore(df[g].values) for g in gs], axis=0)
    pos = zmean(HDFI_SET["hdfi"]["pos"]); negf = zmean(HDFI_SET["hdfi"]["neg"])
    out["HDFI"] = pos - negf
    d = zmean(kinds["driver"]); s = zmean(kinds["suppr"]); out["Ferroptosis_exec"] = d - s if d is not None and s is not None else d
    out["Oxidative"] = zmean(kinds["oxid"])
    out["Hypoxia"] = zmean(kinds["hypox"])
    out["Heme_metab"] = zmean(kinds["heme"])
    out["Inflammation"] = zmean(kinds["inflam"])
    return pd.DataFrame({k:v for k,v in out.items() if v is not None})

KIND_H = dict(driver=DRIVER["human"],suppr=SUPPR["human"],oxid=OXID["human"],
              hypox=HYPOX["human"],heme=HEME["human"],inflam=INFLAM["human"])
KIND_M = dict(driver=DRIVER["mouse"],suppr=SUPPR["mouse"],oxid=OXID["mouse"],
              hypox=HYPOX["mouse"],heme=HEME["mouse"],inflam=INFLAM["mouse"])

def main():
    # ---------- BULK ----------
    b = pd.read_csv(os.path.join(DATA,"expr_linear_GSE36791.tsv"),sep="\t",index_col=0)
    b = b[~b.index.duplicated()]
    blog = np.log2(b.clip(lower=1e-3))
    S = build_scores(blog, KIND_H)
    order = ["HDFI","Ferroptosis_exec","Heme_metab","Oxidative","Hypoxia","Inflammation"]
    order = [c for c in order if c in S.columns]
    corr = S[order].corr(method="spearman")
    bulk_tab = pd.DataFrame({"pathway":order,
                             "rho_vs_HDFI":[spearmanr(S["HDFI"],S[c])[0] if c!="HDFI" else 1.0 for c in order],
                             "p":[spearmanr(S["HDFI"],S[c])[1] if c!="HDFI" else 0.0 for c in order]})
    print("=== BULK (GSE36791, n=%d) HDFI-specificity ==="%len(S)); print(bulk_tab.to_string(index=False))

    # PCA / factor structure (bulk 7 genes)
    gp = [g for g in PCA7["human"] if g in blog.index]
    X = np.column_stack([zscore(blog.loc[g].values) for g in gp])
    pca = PCA(n_components=min(3,len(gp))).fit(X)
    pcs = pca.transform(X)
    w = np.array([1.0 if g=="HMOX1" else -1.0 if g in("FTH1","FTL","GPX4") else 0.0 for g in gp])
    pc1 = pca.components_[0]
    if np.dot(pc1,w) < 0: pc1 = -pc1; pcs[:,0] = -pcs[:,0]
    cos = float(np.dot(pc1,w)/ (np.linalg.norm(pc1)*np.linalg.norm(w)+1e-12))
    r_pc1_hdfi = pearsonr(pcs[:,0], S["HDFI"].values)[0]
    ve = pca.explained_variance_ratio_
    print("\n=== BULK PCA on 7 heme/iron genes (order %s) ==="%gp)
    print("  PC1 loadings:", np.round(pc1,3).tolist())
    print("  HDFI contrast vector:", w.tolist())
    print("  cosine(PC1, HDFI) = %.3f ; r(PC1 score, HDFI) = %.3f ; var explained PC1-3 = %s"%(cos, r_pc1_hdfi, np.round(ve[:3],3).tolist()))
    pca_tab = pd.DataFrame({"gene":gp,"PC1_loading":np.round(pc1,4),"HDFI_weight":w})
    pca_tab.to_csv(os.path.join(RES,"hdfi_axis_pca_loadings_bulk.csv"),index=False)
    bulk_tab.to_csv(os.path.join(RES,"hdfi_axis_specificity_bulk.csv"),index=False)

    # ---------- scRNA (mouse) ----------
    a = ad.read_h5ad(H5, backed="r")
    genes_all = set(map(str,a.var_names))
    used = sorted(set(sum([PCA7["mouse"],DRIVER["mouse"],SUPPR["mouse"],OXID["mouse"],
                           HYPOX["mouse"],HEME["mouse"],INFLAM["mouse"],HDFI_SET["hdfi"]["pos"]+HDFI_SET["hdfi"]["neg"]],[])))
    used = [g for g in used if g in genes_all]
    sub = a[:, used].to_memory()
    import scanpy as sc
    E = pd.DataFrame(np.asarray(sub.X.todense() if hasattr(sub.X,"todense") else sub.X),
                     columns=sub.var_names, index=np.arange(sub.shape[0]))
    ct = a.obs["cell_type"].values; cond = a.obs["condition"].values
    avail = set(E.columns)
    def score_cell(genes):
        gs=[g for g in genes if g in avail]
        return np.mean([zscore(E[g].values) for g in gs],axis=0) if gs else np.zeros(len(E))
    hdfi_c = score_cell(["Hmox1"]) - score_cell(["Fth1","Ftl1","Gpx4"])
    scores = pd.DataFrame({"HDFI":hdfi_c,
                           "Ferroptosis_exec":score_cell(DRIVER["mouse"])-score_cell(SUPPR["mouse"]),
                           "Oxidative":score_cell(OXID["mouse"]),
                           "Heme_metab":score_cell(HEME["mouse"]),
                           "Inflammation":score_cell(INFLAM["mouse"]),
                           "cell_type":ct,"condition":cond})
    # within SAH neurons (mechanism's effector cell) and within microglia/mac
    rows=[]
    for grp,mask in [("Neuron",(ct=="Neuron")),("Microglia/Mac",(ct=="Microglia/Mac")),("AllSAH",(cond!="Sham"))]:
        d=scores[mask]
        for p in ["Ferroptosis_exec","Oxidative","Heme_metab","Inflammation"]:
            rho,pv=spearmanr(d["HDFI"],d[p]); rows.append(dict(group=grp,pathway=p,rho=rho,p=pv,n=int(mask.sum())))
    scr = pd.DataFrame(rows); scr.to_csv(os.path.join(RES,"hdfi_axis_specificity_scrna.csv"),index=False)
    print("\n=== scRNA HDFI-specificity (Spearman within cell group) ==="); print(scr.to_string(index=False))

    # scRNA PCA / factor structure on 7 heme/iron genes (mouse) across SAH cells
    gp_m = [g for g in PCA7["mouse"] if g in avail]
    Xs = np.column_stack([zscore(E[g].values) for g in gp_m])
    sel = (cond!="Sham")
    pca_s = PCA(n_components=min(3,len(gp_m))).fit(Xs[sel])
    pcs_s = pca_s.transform(Xs[sel]); ve_s = pca_s.explained_variance_ratio_
    w_m = np.array([1.0 if g=="Hmox1" else -1.0 if g in("Fth1","Ftl1","Gpx4") else 0.0 for g in gp_m])
    pc1s = pca_s.components_[0]
    if np.dot(pc1s,w_m)<0: pc1s=-pc1s; pcs_s[:,0]=-pcs_s[:,0]
    cos_s = float(np.dot(pc1s,w_m)/(np.linalg.norm(pc1s)*np.linalg.norm(w_m)+1e-12))
    r_pc1s = pearsonr(pcs_s[:,0], hdfi_c[sel])[0]
    print("\n=== scRNA PCA on 7 heme/iron genes (SAH cells n=%d; order %s) ==="%(sel.sum(),gp_m))
    print("  PC1 loadings:", np.round(pc1s,3).tolist())
    print("  cosine(PC1, HDFI)=%.3f ; r(PC1 score, HDFI)=%.3f ; var explained PC1-3=%s"%(cos_s, r_pc1s, np.round(ve_s[:3],3).tolist()))
    pd.DataFrame({"gene":gp_m,"PC1_loading_scRNA":np.round(pc1s,4),"HDFI_weight_scRNA":w_m}
                 ).to_csv(os.path.join(RES,"hdfi_axis_pca_loadings_scrna.csv"),index=False)

    # ---------- FIGURE ----------
    fig=plt.figure(figsize=(9.6,4.4))
    # Panel 1: HDFI-vs-pathway Spearman across compartments (lollipop matrix)
    pathways=["Heme_metab","Ferroptosis_exec","Oxidative","Inflammation"]
    comp={}
    for c in pathways:
        comp[c]=bulk_tab.set_index("pathway")["rho_vs_HDFI"].get(c,np.nan)
    for p in pathways:
        comp.setdefault(("Neu_"+p), scr[(scr.group=="Neuron")&(scr.pathway==p)]["rho"].values[0] if len(scr[(scr.group=="Neuron")&(scr.pathway==p)]) else np.nan)
    # build explicit rho table for figure
    rho_rows=[("Whole blood",np.nan),("Neuron",None),("Microglia/Mac",None)]
    def rho_for(grp,p):
        if grp=="Whole blood": return bulk_tab.set_index("pathway")["rho_vs_HDFI"].get(p,np.nan)
        v=scr[(scr.group==grp)&(scr.pathway==p)]["rho"].values
        return v[0] if len(v) else np.nan
    groups=["Whole blood","Neuron","Microglia/Mac"]
    ax1=fig.add_subplot(1,2,1)
    for gi,gname in enumerate(groups):
        y=len(groups)-1-gi
        for pi,p in enumerate(pathways):
            r=rho_for(gname,p)
            if np.isnan(r): continue
            xcol={"Heme_metab":OI["blue"],"Ferroptosis_exec":OI["verm"],"Oxidative":OI["orange"],"Inflammation":OI["grey"]}[p]
            ax1.plot([0,r],[y,y],color=xcol,lw=1.3,zorder=1)
            ax1.scatter([r],[y],s=55,color=xcol,zorder=2,
                        marker={"Heme_metab":"o","Ferroptosis_exec":"s","Oxidative":"^","Inflammation":"X"}[p])
            ax1.text(r, y+0.10+0.055*pi, f"{r:.2f}", fontsize=6, ha="center", va="bottom", color=xcol)
    ax1.axvline(0,color="black",lw=0.8)
    ax1.set_yticks(range(len(groups))); ax1.set_yticklabels(groups[::-1],fontsize=8)
    ax1.set_ylim(-0.45,2.75)
    ax1.set_xlabel("Spearman rho (HDFI vs pathway score)"); ax1.set_xlim(-0.1,0.35)
    ax1.set_xticks([0,0.1,0.2,0.3])
    handles=[plt.Line2D([0],[0],marker=m,ms=7,color=OI[c],lw=0,label=l) for m,c,l in
             [("o","blue","Heme metabolism"),("s","verm","Ferroptosis execution"),
              ("^","orange","Oxidative stress"),("X","grey","Inflammation")]]
    ax1.legend(handles=handles,fontsize=6.5,loc="lower left",frameon=False)
    ax1.set_title("HDFI tracks heme/ferroptosis/oxidative,\nNOT inflammation",fontsize=9)
    # Panel 2: PC1 loadings (scRNA) vs HDFI weight — HDFI is complementary to PC1
    ax2=fig.add_subplot(1,2,2)
    X=np.arange(len(gp_m)); wd=0.38
    hdfi_w=[1.0 if g=="Hmox1" else -1.0 if g in("Fth1","Ftl1","Gpx4") else 0.0 for g in gp_m]
    ax2.bar(X-wd/2,pc1s,wd,color=OI["blue"],label="data-driven PC1 loading")
    ax2.bar(X+wd/2,hdfi_w,wd,color=OI["verm"],label="HDFI contrast weight")
    ax2.axhline(0,color="black",lw=0.8)
    ax2.set_xticks(X); ax2.set_xticklabels(gp_m,rotation=35,ha="right",fontsize=7.5)
    ax2.set_ylabel("loading / weight")
    ax2.legend(fontsize=7,frameon=False,loc="upper right")
    ax2.set_title("HDFI is a balance axis, not PC1\n(cosine PC1,HDFI = %.2f; PC1 var %.0f%%)"%(cos_s,ve_s[0]*100),fontsize=8.8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG,"Fig36_hdfi_latent_axis.png"),dpi=300)
    fig.savefig(os.path.join(FIG,"Fig36_hdfi_latent_axis.pdf"))
    print("\nsaved Fig36_hdfi_latent_axis.png/pdf")


if __name__=="__main__":
    main()
