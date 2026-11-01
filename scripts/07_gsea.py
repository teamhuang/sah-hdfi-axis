"""07_gsea.py -- GSEA (gseapy prerank) on signed -log10 P ranking of GSE36791 SAH-vs-Control,
against Hallmark + KEGG + Reactome. Reports pathways of interest: HIF-1, ferroptosis, NF-kB,
complement, oxidative stress, hypoxia, inflammation, heme/erythrocyte."""
import os, numpy as np, pandas as pd, gseapy
import logging; logging.disable(logging.WARNING)
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT=os.path.join(BASE,"results"),os.path.join(BASE,"data")

def read_gmt(path):
    d={}
    with open(path,encoding="utf-8") as f:
        for line in f:
            p=line.rstrip("\n").split("\t")
            g=[x for x in p[2:] if x]
            if g: d[p[0]]=set(g)
    return d

sets={}
for fn,pre in [("MSigDB_Hallmark_2020.gmt","H:"),("KEGG_2021_Human.gmt","K:"),("Reactome_2022.gmt","R:")]:
    for k,v in read_gmt(os.path.join(DAT,"gmt",fn)).items(): sets[pre+k]=v

deg=pd.read_csv(os.path.join(RES,"DEG_GSE36791.csv"))
deg=deg.dropna(subset=["logFC"])
stat=np.sign(deg.logFC.values)*(-np.log10(np.clip(deg["P.Value"].values,1e-300,1)))
rnk=pd.Series(stat,index=deg.gene.values).sort_values(ascending=False)
print("ranking range:", round(rnk.min(),2), "to", round(rnk.max(),2))

pre=gseapy.prerank(rnk=rnk, gene_sets=sets, min_size=15, max_size=500,
                   permutation_num=2000, weight=1.0, threads=4, outdir=None, seed=2026, no_plot=True)
res=pre.res2d.copy()
res.columns=[c.strip() for c in res.columns]
res["NES"]=pd.to_numeric(res["NES"]); res["FDR q-val"]=pd.to_numeric(res["FDR q-val"],errors="coerce")
res=res.sort_values("NES")
res.to_csv(os.path.join(RES,"GSEA_GSE36791_all.csv"),index=False)

watch=["HIF-1","Hypoxia","Ferroptosis","NF-kB","NF kappa","Complement","Oxidative","Reactive oxygen",
       "heme","Erythrocyte","Hemoglobin","Inflammation","TNF","IL-6","JAK","STAT","Apoptosis",
       "Necroptosis","Phagosome","Iron"]
def show(sub,label):
    print(f"\n=== {label} ===")
    cols=["Term","NES","p-value","FDR q-val","NOM p-val","Leading edge"]
    print(sub[[c for c in cols if c in sub.columns]].to_string(index=False,max_colwidth=48))

mask=res.Term.str.contains("|".join(watch),case=False,regex=True)
show(res[mask],"PATHWAYS OF INTEREST")
print("\n=== TOP 10 UP-in-SAH (positive NES) ===");show(res[res.NES>0].tail(10),"UP NES")
print("\n=== TOP 10 DOWN-in-SAH (negative NES) ===");show(res[res.NES<0].head(10),"DOWN NES")

# significance summary
sig=res[(res["FDR q-val"]<0.25)]
print(f"\nsets with FDR<0.25: {len(sig)} of {len(res)}")
pd.DataFrame({"pathway":["HIF-1","Ferroptosis","NF-kB","Complement","Oxidative/hypoxia","Heme/Erythrocyte"],
              "best_NES":[res[res.Term.str.contains(t,case=False)].NES.abs().max() if (res.Term.str.contains(t,case=False)).any() else np.nan
                          for t in ["HIF-1","Ferroptosis","NF-kB|kappa B","Complement","Oxidative|Hypoxia|Reactive oxygen","heme|Erythrocyte|Hemoglobin"]]}
            ).to_csv(os.path.join(RES,"GSEA_pathway_focus.csv"),index=False)
print("DONE 07_gsea")
