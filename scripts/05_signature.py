"""05_signature.py -- build HRR / red-cell-dysfunction gene signature, intersect with SAH DEGs,
compute program-activity association with SAH in both cohorts, and ORA on candidate genes."""
import os, numpy as np, pandas as pd
from scipy import stats
import gseapy

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT,RES=os.path.join(BASE,"data"),os.path.join(BASE,"results")

def read_gmt(path):
    d={}
    with open(path,encoding="utf-8") as f:
        for line in f:
            parts=line.rstrip("\n").split("\t")
            d[parts[0]]=set(x for x in parts[2:] if x)
    return d

hall=read_gmt(os.path.join(DAT,"gmt","MSigDB_Hallmark_2020.gmt"))
heme=hall.get("HALLMARK_HEME_METABOLISM", set())
print("HALLMARK_HEME_METABOLISM size:", len(heme))

# ---- curated HRR / erythroid sub-modules (from the two HRR components) ----
curated = {
 "Hb_globin_oxygen_carrying":["HBB","HBA1","HBA2","HBD","HBG1","HBG2","HBE1","HBZ","AHSP"],
 "heme_biosynthesis":["ALAS2","ALAD","FECH","UROS","UROD","CPOX","PPOX","BLVRA","BLVRB","ABCB10","GLRX5","SLC25A37","SLC25A38","HMBS"],
 "iron_handling_oxidative":["TFRC","TF","TFR2","SLC40A1","FTH1","FTL","ACO1","HMOX1","HMOX2","HP","HPX","CP","FTMT"],
 "RBC_membrane_structure":["SLC4A1","EPB42","ANK1","SPTA1","SPTB","GYPA","GYPB","GYPC","RHAG","RHCE","KEL","ERMAP","ADD2","LMNA","CA1","CA2","GATA1"],
 "erythropoiesis_maturation":["GATA1","KLF1","EPOR","MYB","BCL11A","NFE2","ANK2","TFRC","HRI","EIF2AK1","TRNT1","SEC14L5","TNKS"],
}
for k,v in curated.items(): print(f"curated[{k}]={len(v)}")

# ---- GO Biological Process sets via Enrichr (erythroid-specific only) ----
go_terms_want=["oxygen transport","erythrocyte differentiation","erythrocyte maturation",
               "erythrocyte membrane","heme metabolism","heme biosynthesis","erythropoiesis",
               "hemoglobin biosynthesis","hemoglobin binding","heme binding","iron ion transport"]
cache=os.path.join(DAT,"gmt","GO_BP_erythroid_cache.csv")
gobp={}
if os.path.exists(cache):
    df=pd.read_csv(cache)
    for s,gr in df.groupby("set"): gobp[s]=set(gr["gene"])
    print("loaded cached GO erythroid sets:",len(gobp))
else:
    lib=gseapy.get_library(name="GO_Biological_Process_2023", organism="human")
    for setname, genes in lib.items():
        nm=setname.lower()
        if any(t.lower() in nm for t in go_terms_want):
            gobp[setname]=set(genes)
    print("GO BP erythroid matched sets:", len(gobp))
    pd.DataFrame([(s,g) for s,v in gobp.items() for g in v],columns=["set","gene"]).to_csv(cache,index=False)
for k in list(gobp)[:12]: print("  GO:",k, len(gobp[k]))

# ---- union signature ----
sig=set(heme)
for v in curated.values(): sig|=set(v)
for v in gobp.values(): sig|=v
universe=set(pd.read_csv(os.path.join(DAT,"expr_linear_GSE36791.tsv"),sep="\t",index_col=0).index)
sig_all=sig
sig=sorted(sig & universe)
print("\nRBC/HRR signature (union):", len(sig_all), "-> in measurable universe:", len(sig))
open(os.path.join(RES,"RBC_signature_genes.txt"),"w").write("\n".join(sig))

# gene -> membership annotation
rows=[]
def add(gl,src):
    for g in gl: rows.append((g,src))
add(heme,"HALLMARK_HEME_METABOLISM")
for k,v in curated.items(): add(v,"curated_"+k)
for k,v in gobp.items(): add(sorted(v)[:0] or [], "")  # skip dumping huge GO names into annotation
ann=pd.DataFrame(rows,columns=["gene","source"]).drop_duplicates()
ann.to_csv(os.path.join(RES,"RBC_signature_membership.csv"),index=False)

# ---- intersection with SAH DEGs ----
tt36=pd.read_csv(os.path.join(RES,"DEG_GSE36791.csv"))
deg_strict=tt36[(tt36["adj.P.Val"]<0.05)&(tt36["logFC"].abs()>1)]
deg_broad =tt36[tt36["adj.P.Val"]<0.05]
print(f"\nDEG strict (FDR<.05,|lfc|>1)={len(deg_strict)} | DEG broad (FDR<.05)={len(deg_broad)}")
cand=deg_broad[deg_broad.gene.isin(sig)].copy()
cand=cand.merge(ann,on="gene",how="left")
cand=cand.sort_values("P.Value")
cand.to_csv(os.path.join(RES,"candidate_RBC_SAH_genes.csv"),index=False)
print(f"CANDIDATE genes (SAH FDR<.05 DEG  ∩  RBC/HRR signature) = {len(cand)}")
print(cand[["gene","logFC","adj.P.Val","source"]].to_string(index=False))

# also candidate from strict
cs=deg_strict[deg_strict.gene.isin(sig)]
print(f"candidates from strict DEGs = {len(cs)}: {sorted(cs.gene.tolist())}")

# ---- program activity score per sample: SAH vs Control, both cohorts ----
def program_score(tag, geneset):
    lg=np.load(os.path.join(DAT,f"lnorm_{tag}.npy"))
    g=open(os.path.join(DAT,f"genes_{tag}.txt")).read().split()
    gi={x:i for i,x in enumerate(g)}
    idx=[gi[x] for x in geneset if x in gi]
    z=(lg[idx]-lg[idx].mean(axis=1,keepdims=True))/ (lg[idx].std(axis=1,keepdims=True)+1e-9)
    sc=z.mean(axis=0)
    ph=pd.read_csv(os.path.join(DAT,f"pheno_aligned_{tag}.csv")).set_index("gsm").loc[
        pd.read_csv(os.path.join(DAT,f"expr_linear_{tag}.tsv"),sep="\t",index_col=0).columns.tolist()]
    grp=ph["group"].values
    a=sc[grp=="SAH"]; b=sc[grp=="Control"]
    t,p=stats.ttest_ind(a,b); u,pu=stats.mannwhitneyu(a,b)
    d=(a.mean()-b.mean())/np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
    print(f"\nprogram-score {tag}: SAH mean={a.mean():.3f} Ctrl mean={b.mean():.3f} "
          f"t={t:.2f} p={p:.3g} MWU p={pu:.3g} Cohen d={d:.2f} (n={len(idx)} sig genes)")
    return pd.DataFrame({"gsm":ph.index,"group":grp,"score":sc}), (t,p,d)

sc36,st36=program_score("GSE36791", sig)
sc73,st73=program_score("GSE73378", sig)
pd.concat([sc36.assign(cohort="GSE36791"),sc73.assign(cohort="GSE73378")]).to_csv(
    os.path.join(RES,"program_activity_scores.csv"),index=False)

# ---- ORA of candidate genes via Enrichr ----
gl=cand.gene.tolist()
if len(gl)>=3:
    for lib in ["GO_Biological_Process_2023","KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2022"]:
        try:
            enr=gseapy.enrichr(gene_list=gl, gene_sets=lib, organism="human", outdir=None)
            r=enr.results.sort_values("Adjusted P-value")
            r.to_csv(os.path.join(RES,f"ORA_candidates_{lib}.csv"),index=False)
            print(f"\n--- ORA candidates on {lib} (top) ---")
            print(r[["Term","P-value","Adjusted P-value","Overlap"]].head(8).to_string(index=False))
        except Exception as e:
            print(lib,"ORA ERR",repr(e)[:120])
print("\nDONE 05_signature")
