"""06_ml.py -- ML hub-gene selection on the mechanism-informed DEG feature space:
LASSO (L1 logistic) + Random-Forest importance + SVM-RFE; consensus hub genes."""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.inspection import permutation_importance

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT=os.path.join(BASE,"results"),os.path.join(BASE,"data")
rng=np.random.RandomState(2026)

def feats(tag):
    lg=np.load(os.path.join(DAT,f"lnorm_{tag}.npy")); g=np.array(open(os.path.join(DAT,f"genes_{tag}.txt")).read().split())
    ph=pd.read_csv(os.path.join(DAT,f"expr_linear_{tag}.tsv"),sep="\t",index_col=0).columns.tolist()
    pd_=pd.read_csv(os.path.join(DAT,f"pheno_aligned_{tag}.csv")).set_index("gsm").loc[ph]
    y=(pd_["group"].values=="SAH").astype(int)
    return lg,g,y

# feature space = strict DEGs (FDR<.05,|lfc|>1) UNION RBC/SAH candidate genes
deg=pd.read_csv(os.path.join(RES,"DEG_GSE36791.csv"))
strict=set(deg[(deg["adj.P.Val"]<0.05)&(deg["logFC"].abs()>1)].gene)
cand=set(pd.read_csv(os.path.join(RES,"candidate_RBC_SAH_genes.csv")).gene)
feat=sorted(strict|cand)
print(f"feature space: strict={len(strict)} candidate={len(cand)} union={len(feat)}")

lg,g,y=feats("GSE36791"); gi={x:i for i,x in enumerate(g)}
rows=[gi[f] for f in feat if f in gi]; X=lg[rows].T; feat=[g[i] for i in rows]
X=StandardScaler().fit_transform(X)
print(f"X={X.shape} SAH={y.sum()} Control={(1-y).sum()}")
cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=2026)

# ---- 1. LASSO (L1 logistic) ----
lr=LogisticRegressionCV(Cs=np.logspace(-3,1,30),penalty="l1",solver="liblinear",cv=cv,
                        max_iter=5000,random_state=2026,scoring="roc_auc").fit(X,y)
coefs=lr.coef_[0]
lasso_sel=np.array(feat)[coefs!=0]
lasso_rank=pd.Series(np.abs(coefs),index=feat).sort_values(ascending=False)
print(f"LASSO best C={lr.C_[0]:.4g}; nonzero={len(lasso_sel)}; top={list(lasso_rank.head(12).index)}")

# ---- 2. Random Forest importance ----
rf=RandomForestClassifier(n_estimators=2000,class_weight="balanced",random_state=2026,n_jobs=-1).fit(X,y)
gini=pd.Series(rf.feature_importances_,index=feat).sort_values(ascending=False)
perm=permutation_importance(rf,X,y,n_repeats=50,random_state=2026,scoring="roc_auc").importances_mean
perm=pd.Series(perm,index=feat).sort_values(ascending=False)
print("RF gini top:",list(gini.head(12).index))
print("RF permutation top:",list(perm.head(12).index))

# ---- 3. SVM-RFE ----
svc=SVC(kernel="linear",C=1.0,class_weight="balanced",random_state=2026)
rfecv=RFECV(svc,step=1,cv=cv,scoring="roc_auc",min_features_to_select=1).fit(X,y)
svm_rank=pd.Series(rfecv.ranking_,index=feat).sort_values()   # 1=most important
svm_top=list(svm_rank[svm_rank<=10].index)
print(f"SVM-RFE optimal #feat={int(rfecv.support_.sum())}; top10={svm_top}")

# ---- consensus hub genes ----
def topn(s,n): return set(s.head(n).index)
K=12
L=topn(lasso_rank,K); Rf=topn(perm,K); Gf=topn(gini,K); S=set(svm_top[:K])
score=pd.DataFrame({"lasso":lasso_rank.rank(),"rf_perm":perm.rank(),"rf_gini":gini.rank(),"svm_rfe":1/svm_rank.rank(ascending=False)})
# average percentile rank (lower=better) across the four rankings
allg=sorted(feat)
def pct(series):
    r=series.reindex(allg)
    return r.rank(ascending=False, method="average", pct=True)
combined=pd.DataFrame({"lasso":pct(lasso_rank),"rf_perm":pct(perm),"rf_gini":pct(gini),"svm_rfe":pct(1/svm_rank)}).dropna(how="all")
combined["mean_pct"]=combined.mean(axis=1)
combined=combined.sort_values("mean_pct")
combined.to_csv(os.path.join(RES,"ml_ranking_all.csv"))
# HUB = genes in top-K of >=3 of the 4 rankings (robust consensus)
votes=pd.Series({gene: (gene in L)+(gene in Rf)+(gene in Gf)+(gene in S) for gene in allg}).sort_values(ascending=False)
hub=sorted(votes[votes>=3].index)
print(f"\nCONSENSUS HUB genes (top-{K} in >=3 of {['LASSO','RF-perm','RF-gini','SVM-RFE']}): {hub}")
print("vote table:\n",votes.head(15).to_string())
relaxed=list(combined.head(8).index)
print("Top-8 by combined percentile rank:",relaxed)

pd.DataFrame({"hub_gene":hub}).to_csv(os.path.join(RES,"hub_genes.csv"),index=False)
pd.DataFrame({"relaxed_top8":relaxed}).to_csv(os.path.join(RES,"hub_genes_relaxed.csv"),index=False)

# annotated hub table
src=pd.read_csv(os.path.join(RES,"RBC_signature_membership.csv")).drop_duplicates("gene")
ann=deg[deg.gene.isin(hub)][["gene","logFC","adj.P.Val"]].merge(
        pd.DataFrame({"gene":hub,"votes":[int(votes[g]) for g in hub]}),on="gene",how="left")
ann=ann.merge(src[["gene","source"]],on="gene",how="left")
ann["in_candidate_RBC"]=ann.gene.isin(cand)
ann=ann.sort_values(["votes","logFC"],ascending=[False,False])
ann.to_csv(os.path.join(RES,"hub_genes_annotated.csv"),index=False)
print("Hub annotated:\n",ann.to_string(index=False))
votes.to_csv(os.path.join(RES,"ml_vote_table.csv"),header=["votes"])

# mechanism-informed overlap: how many hub/relaxed genes are RBC candidates?
print("hub ∩ candidate-RBC:", sorted(set(hub)&cand), "| relaxed ∩ candidate:", sorted(set(relaxed)&cand))

# save coefficient table for downstream risk score
coef=pd.DataFrame({"gene":feat,"lasso_coef":coefs,
                   "rf_gini":rf.feature_importances_,
                   "rf_perm":perm.reindex(feat).values,
                   "svm_rank":svm_rank.reindex(feat).values})
coef.to_csv(os.path.join(RES,"ml_feature_table.csv"),index=False)
print("\nDONE 06_ml")
