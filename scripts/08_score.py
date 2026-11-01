"""08_score.py -- build a parsimonious SAH molecular score (SMS), evaluate with repeated CV in
discovery and independent external validation in GSE73378 (past/remote aSAH) via within-cohort
gene-wise standardisation + permutation; report honest AUCs."""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_score, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT=os.path.join(BASE,"results"),os.path.join(BASE,"data")

def cohort(tag):
    lg=np.load(os.path.join(DAT,f"lnorm_{tag}.npy")); g=np.array(open(os.path.join(DAT,f"genes_{tag}.txt")).read().split())
    ph=pd.read_csv(os.path.join(DAT,f"expr_linear_{tag}.tsv"),sep="\t",index_col=0).columns.tolist()
    pd_=pd.read_csv(os.path.join(DAT,f"pheno_aligned_{tag}.csv")).set_index("gsm").loc[ph]
    return lg,g,(pd_["group"].values=="SAH").astype(int)

lg36,g36,y36=cohort("GSE36791"); lg73,g73,y73=cohort("GSE73378")
i36={x:i for i,x in enumerate(g36)}; i73={x:i for i,x in enumerate(g73)}

def zmat(lg,idx):
    M=lg[idx].T
    return (M-M.mean(0))/(M.std(0)+1e-9)

def cv_auc(X,y,C=0.1):
    lr=LogisticRegression(penalty="l1",solver="liblinear",C=C,max_iter=5000,random_state=2026)
    rskf=RepeatedStratifiedKFold(n_splits=5,n_repeats=20,random_state=2026)
    aucs=cross_val_score(lr,X,y,cv=rskf,scoring="roc_auc")
    return float(np.mean(aucs))

models={}
# full 66-feature union model
deg=pd.read_csv(os.path.join(RES,"DEG_GSE36791.csv"))
strict=set(deg[(deg["adj.P.Val"]<0.05)&(deg["logFC"].abs()>1)].gene)
cand=set(pd.read_csv(os.path.join(RES,"candidate_RBC_SAH_genes.csv")).gene)
feat=sorted(strict|cand)
f36=[i36[f] for f in feat if f in i36]; ffeat=[g36[i] for i in f36]
f73=[i73[f] for f in ffeat if f in i73]; ffeat73=[g73[i] for i in f73]
common=[f for f in ffeat if f in set(ffeat73)]
models["M_full_l1"]=(ffeat,y36,zmat(lg36,f36))

# parsimonious scores
hub=pd.read_csv(os.path.join(RES,"hub_genes.csv")).hub_gene.tolist()
mech=["HP","SLC40A1","FTH1","SLC25A37","GYPB","NFE2","BCL11A","MARC1"]
hub_idx=[i36[h] for h in hub if h in i36]
mech_idx=[i36[m] for m in mech if m in i36]
models["M_hub4"]=( [g36[i] for i in hub_idx], y36, zmat(lg36,hub_idx) )
models["M_mech"]=( [g36[i] for i in mech_idx], y36, zmat(lg36,mech_idx) )

summary=[]
for name,(feat36,y,X) in models.items():
    auc_cv=cv_auc(X,y)
    lr=LogisticRegression(penalty="l1",solver="liblinear",C=0.1,max_iter=5000,random_state=2026).fit(X,y)
    auc_app=roc_auc_score(y,lr.predict_proba(X)[:,1])
    # external validation on genes shared by both cohorts (order preserved)
    pos=[k for k,f in enumerate(feat36) if f in i73]
    featC=[feat36[k] for k in pos]
    Xc=X[:,pos]
    Xv=zmat(lg73,[i73[f] for f in featC])           # same gene order as Xc
    lrX=LogisticRegression(penalty="l1",solver="liblinear",C=0.1,max_iter=5000,random_state=2026).fit(Xc,y)
    pv=lrX.predict_proba(Xv)[:,1]
    auc_ext=roc_auc_score(y73,pv)
    rs=np.random.RandomState(2026)
    nulls=[roc_auc_score(y73, lrX.predict_proba(Xv[rs.permutation(len(y73))])[:,1]) for _ in range(50)]
    perm_p=(np.sum([n>=auc_ext for n in nulls])+1)/(len(nulls)+1)
    summary.append({"model":name,"n_feat":len(feat36),"AUC_discovery_apparent":auc_app,
                    "AUC_discovery_CV20x5":auc_cv,"AUC_external_GSE73378":auc_ext,
                    "external_perm_p":perm_p,"null_AUC_mean":np.mean(nulls)})
    print(f"{name}: feats={len(feat36)} apparent={auc_app:.3f} CV={auc_cv:.3f} "
          f"EXTERNAL={auc_ext:.3f} (perm p={perm_p:.3f}, null mean={np.mean(nulls):.3f})")

S=pd.DataFrame(summary); S.to_csv(os.path.join(RES,"score_validation_summary.csv"),index=False)

# save ROC curve for the parsimonious hub model (discovery CV predictions) for plotting
lr=LogisticRegression(penalty="l1",solver="liblinear",C=0.1,max_iter=5000,random_state=2026)
X=models["M_hub4"][2]
pr=cross_val_predict(lr,X,y36,cv=StratifiedKFold(n_splits=5,shuffle=True,random_state=2026),method="predict_proba")[:,1]
fpr,tpr,_=roc_curve(y36,pr); pd.DataFrame({"fpr":fpr,"tpr":tpr}).to_csv(os.path.join(RES,"roc_hub_discovery.csv"),index=False)
# also external scores per-sample for the hub model (genes present in both cohorts)
h=[x for x in hub if x in i36 and x in i73]
Xd=zmat(lg36,[i36[x] for x in h]); lrf=LogisticRegression(penalty="l1",solver="liblinear",C=0.1,max_iter=5000,random_state=2026).fit(Xd,y36)
coefs=dict(zip(h,lrf.coef_[0]))
hv=[i73[x] for x in h]; Xv=zmat(lg73,hv); sc_ext=Xv@list(coefs.values())
pd.DataFrame({"gsm":pd.read_csv(os.path.join(DAT,"expr_linear_GSE73378.tsv"),sep="\t",index_col=0).columns,
              "group":y73,"SMS":sc_ext}).to_csv(os.path.join(RES,"SMS_GSE73378.csv"),index=False)
pd.Series(coefs,name="coef").to_csv(os.path.join(RES,"SMS_hub_coefficients.csv"))
print("\nHub SMS coefficients:",coefs)
print("DONE 08_score")
