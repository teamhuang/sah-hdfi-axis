"""13_ml_xgb.py -- gradient-boosting model comparison across feature sets:
A=DEGs, B=RBC/HRR candidates, C=heme-iron panel, D=single-cell-derived macrophage/heme markers.
Repeated-CV AUC (discovery) + external AUC (GSE73378) for each."""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
try:
    from xgboost import XGBClassifier; HAVE_XGB=True
except Exception:
    HAVE_XGB=False

BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES,DAT=os.path.join(BASE,"results"),os.path.join(BASE,"data")

def cohort(tag):
    lg=np.load(os.path.join(DAT,f"lnorm_{tag}.npy")); g=np.array(open(os.path.join(DAT,f"genes_{tag}.txt")).read().split())
    ph=pd.read_csv(os.path.join(DAT,f"expr_linear_{tag}.tsv"),sep="\t",index_col=0).columns.tolist()
    pd_=pd.read_csv(os.path.join(DAT,f"pheno_aligned_{tag}.csv")).set_index("gsm").loc[ph]
    return lg,{x:i for i,x in enumerate(g)},(pd_["group"].values=="SAH").astype(int)
lg36,i36,y36=cohort("GSE36791"); lg73,i73,y73=cohort("GSE73378")
def z(lg,idx): M=lg[idx].T; return (M-M.mean(0))/(M.std(0)+1e-9)

deg=pd.read_csv(os.path.join(RES,"DEG_GSE36791.csv"))
A=list(deg[(deg["adj.P.Val"]<0.05)&(deg["logFC"].abs()>1)].gene)
B=list(pd.read_csv(os.path.join(RES,"candidate_RBC_SAH_genes.csv")).gene)
C=[g for g in ["HP","HMOX1","FTH1","FTL","SLC40A1","SLC25A37","NFE2","ALAS2","FECH","BLVRB",
               "TFRC","GPX4","ACSL4","MARC1","AHSP","GYPB","BCL11A"] ]
# D: single-cell-derived macrophage/erythrophagocytosis-heme markers (human orthologs)
Dmeas=[g for g in ["HP","CD163","MRC1","HMOX1","FTH1","SLC40A1","FTMT","APOE","TREM2","LPL","CST7","C1QA"] if g in i36]
if len(Dmeas)<4: Dmeas=[g for g in ["HP","HMOX1","FTH1","SLC40A1","CST7","APOE"] if g in i36]
D=Dmeas
sets={"A_DEGs":A,"B_RBC_candidates":B,"C_heme_iron_panel":C,"D_singlecell_markers":D}

def model():
    if HAVE_XGB:
        return XGBClassifier(n_estimators=200,max_depth=3,learning_rate=0.05,subsample=0.8,
                             colsample_bytree=0.8,eval_metric="auc",random_state=2026,verbosity=0)
    return HistGradientBoostingClassifier(max_iter=200,max_depth=3,learning_rate=0.05,random_state=2026)

rows=[]
for name,gl in sets.items():
    gl=[g for g in gl if g in i36]
    X=z(lg36,[i36[g] for g in gl])
    cv=RepeatedStratifiedKFold(n_splits=5,n_repeats=20,random_state=2026)
    aucs=cross_val_score(model(),X,y36,cv=cv,scoring="roc_auc")
    # external
    gl73=[g for g in gl if g in i73]
    ext=np.nan
    if len(gl73)>=2:
        Xtr=z(lg36,[i36[g] for g in gl73]); Xte=z(lg73,[i73[g] for g in gl73])
        m=model().fit(Xtr,y36); ext=roc_auc_score(y73,m.predict_proba(Xte)[:,1])
    rows.append({"feature_set":name,"n_genes":len(gl),"CV_AUC_mean":aucs.mean(),
                 "CV_AUC_sd":aucs.std(),"external_AUC_GSE73378":ext})
    print(f"{name}: n={len(gl)} CV_AUC={aucs.mean():.3f}±{aucs.std():.3f} external={ext if np.isnan(ext) else round(ext,3)}")
S=pd.DataFrame(rows); S.to_csv(os.path.join(RES,"xgb_model_comparison.csv"),index=False)
print("\nbest internal:",S.sort_values('CV_AUC_mean',ascending=False).iloc[0].feature_set)
print("DONE 13_ml_xgb")
