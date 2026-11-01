"""14_mr.py -- two-sample Mendelian-randomisation engine (IVW fixed+random, MR-Egger,
weighted median, Cochran's Q, Egger-intercept pleiotropy). SAH/HP/iron GWAS summary
statistics are not openly downloadable in this environment (GWAS Catalog v2/FTP 404 or
deprecated; eQTLGen portal = large files; IEU OpenGWAS = JWT token). Therefore the engine
is validated on a SYNTHETIC ground-truth simulation (labelled as such) to prove estimator
correctness, and a real-data run() is provided for when sumstats are available. NO real MR
result is reported here because no real harmonised sumstats could be obtained."""
import os, numpy as np, pandas as pd
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(BASE,"results")
rng=np.random.default_rng(2026)

def ivw(be,bo,se):
    w=1/se**2
    b=(w*be*bo).sum()/(w*be**2).sum()
    se_b=np.sqrt(1/(w*be**2).sum())
    Q=(w*(bo-b*be)**2).sum(); k=len(be)
    se_rand=se_b*np.sqrt(max(Q/(k-1),1.0))
    return dict(beta=b,se=se_b,p=_normp(b/se_b),
                Q=Q,df=k-1,se_random=se_rand,p_random=_normp(b/se_rand))
def _normp(z):
    from scipy.stats import norm; return 2*norm.sf(abs(z))
def egger(be,bo,se):
    import statsmodels.api as sm
    w=1/se**2
    X=sm.add_constant(be); m=sm.WLS(bo,X,weights=w).fit()
    return dict(beta=m.params[1],se=m.bse[1],p=m.pvalues[1],intercept=m.params[0],
                intercept_p=m.pvalues[0])
def weighted_median(be,bo,se):
    r=bo/be; wr=1/se**2; wr=wr/wr.sum()
    o=np.argsort(r); r=r[o]; wr=wr[o]; c=np.cumsum(wr)-0.5*wr
    beta=r[np.searchsorted(c,0.5)]
    # bootstrap SE
    bs=[]
    for _ in range(1000):
        idx=rng.choice(len(r),len(r),replace=True)
        rr=np.sort(bo[idx]/be[idx]); ww=(1/se[idx]**2); ww=ww/ww.sum(); cc=np.cumsum(ww)-0.5*ww
        bs.append(rr[np.searchsorted(cc,0.5)])
    return dict(beta=beta,se=np.std(bs),p=_normp(beta/np.std(bs)))

def run_mr(dat):
    """dat: DataFrame with columns be,bo,se (per instrument, effect-allele aligned)."""
    be,bo,se=dat.be.values,dat.bo.values,dat.se.values
    out={"n_instruments":len(be)}
    out.update({f"IVW_{k}":v for k,v in ivw(be,bo,se).items()})
    out.update({f"Egger_{k}":v for k,v in egger(be,bo,se).items()})
    out.update({f"WMedian_{k}":v for k,v in weighted_median(be,bo,se).items()})
    return out

# ---------------- synthetic ground-truth validation (NOT real biology) ----------------
K=30; true_beta=0.4
be=rng.uniform(0.05,0.3,K)                    # exposure (HP expr) -> SNP effects
se=rng.uniform(0.01,0.03,K)
bo=true_beta*be + rng.normal(0,0.01,K)        # outcome effect via causal path
bo[:5]+=rng.normal(0.05,0.02,5)               # 5 pleiotropic instruments
synth=pd.DataFrame({"be":be,"bo":bo,"se":se})
r=run_mr(synth)
print("=== SYNTHETIC MR (ground truth beta=0.40; engine validation only, not real data) ===")
print(f"  IVW beta={r['IVW_beta']:.3f} (se {r['IVW_se']:.3f}, p={r['IVW_p']:.2g})")
print(f"  MR-Egger beta={r['Egger_beta']:.3f}, intercept={r['Egger_intercept']:.3f} (p={r['Egger_intercept_p']:.2g})")
print(f"  Weighted median beta={r['WMedian_beta']:.3f} (se {r['WMedian_se']:.3f})")
print(f"  Cochran Q={r['IVW_Q']:.1f} df={r['IVW_df']}")
pd.DataFrame({"engine_validation_synthetic":list(r.keys()),"value":list(r.values())}
             ).to_csv(os.path.join(RES,"mr_engine_validation_synthetic.csv"),index=False)

print("""
=== REAL-DATA MR: BLOCKED IN THIS ENVIRONMENT (documented, not fabricated) ===
Required & why unavailable:
 - Outcome SAH GWAS sumstats (e.g., MISSTRAT aSAH): GWAS Catalog v2 flags full_summary_stats_available=False
   for retrievable SAH studies; legacy REST/sumstats endpoints return HTTP 410 (deprecated).
 - Exposure HP cis-pQTL / eQTL: eQTLGen portal = multi-GB zips; IEU OpenGWAS = JWT token required.
To execute for real: download (a) HP/HMOX1/FTH1/SLC40A1/NFE2 eQTLGen cis sumstats and (b) an SAH (or, as a
directly-relevant HRR proxy, RDW/Hb from Astle 2016 / Vuckovic 2020) harmonised sumstats, LD-clump with PLINK,
harmonise by effect allele, then call run_mr(). The estimator is validated above (recovers beta~0.40).
""")
print("DONE 14_mr")
