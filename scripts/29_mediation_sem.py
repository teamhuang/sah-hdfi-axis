"""
29_mediation_sem.py
Reviewer round-7 #4: build the clinical<->mechanistic BRIDGE as a mediation /
structural-equation model, NOT a hard gene->HRR link:

    HRR (X) -> heme-detoxification failure (M1) -> brain injury (M2) -> death (Y)

Implemented: (a) simple mediation (Baron-Kenny paths a,b,c,c' + ACIE=a*b +
Sobel + bootstrap CI); (b) serial mediation X->M1->M2->Y with bootstrap
indirect effects. Both engines are validated on synthetic data where the true
indirect effect is known.

Real-data execution is BLOCKED: no public cohort pairs admission Hb/RDW with
the molecular mediators, and MIMIC-IV (which has HRR + outcome) is not present
in this project folder. We therefore report ONLY the validated engine + the
exact analysis recipe, and DO NOT fabricate a real mediation estimate.
Outputs: results/mediation_engine_validation_synthetic.csv, mediation_engine_recipe.json
"""
import os, json
import numpy as np, pandas as pd
import statsmodels.api as sm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
rng = np.random.default_rng(2026)

def _ols(y, X):
    return sm.OLS(y, sm.add_constant(X)).fit()

def simple_mediation(X, M, Y, nboot=2000, seed=2026):
    """ACIE = a*b (X->M->Y); ADE = c'; total = c. Sobel + bootstrap CI."""
    r = _ols(Y, X); c = r.params[1]
    ra = _ols(M, X); a = ra.params[1]
    rb = _ols(Y, np.column_stack([X, M])); b = rb.params[2]; cp = rb.params[1]
    acie = a*b
    se_a, se_b = ra.bse[1], rb.bse[2]
    sobel_z = acie/np.sqrt(b*b*se_a*se_a + a*a*se_b*se_b)
    from scipy.stats import norm
    sobel_p = 2*(1-norm.cdf(abs(sobel_z)))
    g = np.random.default_rng(seed); boots=[]
    n=len(X)
    for _ in range(nboot):
        i=g.integers(0,n,n)
        Xs,Ms,Ys=X[i],M[i],Y[i]
        aa=_ols(Ms,Xs).params[1]; bb=_ols(Ys,np.column_stack([Xs,Ms])).params[2]
        boots.append(aa*bb)
    lo,hi=np.percentile(boots,[2.5,97.5])
    return dict(total_c=float(c), a=float(a), b=float(b), cprime=float(cp),
                ACIE=float(acie), sobel_z=float(sobel_z), sobel_p=float(sobel_p),
                boot_lo=float(lo), boot_hi=float(hi), sig=bool(lo>0 or hi<0))

def serial_mediation(X, M1, M2, Y, nboot=2000, seed=2026):
    """X->M1->M2->Y. Indirect paths: X->M1->M2->Y (a1*d21*b2), X->M1->Y, X->M2->Y."""
    a1=_ols(M1,X).params[1]                       # X->M1
    d21=_ols(M2,np.column_stack([X,M1])).params[2] # M1->M2 (adj X)
    ym=_ols(Y,np.column_stack([X,M1,M2]))
    b2=ym.params[3]; b1=ym.params[2]; cp=ym.params[1]
    ind_full=a1*d21*b2; ind_m1=a1*b1; ind_m2=_ols(M2,X).params[1]*b2 - d21*ind_m1*0  # X->M2 direct-ish
    # cleaner X->M2->Y: regress M2 on X then *b2
    a2=_ols(M2,X).params[1]; ind_m2=a2*b2
    total=_ols(Y,X).params[1]
    g=np.random.default_rng(seed); n=len(X); boots=[]
    for _ in range(nboot):
        i=g.integers(0,n,n)
        Xs,M1s,M2s,Ys=X[i],M1[i],M2[i],Y[i]
        aa1=_ols(M1s,Xs).params[1]; dd=_ols(M2s,np.column_stack([Xs,M1s])).params[2]
        yy=_ols(Ys,np.column_stack([Xs,M1s,M2s])); bb2=yy.params[3]
        boots.append(aa1*dd*bb2)
    lo,hi=np.percentile(boots,[2.5,97.5])
    return dict(path_XM1M2Y=float(ind_full), path_XM1Y=float(ind_m1), path_XM2Y=float(ind_m2),
                direct_cprime=float(cp), total=float(total), serial_ACIE=float(ind_full),
                serial_boot_lo=float(lo), serial_boot_hi=float(hi), sig=bool(lo>0 or hi<0))

# ---------------- synthetic validation ----------------
n=800
X = rng.normal(0,1,n)
# simple mediation, true a=0.8,b=0.6 -> ACIE=0.48 ; direct 0.1
M = 0.8*X + rng.normal(0,1,n)
Y = 0.6*M + 0.1*X + rng.normal(0,1,n)
sm1 = simple_mediation(X,M,Y)
# serial, true X->M1(0.7)->M2(0.7)->Y(0.7) => serial indirect ~0.343
M1 = 0.7*X + rng.normal(0,1,n)
M2 = 0.7*M1 + rng.normal(0,1,n)
Y2 = 0.7*M2 + 0.1*M1 + 0.1*X + rng.normal(0,1,n)
sr = serial_mediation(X,M1,M2,Y2)
# negative control: no mediation (X not related to M)
Mnull = rng.normal(0,1,n); Ynull = 0.5*X + rng.normal(0,1,n)
sm_null = simple_mediation(X,Mnull,Ynull)

rows=[dict(engine="simple", true_ACIE=0.48, est_ACIE=round(sm1["ACIE"],3), ci=f"[{sm1['boot_lo']:.3f},{sm1['boot_hi']:.3f}]", sig=sm1["sig"],
            verdict="PASS" if (0.4<sm1["ACIE"]<0.6 and sm1["sig"]) else "FAIL"),
      dict(engine="serial", true_ACIE=0.343, est_ACIE=round(sr["serial_ACIE"],3), ci=f"[{sr['serial_boot_lo']:.3f},{sr['serial_boot_hi']:.3f}]", sig=sr["sig"],
            verdict="PASS" if (0.25<sr["serial_ACIE"]<0.45 and sr["sig"]) else "FAIL"),
      dict(engine="negative_control", true_ACIE=0.0, est_ACIE=round(sm_null["ACIE"],3), ci=f"[{sm_null['boot_lo']:.3f},{sm_null['boot_hi']:.3f}]", sig=sm_null["sig"],
            verdict="PASS" if (not sm_null["sig"]) else "FAIL")]
df=pd.DataFrame(rows); df.to_csv(os.path.join(RES,"mediation_engine_validation_synthetic.csv"), index=False)

recipe=dict(
 model="HRR -> heme-detoxification-failure (HDFI / plasma HP, free-Hb, HMOX1) -> brain injury (infarct volume / NIHSS / ICP) -> death/mRS",
 engines=["Baron-Kenny + bootstrap ACIE (simple)", "serial mediation X->M1->M2->Y (bootstrap)", "SEM with lavaan/semopy (full latent path, robust SE)"],
 binary_outcome_note="For death/mRS use natural-effects (AIPW) mediation or logistic SEM; linear ACIE here is for engine validation only.",
 real_data_inputs_needed=["MIMIC-IV SAH: admission Hb, RDW -> HRR; in-hospital death, mRS",
                          "matched plasma at admission: HP, free Hb, HMOX1, ferritin, MDA/4-HNE (proxy for HDFI) OR an HDFI from paired blood RNA",
                          "imaging severity (Fisher/modified Fisher, DCI, infarct) as M2"],
 steps=["assemble cohort with HRR + mediator + severity + outcome", "standardise/one-hot covariates (age, sex, WFNS)",
        "fit mediator~HRR+cov; outcome~mediator+HRR+cov; bootstrap ACIE", "sensitivity: E-value / natural-effects for unmeasured mediator-outcome confounding",
        "report % mediated = ACIE/total"],
 blocker="no public cohort pairs HRR with the molecular mediators; MIMIC-IV not present in this folder. No real mediation estimate is reported — engine only.",
 validation=rows, honesty="Only synthetic validation is a result; real mediation is NOT computed and NOT fabricated.")
json.dump(recipe, open(os.path.join(RES,"mediation_engine_recipe.json"),"w"), indent=2, default=str)
print(df.to_string(index=False))
print("wrote mediation_engine_validation_synthetic.csv + mediation_engine_recipe.json")
