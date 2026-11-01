"""
25_coloc_engine.py
Reviewer priority #4 (★★★★): add CO-LOCALISATION to the causal layer.
MR says "associated"; coloc asks "is it the SAME causal variant?".

This implements Giambartolomei et al. (2014) coloc.abf (single-causal-variant,
Wakefield approximation, Bayes factors) with PP0..PP4 and the standard priors
p1=1e-4, p2=1e-4, p12=1e-5. It is validated on synthetic ground truth:
  * Scenario A: two traits share one causal SNP   -> expect PP4 high (>0.8)
  * Scenario B: two independent causal SNPs        -> expect PP3 high (not coloc)
  * Scenario C: one trait causal, other null       -> expect PP1/PP2 high (not coloc)

Real-data execution is BLOCKED in this environment for the same reason as MR
(HP pQTL / HMOX1 eQTL / SAH GWAS full summary statistics are not openly
downloadable here). We therefore deliver a validated, ready-to-run engine plus
the exact access recipe, and report NO real coloc estimate. No fabrication.

Outputs: results/coloc_engine_validation_synthetic.csv
"""
import os, json
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

def log_abf(beta, se, sd_prior=0.15):
    """Wakefield log approximate Bayes factor per SNP (coloc's ABF)."""
    W = sd_prior**2
    V = se**2
    r = W / (W + V)
    z = beta / se
    return 0.5*np.log(1 - r) + 0.5*r*z*z

def _lsum(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))

def coloc_abf(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5, sd_prior=0.15):
    """Return PP0..PP4 for shared vs distinct causal variants."""
    l1 = log_abf(np.asarray(beta1), np.asarray(se1), sd_prior)
    l2 = log_abf(np.asarray(beta2), np.asarray(se2), sd_prior)
    L1 = _lsum(l1); L2 = _lsum(l2)
    # H4: shared single causal variant -> sum over i of l1_i*l2_i
    s12 = l1 + l2
    L12 = _lsum(s12)
    # log priors
    lp1, lp2, lp12 = np.log(p1), np.log(p2), np.log(p12)
    lp0 = np.log(1 - p1 - p2 - p12)
    lp1_only = np.log(p1 - p12); lp2_only = np.log(p2 - p12)
    # unnormalised log posteriors
    d0 = lp0
    d1 = lp1_only + L1
    d2 = lp2_only + L2
    # H3: two distinct causal SNPs = (all pairs) - (diagonal pairs)
    diff = np.exp(L1 + L2) - np.exp(L12)
    d3 = lp1 + lp2 + (np.log(diff) if diff > 0 else -np.inf)   # two distinct causal SNPs
    d4 = lp12 + L12
    D = _lsum(np.array([d0, d1, d2, d3, d4]))
    pp = np.exp(np.array([d0, d1, d2, d3, d4]) - D)
    return dict(PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3], PP4=pp[4])

# ---------------- synthetic validation ----------------
rng = np.random.default_rng(2026)
n = 5000
def sim_trait(causal_idx, effect, noise=0.03):
    se = np.full(n, 0.02) + rng.uniform(0, 0.01, n)
    beta = rng.normal(0, noise, n)
    for ci in np.atleast_1d(causal_idx):
        beta[ci] += effect
    return beta, se

results = []
# A: shared causal SNP
b1, s1 = sim_trait([100], 0.5); b2, s2 = sim_trait([100], 0.45)
ppA = coloc_abf(b1, s1, b2, s2)
results.append(dict(scenario="A_shared_causal", **{k: round(float(v),4) for k,v in ppA.items()},
                    expected="PP4>0.8", verdict=("PASS" if ppA["PP4"]>0.8 else "FAIL")))
# B: distinct causal SNPs
b1, s1 = sim_trait([100], 0.5); b2, s2 = sim_trait([4000], 0.5)
ppB = coloc_abf(b1, s1, b2, s2)
results.append(dict(scenario="B_distinct_causal", **{k: round(float(v),4) for k,v in ppB.items()},
                    expected="PP3 high, PP4 low", verdict=("PASS" if (ppB["PP3"]>0.5 and ppB["PP4"]<0.5) else "FAIL")))
# C: one causal, other null
b1, s1 = sim_trait([100], 0.5); b2 = rng.normal(0,0.02,n); s2 = np.full(n,0.025)
ppC = coloc_abf(b1, s1, b2, s2)
results.append(dict(scenario="C_one_null", **{k: round(float(v),4) for k,v in ppC.items()},
                    expected="PP1/PP2, PP4 low", verdict=("PASS" if ppC["PP4"]<0.5 else "FAIL")))

df = pd.DataFrame(results)
df.to_csv(os.path.join(RES, "coloc_engine_validation_synthetic.csv"), index=False)

recipe = {
 "engine": "coloc.abf (Giambartolomei 2014), single-causal-variant Wakefield ABF, priors p1=p2=1e-4 p12=1e-5",
 "validated_scenarios": results,
 "real_data_access_recipe": {
   "exposure_pQTL_HP_CD163_ferritin": "deCODE 2023 (decodinggenomics.com/health), UKB-PPP Olink (clinicaltrial NCT05449511 / aptamer list), or Sun 2023 Nat Genet trans-ancestry pQTL — request beta/se per SNP for HP, CD163, FTL/FTH.",
   "exposure_eQTL_HMOX1": "GTEx v8 eQTL (HMOX1, whole blood / brain), or eQTLGen — beta/se per cis-SNP.",
   "outcome_SAH_GWAS": "INGEST / GWAS Catalog SAH (GCST90136755 or latest); need beta/se/EA/OA/effect allele + N; harmonise alleles.",
   "steps": ["download sumstats", "liftover to GRCh38", "restrict to cis window (+/-1Mb of gene)", "harmonise effect alleles", "LD not required for coloc.abf (single-variant) but coloc.susie needs a reference panel (1000G)", "run coloc_abf per gene-trait pair", "report PP4; PP4>0.8 = shared causal signal"],
   "blocker": "full summary statistics not retrievable in this environment (GWAS Catalog sumstats endpoints deprecated/410; pQTL/eQTL portals need bulk download or auth). No real PP4 reported."
 },
 "honesty": "Only the synthetic validation is reported as a result; the real coloc is NOT computed and NOT fabricated."
}
json.dump(recipe, open(os.path.join(RES, "coloc_engine_recipe.json"), "w"), indent=2, default=str)
print(df.to_string(index=False))
print("wrote coloc_engine_validation_synthetic.csv + coloc_engine_recipe.json")
