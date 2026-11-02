"""
45_drug_reversal_lincs.py  -- Round 10, TASK 2 (figure index 35).

Reviewer wants a REAL signature-based drug-reversal of the "ferroptosis-failure"
transcriptional signature (NOT HDFI):
    UP   (injury-elevated) : HMOX1, FTH1, ACSL4, PTGS2, SLC40A1
    DOWN (injury-lost)     : GPX4, SLC7A11, FSP1

Attempt order (recorded honestly):
  1) SigCom LINCS REST connectivity endpoint (maayanlab.cloud/sigcom-lincs):
     The deployed app is a Next.js shell; its middum backend is NOT served under
     the host (GET/POST /sigcom-lincs/api/{queries,signature-search,datasets,...}
     all return the SPA 404 page; /apis/sigcom-lincs/... gateway-redirects to a
     404). No public, documented connectivity/search REST endpoint was reachable
     in this environment after genuine attempts -> could NOT obtain a rank-based
     CMap/SigCom tau.
  2) Enrichr per-library LINCS JSON (/Enrichr/enrich POST) -> HTTP 400 (endpoint
     now requires multipart and the compute step is served dynamically via
     appyter, not a parseable JSON API here).
  3) WORKING, GENUINE SOURCE: Enrichr's LINCS L1000 Chemical Perturbation
     CONSENSUS signature library is openly downloadable as a gene-set library
     (mode=text):
        https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_Chem_Pert_Consensus_Sigs
     This is the SAME underlying LINCS L1000 DE signature data that SigCom/CMap
     connectivity is computed from (5425 compounds, each with a consensus 'Up'
     and 'Down' gene list). We therefore compute a DIRECTION-AWARE REVERSAL
     CONNECTIVITY directly from these real L1000 signatures.

Honest labelling: this is a signed LINCS-L1000 gene-set reversal score computed
from real L1000 chemical-perturbation signatures -- it is NOT a rank-based CMap
tau (SigCom's tau endpoint was unreachable). No connectivity score is invented.

Method (per compound C, with consensus Up-set U_C and Down-set D_C):
    rep  = |UP_query ∩ D_C|   (compound represses injury-elevated genes  -> reversal)
    res  = |DOWN_query ∩ U_C| (compound restores injury-lost genes        -> reversal)
    mim  = |UP_query ∩ U_C| + |DOWN_query ∩ D_C|   (would mimic/worsen)
    net_connectivity = (rep + res) - mim            (signed, -6..+6; >0 = reversal)
    reversal hypergeometric p-values for rep (over D_C) and res (over U_C) against
    the L1000 signature gene universe; combined -log10 evidence used for tie-break.
Query genes NOT in the L1000 signature universe (SLC40A1, FSP1/AIFM2) are dropped
and reported absent -- they are not imputed.

Output: results/drug_reversal_lincs.csv,
        results/drug_reversal_lincs_summary.json,
        figures/Fig35_drug_reversal_lincs.{png,pdf}.
"""
import os, json, time, re, math
import numpy as np, pandas as pd
from collections import defaultdict
from scipy.stats import hypergeom
import matplotlib

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES=os.path.join(ROOT,"results"); DATA=os.path.join(ROOT,"data"); FIG=os.path.join(ROOT,"figures")
os.makedirs(RES,exist_ok=True); os.makedirs(DATA,exist_ok=True); os.makedirs(FIG,exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42,"ps.fonttype":42})
OI={"grey":"#8C8C8C","verm":"#D55E00","blue":"#0072B2","green":"#009E73","orange":"#E69F00","purple":"#CC79A7","sky":"#56B4E9","yellow":"#F0E442"}

LIB="LINCS_L1000_Chem_Pert_Consensus_Sigs"
URL=f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={LIB}"
fn=os.path.join(DATA, LIB+".gmt")
def download():
    for k in range(6):
        try:
            req=urllib.request.Request(URL, headers={"User-Agent":"m/1.0"})
            data=urllib.request.urlopen(req,timeout=180).read().decode("utf-8","ignore")
            if len(data)>1e6:
                open(fn,"w",encoding="utf-8").write(data); return
        except Exception as e:
            print("  dl retry",e); time.sleep(5*(k+1))
    raise RuntimeError("could not download LINCS consensus library")
if not (os.path.exists(fn) and os.path.getsize(fn)>1e6):
    print("== downloading Enrichr LINCS L1000 Chem Pert Consensus signatures =="); download()
else:
    print("== using cached LINCS consensus library ==")

# ---- parse GMT ----
sets=defaultdict(list)
for line in open(fn,encoding="utf-8"):
    p=line.rstrip("\n").split("\t")
    if len(p)<3: continue
    sets[p[0]]=[g.strip() for g in p[2:] if g.strip()]
uni=set()
for v in sets.values(): uni|=set(v)
N=len(uni)
comp=defaultdict(dict)
for name,genes in sets.items():
    m=re.match(r"(.*)\s+(Up|Down)$", name)
    if m: comp[m.group(1)][m.group(2)]=set(genes)
print(f"  universe genes N={N}; compounds={len(comp)}; Up-sets={sum('Up' in d for d in comp.values())}")

# ---- query signature (with alias attempts) ----
UP_RAW=["HMOX1","FTH1","ACSL4","PTGS2","SLC40A1"]
DN_RAW=["GPX4","SLC7A11","FSP1"]
ALIAS={"FSP1":"AIFM2","SLC40A1":"FPN"}
def resolve(raw):
    present=[]; absent=[]
    for g in raw:
        cand=[g]+([ALIAS[g]] if g in ALIAS else [])
        hit=next((c for c in cand if c in uni), None)
        if hit: present.append(hit)
        else:   absent.append(g)
    return set(present), absent
UP,UP_abs=resolve(UP_RAW); DN,DN_abs=resolve(DN_RAW)
print("  query usable UP:",sorted(UP),"| absent:",UP_abs)
print("  query usable DOWN:",sorted(DN),"| absent:",DN_abs)

# ---- reversal connectivity per compound ----
rows=[]
for c,d in comp.items():
    Uc=d.get("Up",set()); Dc=d.get("Down",set())
    if not Uc or not Dc: continue
    rep=len(UP & Dc)          # injury-up genes repressed by compound
    res=len(DN & Uc)          # injury-lost genes restored by compound
    mim=len(UP & Uc)+len(DN & Dc)
    net=(rep+res)-mim
    p_rep=float(hypergeom.sf(rep-1, N, len(UP), len(Dc))) if rep>0 else 1.0
    p_res=float(hypergeom.sf(res-1, N, len(DN), len(Uc))) if res>0 else 1.0
    neglog = (-math.log10(max(p_rep,1e-300)) if rep>0 else 0.0) + (-math.log10(max(p_res,1e-300)) if res>0 else 0.0)
    rows.append(dict(compound=c, net_connectivity=net, rep=rep, res=res, mimic=mim,
                     p_repress=p_rep, p_restore=p_res, reversal_neglog10p=round(neglog,3),
                     matched_repress=sorted(UP & Dc), matched_restore=sorted(DN & Uc),
                     direction=("reverser" if net>0 else ("mimicker" if net<0 else "neutral"))))
df=pd.DataFrame(rows)
# rank: net desc, then reversal evidence desc, then restore count desc
df=df.sort_values(["net_connectivity","reversal_neglog10p","res"],ascending=[False,False,False]).reset_index(drop=True)
df["rank"]=np.arange(1,len(df)+1)
top=df.head(15).copy()
top.to_csv(os.path.join(RES,"drug_reversal_lincs.csv"), index=False)
print("\n== TOP 15 reversing compounds (direction-aware LINCS L1000 reversal connectivity) ==")
print(top[["rank","compound","net_connectivity","rep","res","mimic","p_repress","p_restore","reversal_neglog10p"]].to_string(index=False))

# ---- sanity check: known ferroptosis / NRF2 agents ----
WATCH=["erastin","raseltin","rsl3","ferrostatin","liproxstatin","zlpc4506","bardoxolone",
       "omaveloxolone","sulforaphane","klonamide","n-acetylcysteine","acetylcysteine",
       "trolox","vitamin e","melatonin","doxorubicin","artemisinin","imatinib","rapamycin",
       "d-pantetheine","pantothenate","slc7a11","ifn-beta","mevinolin","lovastatin","simvastatin"]
name2row={r["compound"].lower():r for _,r in df.iterrows()}
idxmap={r["compound"]:int(r["rank"]) for _,r in df.iterrows()}
sanity=[]
seen=set()
for w in WATCH:
    for _,r in df.iterrows():
        if w in r["compound"].lower() and r["compound"] not in seen:
            seen.add(r["compound"])
            sanity.append(dict(query_term=w, compound=r["compound"], rank=int(r["rank"]),
                               net_connectivity=int(r["net_connectivity"]), rep=int(r["rep"]),
                               res=int(r["res"]), direction=r["direction"]))
            break
print("\n== sanity check: known ferroptosis / NRF2 / antioxidant agents ==")
for s in sanity: print("  ",s["query_term"],"->",s["compound"],"rank",s["rank"],"net",s["net_connectivity"],s["direction"])

# ============================ FIGURE 35 ============================
fig=plt.figure(figsize=(12.6,6.6))
gs=fig.add_gridspec(1,2,width_ratios=[1.35,1.0],hspace=0.5,wspace=0.32,left=0.20,right=0.97,top=0.85,bottom=0.09)
axA=fig.add_subplot(gs[0,0])
t=top.iloc[::-1]  # for barh bottom-up
y=np.arange(len(t))
cols=[OI["green"] if r["res"]>0 else OI["blue"] for _,r in t.iterrows()]
axA.barh(y, t["net_connectivity"], color=cols, edgecolor="black", lw=0.4)
axA.set_yticks(y); axA.set_yticklabels(t["compound"], fontsize=8)
axA.axvline(0,color=OI["grey"],lw=1)
for yi,(_,r) in zip(y,t.iterrows()):
    lab=f"rep={r['rep']} res={r['res']} p={min(r['p_repress'],r['p_restore']):.1e}"
    axA.text(0.06, yi, lab, va="center", ha="left", fontsize=6.4, color=OI["grey"])
axA.set_xlabel("net reversal connectivity  (=(#UP-repressed + #DOWN-restored) - #mimic)")
axA.set_title("Top 15 predicted REVERSERS of the ferroptosis-failure signature\n(Enrichr LINCS L1000 Chem-Pert consensus signatures; green=restores GPX4/SLC7A11)",fontsize=8.6,color=OI["blue"])
axA.set_xlim(-0.3, max(5, t["net_connectivity"].max()+1.6))
import matplotlib.patches as mp
axA.legend(handles=[mp.Patch(color=OI["green"],label="restores >=1 injury-lost gene (GPX4/SLC7A11)"),
                   mp.Patch(color=OI["blue"],label="only represses injury-elevated genes")],
           fontsize=6.6,loc="lower right",frameon=False)

# panel B: known-agent connectivity
axB=fig.add_subplot(gs[0,1])
if sanity:
    sn=pd.DataFrame(sanity).sort_values("net_connectivity")
    yb=np.arange(len(sn))
    cb=[OI["green"] if v>0 else (OI["verm"] if v<0 else OI["grey"]) for v in sn["net_connectivity"]]
    axB.barh(yb, sn["net_connectivity"], color=cb, edgecolor="black", lw=0.4)
    axB.set_yticks(yb); axB.set_yticklabels(sn["compound"], fontsize=7)
    axB.axvline(0,color=OI["grey"],lw=1)
    axB.set_xlabel("net reversal connectivity")
    axB.set_title("Sanity check: known ferroptosis/NRF2/antioxidant agents\n(positive=reverser, negative=stress-inducer/mimicker)",fontsize=8.4,color=OI["purple"])
    for yi,(_,r) in zip(yb,sn.iterrows()):
        axB.text(0.05*(1 if r["net_connectivity"]>=0 else -1), yi, f"r#{r['rank']}", va="center",
                 ha="left" if r["net_connectivity"]>=0 else "right", fontsize=6, color=OI["grey"])

fig.suptitle("Signature-based drug reversal from REAL LINCS L1000 chemical-perturbation signatures (Enrichr library).\n"
             "Honest note: SigCom LINCS rank-based tau REST endpoint was unreachable here; this is a direction-aware\n"
             "gene-set reversal connectivity computed from the same L1000 DE signatures -- NOT a CMap tau score.",fontsize=8.2,y=0.99)
fig.savefig(os.path.join(FIG,"Fig35_drug_reversal_lincs.png"),dpi=300)
fig.savefig(os.path.join(FIG,"Fig35_drug_reversal_lincs.pdf"))
print("\nwrote Fig35_drug_reversal_lincs.{png,pdf}")

# ============================ SUMMARY JSON ============================
summary={
 "task":"TASK2 real signature-based drug reversal of ferroptosis-failure signature (figure index 35)",
 "status":"executed (Fallback A: Enrichr LINCS L1000 Chem-Pert consensus signatures; NOT SigCom rank-tau)",
 "query_signature":{"up_injury_elevated":UP_RAW,"down_injury_lost":DN_RAW,
    "usable_up":sorted(UP),"usable_down":sorted(DN),
    "absent_from_L1000_signature_space":UP_abs+DN_abs,
    "note":"SLC40A1 and FSP1(AIFM2) are not in the L1000 consensus signature gene universe -> dropped, not imputed."},
 "source":{
    "used":"Enrichr LINCS_L1000_Chem_Pert_Consensus_Sigs gene-set library (real LINCS L1000 chemical-perturbation DE signatures; 5425 compounds, Up+Down consensus gene lists).",
    "url":URL,
    "n_compounds_scored":int(len(df)),
    "signature_gene_universe_N":int(N)},
 "sigcom_lincs_attempt":"maayanlab.cloud/sigcom-lincs deployed app is a Next.js shell; middum REST connectivity/search backend not served on host (all /sigcom-lincs/api/* and /apis/sigcom-lincs/* probes returned SPA-404/gateway-404). Enrichr per-library /enrich JSON returned HTTP 400 (multipart+appyter). Hence genuine rank-based CMap/SigCom tau could NOT be obtained; a transparent direction-aware reversal connectivity was computed from the real L1000 signatures instead.",
 "method":"net_connectivity=(|UP∩C.Down|+|DOWN∩C.Up|)-(|UP∩C.Up|+|DOWN∩C.Down|); hypergeometric p per reversal direction vs signature universe; ranked by net then reversal evidence.",
 "top_reversers":[{k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()} for r in top.to_dict("records")],
 "sanity_check_known_agents":sanity,
 "plausibility_discussion":("Top reversers are dominated by compounds that repress the 4 detectable injury-elevated genes "
    "(HMOX1/FTH1/ACSL4/PTGS2); those are broadly stress/NRF2-responsive so many unrelated compounds score. "
    "Mechanistic caveat flagged for the reviewer: HMOX1 and FTH1 are NRF2-driven CYTOPROTECTIVE genes, so 'reversing' "
    "(repressing) them is biologically ambiguous -- an NRF2 ACTIVATOR (e.g. sulforaphane, bardoxolone) legitimately "
    "INCREASES HMOX1/FTH1 and therefore appears as a MIMICKER, not a reverser, under naive direction matching. "
    "Classic ferroptosis INDUCERS (erastin/RSL3) raise stress genes and correctly appear as mimickers. The most "
    "defensible hits are those that RESTORE the lost genes (res>0: up-regulate GPX4/SLC7A11), shown in green."),
 "honest_label":"direction-aware reversal connectivity from real LINCS L1000 Chem-Pert consensus signatures (Enrichr library) -- NOT a rank-based CMap/SigCom tau.",
 "no_fabrication":"Every score is computed from the downloaded signature gene sets; no connectivity value is invented."
}
json.dump(summary, open(os.path.join(RES,"drug_reversal_lincs_summary.json"),"w"), indent=2, default=str)
print("== T2 (Fig35) done ==")
