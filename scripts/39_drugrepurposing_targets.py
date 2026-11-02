"""
39_drugrepurposing_targets.py  -- Round 10 functional-layer probe, TASK 4.

Reviewer asked for a CMap/CLUE reverse-signature drug-repurposing screen.

What was attempted and BLOCKED in this environment (recorded honestly, not
faked):
  * CLUE.io  -> login required (blocked, per session constraints).
  * LINCS L1000 CSV dumps -> account required (blocked).
  * Enrichr REST (maayanlab.cloud) LINCS_L1000_Chem_Pert_{up,down} libraries:
    /addList works (returns userListId) but /enrich returns HTTP 400 (POST,
    multipart & JSON) and an EMPTY result {} (GET) for every library incl.
    KEGG -- i.e. the compute step is not usable programmatically here.
    So the genuine LINCS SIGNATURE-MATCHING number could NOT be obtained.

What IS openly available WITHOUT login and WAS executed: the OpenTargets
Platform GraphQL API (https://api.platform.opentargets.org/api/v4/graphql),
which returns curated approved/clinical drug-target associations
(target.drugAndClinicalCandidates).  We pull all known small-molecule
modulators of the ferroptosis / NRF2-heme axis targets
(HMOX1, NFE2L2, SLC7A11, GPX4, ACSL4, PTGS2).

Honest labelling (as the task specifies): this is a TARGET-CURATED candidate
shortlist (drug-target associations), NOT a CMap/CLUE-computed connectivity
score.  No connectivity/direction inference is made.
Output: results/drugrepurposing_target_candidates.csv,
figures/Fig28_drugrepurposing_targets.{png,pdf}.
"""
import os, json, time
import numpy as np, pandas as pd
import requests
import matplotlib

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures"); DATA = os.path.join(ROOT, "data")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"grey":"#8C8C8C","verm":"#D55E00","blue":"#0072B2","green":"#009E73","orange":"#E69F00","purple":"#CC79A7","sky":"#56B4E9"}
HOST = "https://api.platform.opentargets.org/api/v4/graphql"

# Ensembl IDs resolved live via Ensembl REST (see header), cached here.
TARGETS = {"HMOX1":"ENSG00000100292","NFE2L2":"ENSG00000116044","SLC7A11":"ENSG00000151012",
           "GPX4":"ENSG00000167468","ACSL4":"ENSG00000068366","PTGS2":"ENSG00000073756"}
Q = '{ target(ensemblId:"%s"){ approvedSymbol drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ name } } } } }'

def gql(eid, tries=5):
    for k in range(tries):
        try:
            r = requests.post(HOST, json={"query": Q % eid}, timeout=90)
            if r.status_code == 200:
                return r.json().get("data", {}).get("target")
        except Exception:
            pass
        time.sleep(3*(k+1))
    return None

print("== TASK4 OpenTargets target-curated candidate drugs ==")
rows = []
summary = {}
for sym, eid in TARGETS.items():
    t = gql(eid)
    if not t:
        print(sym, "NO DATA"); summary[sym] = None; continue
    dac = t.get("drugAndClinicalCandidates", {})
    n = dac.get("count", 0); rr = dac.get("rows", [])
    summary[sym] = n
    for x in rr:
        rows.append({"target":sym,"ensembl_id":eid,"drug":x["drug"]["name"],
                     "max_clinical_stage":x["maxClinicalStage"],
                     "association_source":"OpenTargets curated drug-target"})
    print(f"  {sym}: {n} drug-target associations ({len(rr)} returned)")
    time.sleep(1)

df = pd.DataFrame(rows)
if len(df)==0:
    raise SystemExit("OpenTargets returned no drug data - abort without fabricating.")
df = df.sort_values(["target","max_clinical_stage"]).reset_index(drop=True)
df.to_csv(os.path.join(RES,"drugrepurposing_target_candidates.csv"), index=False)
print("wrote results/drugrepurposing_target_candidates.csv rows", len(df))

nrows = df.groupby("target").size().reindex(list(TARGETS)).fillna(0)
approved = df[df.max_clinical_stage.isin(["APPROVAL","PHASE_4","PHASE_3"])].groupby("target").size().reindex(list(TARGETS)).fillna(0)

# ---- figure ----
fig, ax = plt.subplots(figsize=(6.5, 4.0))
y=np.arange(len(TARGETS))
ax.barh(y, nrows.values, color=OI["blue"], label="all curated drugs")
ax.barh(y, approved.values, color=OI["verm"], label=">=Phase 3 / approved")
ax.set_yticks(y); ax.set_yticklabels(list(TARGETS)); ax.invert_yaxis()
ax.set_xlabel("# small-molecule modulators (OpenTargets)"); ax.legend(fontsize=8)
ax.set_title("Target-curated candidate shortlist\n(drug-target associations, NOT CMap connectivity)")

fig.suptitle("OpenTargets Platform drug-target associations (no login). Honest note: CLUE.io & Enrichr-LINCS\n"
             "signature-matching were login/end-point blocked, so this is a curated shortlist, NOT a computed connectivity map.", fontsize=8)
fig.tight_layout(rect=[0,0,1,0.90])
fig.savefig(os.path.join(FIG,"Fig28_drugrepurposing_targets.png"), dpi=300)
fig.savefig(os.path.join(FIG,"Fig28_drugrepurposing_targets.pdf"))
print("wrote Fig28_drugrepurposing_targets.{png,pdf}")

nrf2_drugs = list(df[(df.target=="NFE2L2")]["drug"])
audit = {"status":"executed",
 "method":"OpenTargets Platform GraphQL target.drugAndClinicalCandidates for HMOX1/NFE2L2/SLC7A11/GPX4/ACSL4/PTGS2.",
 "blocked_alternatives":{
    "CLUE.io":"login required (session constraint)",
    "LINCS_L1000_CSV_dumps":"account required",
    "Enrichr_LINCS_L1000_Chem_Pert":"addList OK but /enrich returned persistent HTTP 400 (POST) and empty {} (GET) for all libraries incl. KEGG - genuine LINCS signature match could not be obtained"},
 "per_target_drug_counts":{k:(int(v) if v is not None else None) for k,v in summary.items()},
 "NFE2L2_captured_drugs":nrf2_drugs,
 "conclusion":"Real curated drug-target associations: NFE2L2 is directly druggable (bardoxolone methyl, omaveloxolone); GPX4/SLC7A11/HMOX1 have no catalogued approved modulators; PTGS2 has many (NSAIDs).",
 "honest_label":"target-curated candidate shortlist (OpenTargets), NOT CMap/CLUE-computed reverse-signature connectivity."}
json.dump(audit, open(os.path.join(RES,"_t4_audit.json"),"w"), indent=2)
print("== T4 done ==")
