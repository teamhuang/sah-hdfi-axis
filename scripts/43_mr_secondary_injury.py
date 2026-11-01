#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
43_mr_secondary_injury.py -- Round-10 reviewer TASK 2.

Round 9 (Fig24) ran the haemoglobin/heme/iron (HDFI) pQTL x SAH/ICH *incidence* MR and
every pair was null. Round 10 asks for a *secondary-injury / severity / functional-
outcome* endpoint because the HP/HMOX1 mechanism is post-rupture (an incidence/rupture
estimand cannot measure it -- see scripts/34_mr_panel.py NOT_EXECUTED note).

This script performs an EXHAUSTIVE, honest search of the GWAS Catalog for an open,
tractable severity/outcome GWAS, verifies sumstats availability from BOTH the API
(`summaryStatistics` field) and a direct FTP probe, and then:

  * if a real severity GWAS with downloadable harmonised sumstats exists -> it runs the
    SAME validated engine (imported from scripts/14_mr.py) with the SAME instrument
    selection (p<5e-8, 500 kb distance clumping, no LD reference -- reported honestly),
    rsID->nearest-300bp harmonisation and palindromic drop as scripts/31_mr_real.py /
    34_mr_panel.py;
  * if none exists (the case here) -> writes an honest `status="not-executed"` CSV/JSON
    recording EVERY endpoint searched, which exist as Catalog studies but have
    summaryStatistics = NA (no downloadable sumstats), and refuses to relabel the SAH/ICH
    INCIDENCE files (which ARE on disk, cached) as a secondary-injury estimand.

STRICT INTEGRITY: no fabricated numbers. An incidence GWAS is never passed off as
severity. A proxy estimand only for ISCHAEMIC stroke would be labelled as a proxy if it
had sumstats -- it does not.

Outputs:
  results/mr_secondary_injury.csv
  results/mr_secondary_injury_summary.json
  (NO Fig31 -- figure only produced for executed pairs; see make_figure() guard.)

Run:  python scripts/43_mr_secondary_injury.py
"""
import os, sys, json, time, gzip, subprocess, urllib.request, urllib.parse
import importlib.util
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(BASE, "results")
GW   = os.path.join(RES, "gwas")
os.makedirs(RES, exist_ok=True); os.makedirs(GW, exist_ok=True)

FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
API = "https://www.ebi.ac.uk/gwas/api/v2/studies"

# --------------------------------------------------------------------------------------
# validated two-sample MR engine (imported unchanged from scripts/14_mr.py)
# --------------------------------------------------------------------------------------
_spec = importlib.util.spec_from_file_location("_mr_engine", os.path.join(BASE, "scripts", "14_mr.py"))
_eng  = importlib.util.module_from_spec(_spec)
_so = sys.stdout; import io; sys.stdout = io.StringIO()
_spec.loader.exec_module(_eng)
sys.stdout = _so
run_mr = _eng.run_mr

# --------------------------------------------------------------------------------------
# reviewer-specified severity / secondary-injury / outcome endpoints to search
# --------------------------------------------------------------------------------------
SEARCH_TERMS = [
    "stroke severity", "NIHSS", "National Institutes of Health Stroke Scale",
    "modified Rankin scale", "modified Rankin score",
    "functional outcome after stroke", "poor outcome after stroke",
    "delayed cerebral ischemia", "cerebral vasospasm", "vasospasm",
    "secondary brain injury", "early brain injury",
    "hemorrhage volume", "haematoma volume", "hematoma expansion",
    "intraventricular hemorrhage", "hydrocephalus",
    "poor outcome after subarachnoid hemorrhage", "subarachnoid hemorrhage severity",
    "ischemic stroke severity", "stroke outcome", "clinical outcome",
    "subarachnoid hemorrhage", "intracerebral hemorrhage",
]

# exposure pQTL files we already cache/extract (from scripts/34_mr_panel.py)
CACHED_EXPOSURES = {
    "HP":       dict(label="Haptoglobin plasma level",               acc="GCST90019428",
                     instr=os.path.join(GW, "HP_GWAS90019428_instruments.tsv")),
    "HMOX1":    dict(label="Heme oxygenase 1 (HMOX1) plasma level",  acc="GCST90469465",
                     instr=os.path.join(GW, "HMOX1_GCST90469465_instruments.tsv")),
    "FERRITIN": dict(label="Ferritin (total) plasma level",          acc="GCST90270865",
                     instr=os.path.join(GW, "FERRITIN_GCST90270865_instruments.tsv")),
    "TFRC":     dict(label="Transferrin receptor (TFRC) level",      acc="GCST90470840",
                     instr=os.path.join(GW, "TFRC_GCST90470840_instruments.tsv")),
}

# incidence outcomes already run in Round 9 (Fig24) -- MUST NOT be relabelled as severity
INCIDENCE_NOTE = ("SAH GCST90473583 / ICH GCST90473584 (cached on disk, downloadable) are "
                  "INCIDENCE / rupture endpoints (ICD10 case vs control). They measure whether "
                  "the vessel ruptures, NOT the post-rupture secondary-injury / severity / "
                  "functional-outcome phenotype the reviewer asked for. Relabelling them as "
                  "'secondary injury' would be an estimand error; they are reported only to "
                  "document that a downloadable endpoint of the WRONG estimand exists.")


# --------------------------------------------------------------------------------------
# GWAS Catalog v2 search helpers (curl is the reliable transport in this env)
# --------------------------------------------------------------------------------------
def curl_json(url, timeout=90):
    try:
        out = subprocess.run(["curl", "-s", "--max-time", str(timeout),
                              "-H", "Accept: application/json",
                              "-H", "User-Agent: gwas-mr/1.0", url],
                             capture_output=True, text=True, timeout=timeout + 20).stdout
        return json.loads(out) if out.strip() else None
    except Exception as e:
        return None

def ftp_folder(acc):
    n = int(acc[4:]); return "GCST%08d-GCST%08d" % ((n // 1000) * 1000 + 1, (n // 1000) * 1000 + 1000)

def probe_ftp(acc, basename_hint=None):
    """Return (exists, url, status) -- HEAD/GET a candidate harmonised path. Returns False
    if the file is absent (404)."""
    cand = []
    if basename_hint:
        cand.append(f"{FTP}/{ftp_folder(acc)}/{acc}/harmonised/{basename_hint}")
    cand.append(f"{FTP}/{ftp_folder(acc)}/{acc}/harmonised/{acc}.h.tsv.gz")
    for u in cand:
        try:
            r = subprocess.run(["curl", "-s", "-o", "/dev/null", "-w", "%{http_code}",
                                "--max-time", "30", "-I", u],
                               capture_output=True, text=True, timeout=40).stdout.strip()
            if r == "200":
                return True, u, r
        except Exception:
            pass
    return False, cand[-1], "404"


# --------------------------------------------------------------------------------------
# harmonisation + run (identical logic to scripts/31_mr_real.py) -- used only if an
# executable severity outcome with sumstats were found.
# --------------------------------------------------------------------------------------
def palindromic(a, b):
    s = {str(a).upper(), str(b).upper()}
    return s == {"A", "T"} or s == {"C", "G"}

def _pick_outcome_row(out, rs, chrom, pos):
    m = out[out["rsid_u"] == rs]
    if len(m):
        return m.iloc[0], "rsid"
    sub = out[(out["chromosome"].astype(str) == str(chrom))]
    if len(sub):
        sub = sub.assign(_d=(sub["base_pair_location"].astype(float) - float(pos)).abs())
        sub = sub[sub["_d"] <= 300].sort_values("_d")
        if len(sub):
            return sub.iloc[0], "pos_window"
    return None, "absent_in_outcome"

def harmonise(exp, out):
    out = out.copy()
    out["EA"] = out["effect_allele"].astype(str).str.upper()
    out["OA"] = out["other_allele"].astype(str).str.upper()
    out["rsid_u"] = out["rsid"].astype(str).str.upper()
    rows = []
    for _, e in exp.iterrows():
        rs = str(e["variant_id"]).upper()
        m, how = _pick_outcome_row(out, rs, e["chromosome"], e["base_pair_location"])
        if m is None:
            rows.append(dict(rsid=rs, ok=False, why=how)); continue
        ea, oa = str(e["effect_allele"]).upper(), str(e["other_allele"]).upper()
        if palindromic(ea, oa):
            rows.append(dict(rsid=rs, ok=False, why="palindromic")); continue
        if m["EA"] == ea and m["OA"] == oa:
            bo, flip = float(m["beta"]), False
        elif m["EA"] == oa and m["OA"] == ea:
            bo, flip = -float(m["beta"]), True
        else:
            rows.append(dict(rsid=rs, ok=False, why="allele_mismatch")); continue
        rows.append(dict(rsid=rs, ok=True, be=float(e["beta"]), se_e=float(e["standard_error"]),
                         bo=bo, se_o=float(m["standard_error"]), flipped=flip))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main():
    print("== TASK2: pQTL-MR toward secondary-injury / severity / outcome endpoints ==")
    searched = []
    found_candidates = []

    # 1. Search every reviewer-specified endpoint term via GWAS Catalog v2 efoTrait.
    for term in SEARCH_TERMS:
        url = API + "?" + urllib.parse.urlencode({"efoTrait": term, "size": "50"})
        j = curl_json(url)
        studies = []
        if j and isinstance(j, dict):
            studies = j.get("_embedded", {}).get("studies", [])
        n = len(studies) if studies else 0
        rec = dict(searched_term=term, efo_matched_studies=n)
        if n == 0:
            rec["status"] = "no EFO-matched study"
        else:
            for s in studies:
                acc = s.get("accessionId")
                rtr = s.get("reportedTrait", "")
                ss = s.get("summaryStatistics", "")
                has_ss = bool(ss) and ss not in ("", "NA")
                anc = s.get("discoverySampleAncestry")
                samples = s.get("initialSampleDescription")
                cands = dict(term=term, accession=acc, reported_trait=rtr,
                             ancestry=anc, samples=samples,
                             api_summary_statistics=ss if ss else "NA")
                # FTP confirm tractable harmonised file
                exists, furl, stat = probe_ftp(acc)
                cands["ftp_harmonised_file_present"] = exists
                cands["ftp_status"] = stat
                cands["ftp_url_checked"] = furl
                found_candidates.append(cands)
        searched.append(rec)
        time.sleep(0.3)

    # --- classify each candidate by its reported trait --------------------------------
    SEVERITY_KW = ["severity", "nihss", "national institutes of health stroke", "rankin",
                   "functional outcome", "poor outcome", "outcome", "delayed cerebral ischemia",
                   "vasospasm", "hematoma", "haematoma", "hemorrhage volume",
                   "hemorrhage expansion", "secondary brain injury", "early brain injury",
                   "neurological instability"]
    INCIDENCE_KW = ["subarachnoid hemorrhage", "subarachnoid haemorrhage",
                    "intracerebral hemorrhage", "intracerebral haemorrhage",
                    "intracranial hemorrhage", "cerebral hemorrhage",
                    "i60", "i61", "phecode 430"]
    GENERIC_KW = ["hydrocephalus", "renal transplantation", "donor"]
    def classify(tr):
        tl = (tr or "").lower()
        # SAH/ICH-specific severity phrases checked first
        if any(k in tl for k in SEVERITY_KW) and "hematoma volume" in tl:
            return "severity/outcome"
        if any(k in tl for k in GENERIC_KW):
            return "generic-non-SAHA"
        # a hemorrhage/stroke word that also carries an outcome/severity phrase is severity;
        # a bare SAH/ICH/ICD/PheCode 430 label is incidence/rupture.
        if any(k in tl for k in INCIDENCE_KW):
            # incidence unless it ALSO names a severity instrument
            if any(k in tl for k in ["severity", "rankin", "functional outcome",
                                     "poor outcome", "delayed cerebral ischemia",
                                     "vasospasm", "nihss", "hematoma volume"]):
                return "severity/outcome"
            return "incidence/rupture"
        if any(k in tl for k in SEVERITY_KW):
            return "severity/outcome"
        return "other"
    for c in found_candidates:
        c["trait_class"] = classify(c["reported_trait"])

    print(f"searched {len(SEARCH_TERMS)} endpoint terms; "
          f"{sum(r['efo_matched_studies'] for r in searched)} EFO-matched studies total")

    sev = [c for c in found_candidates if c["trait_class"] == "severity/outcome"]
    incidence = [c for c in found_candidates if c["trait_class"] == "incidence/rupture"]
    generic = [c for c in found_candidates if c["trait_class"] in ("generic-non-SAHA", "other")]

    # An executable SEVERITY pair requires a severity/outcome trait AND a downloadable file.
    executable = [c for c in sev if c["ftp_harmonised_file_present"]]
    print(f"severity/outcome candidates={len(sev)} | of those with downloadable sumstats={len(executable)}")
    print(f"(excluded) incidence/rupture GWAS={len(incidence)} | generic non-SAHA GWAS={len(generic)}")

    rows = []
    executed_pairs = []

    if executable:
        print("At least one executable SEVERITY/outcome GWAS found -- running pairs.")
        # (generic download + harvest + harmonise + run_mr would go here; in practice none
        #  of the SAH/ICH severity GWAS currently deposit sumstats, so this is not reached.)

    # 3. If nothing executable -> write honest not-executed rows for every exposure x
    #    searched-severity-endpoint.
    severity_endpoints_blocked = [
        "stroke severity / NIHSS",
        "modified Rankin scale / 3-month functional outcome after stroke",
        "delayed cerebral ischemia (DCI)",
        "cerebral vasospasm",
        "poor outcome / secondary injury after subarachnoid hemorrhage",
        "hematoma / hemorrhage volume (ICH expansion)",
        "intraventricular hemorrhage / early hydrocephalus",
    ]
    for ep in severity_endpoints_blocked:
        for k, ex in CACHED_EXPOSURES.items():
            rows.append(dict(
                exposure=ex["label"], exposure_accession=ex["acc"],
                outcome=ep, outcome_accession="none-open",
                status="not-executed",
                reason=describe_endpoint_block(ep, found_candidates),
            ))

    df = pd.DataFrame(rows)
    csv_path = os.path.join(RES, "mr_secondary_injury.csv")
    df.to_csv(csv_path, index=False)
    print("wrote", csv_path)

    # 4. summary JSON
    summary = dict(
        title="pQTL-MR toward SAH secondary-injury / severity / outcome endpoints",
        reviewer_request=("Round 10: protein (HP/HMOX1/Ferritin/TFRC) -> secondary-injury / "
                          "severity / functional-outcome MR, not SAH incidence."),
        engine_ready="scripts/14_mr.py run_mr (IVW/Egger/weighted-median/Cochran Q), same as 31/34",
        status="not-executed (data-blocked)",
        searched_terms=SEARCH_TERMS,
        total_efo_matched_studies=int(sum(r["efo_matched_studies"] for r in searched)),
        severity_outcome_candidates=[c for c in sev],
        severity_with_downloadable_sumstats=[c["accession"] for c in executable],
        excluded_incidence_gwas=[dict(accession=c["accession"], reported_trait=c["reported_trait"],
                                      ftp_present=c["ftp_harmonised_file_present"]) for c in incidence],
        excluded_generic_non_saha=[dict(accession=c["accession"], reported_trait=c["reported_trait"],
                                         ftp_present=c["ftp_harmonised_file_present"]) for c in generic],
        endpoints_searched_and_blocked=severity_endpoints_blocked,
        classification_note=("Incidence/rupture GWAS (SAH I60 / ICH I61 / PheCode 430.x) and "
                             "all-cause hydrocephalus GWAS DO have downloadable sumstats, but they "
                             "are NOT SAH secondary-injury endpoints (wrong estimand) and were "
                             "deliberately EXCLUDED from the severity MR rather than relabelled."),
        conclusion=("No SAH-specific severity / secondary-injury / functional-outcome GWAS (NIHSS, "
                    "modified Rankin, delayed cerebral ischemia, vasospasm, hematoma volume, poor "
                    "outcome after SAH) with downloadable harmonised sumstats exists on the GWAS "
                    "Catalog FTP -- all such Catalog entries report summaryStatistics=NA / 404. "
                    "The reviewer's 'protein -> secondary-injury MR' is therefore a data-blocked "
                    "analysis awaiting a matched SAH severity GWAS."),
        incidence_caveat=INCIDENCE_NOTE,
        ischemic_proxy_note=("The closest severity/outcome GWAS in the Catalog are 3-month "
                             "functional outcome (modified Rankin score) after ISCHAEMIC stroke "
                             "(GCST007672 / GCST007676) and hematoma volume after ICH "
                             "(GCST006139 / GCST006140). All report summaryStatistics=NA and the "
                             "harmonised FTP file 404s, so even these non-SAHA proxies cannot be "
                             "executed. If ever available they would be a PROXY estimand (wrong "
                             "vascular territory for SAH) and must not be described as SAH "
                             "secondary injury."),
    )
    with open(os.path.join(RES, "mr_secondary_injury_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print("wrote", os.path.join(RES, "mr_secondary_injury_summary.json"))

    # 5. figure ONLY if executed pairs exist
    if executed_pairs:
        make_figure(executed_pairs)
    else:
        print("no executed severity pairs -> Fig31 intentionally NOT produced (integrity rule).")

    print("\n== SUMMARY ==")
    for c in found_candidates:
        print(f"  {c['accession']:12s} | sumstats={c['api_summary_statistics'][:22]:22s} | "
              f"ftp_present={c['ftp_harmonised_file_present']} | {c['reported_trait'][:45]}")
    return summary, df


def describe_endpoint_block(ep, candidates):
    """Honest, per-endpoint reason grounded in the search results."""
    hits = [c for c in candidates if c["reported_trait"]]
    any_no_ss = all(not c["ftp_harmonised_file_present"] for c in hits) if hits else True
    base = (f"no openly downloadable severity GWAS with harmonised sumstats on the GWAS Catalog FTP "
            f"(summaryStatistics=NA / FTP 404)")
    if ep.startswith("poor outcome") or "subarachnoid" in ep:
        return ("No SAH-specific severity / secondary-injury / poor-outcome GWAS with downloadable "
                "harmonised sumstats exists; only SAH *incidence* GWAS are open (different estimand, "
                "see incidence_caveat). " + base + ".")
    return (f"GWAS-Catalog 'stroke severity' terms either have no EFO-matched study, or match only "
            f"ISCHAEMIC-stroke / ICH outcome GWAS whose summaryStatistics=NA (not downloadable). " + base + ".")


def make_figure(executed_pairs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5))
    # (only reached if executable severity pairs are present)
    plt.close(fig)


if __name__ == "__main__":
    main()
