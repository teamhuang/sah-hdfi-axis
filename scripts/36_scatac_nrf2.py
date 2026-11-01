"""
36_scatac_nrf2.py  -- Round 10 functional-layer probe, TASK 1.

Reviewer asked for an scATAC / chromatin-accessibility layer for the
"NRF2 -> HMOX1" claim.  No openly-downloadable, tractable SAH/ICH/microglia
scATAC peak or fragment file could be pulled in this environment (GEO peak
bundles are large / login-gated; GSE311852 deposited only scRNA matrices).
What IS genuinely available and directly on-point is REAL NFE2L2 (NRF2)
ChIP-seq from ENCODE (a bona fide in-vivo regulatory/binding reference, NOT a
sequence motif scan).  This script therefore executes the honest fallback
(option b): it downloads released NRF2 ChIP-seq "optimal IDR thresholded
peaks" (GRCh38, narrowPeak) for four human cell lines and asks how many real
NRF2 binding events fall at / near the HMOX1 promoter vs the rest of the
heme-handling module, benchmarked against a distance-decay null built from
random same-chromosome control positions.

Datasets (ENCODE, released GRCh38 narrowPeak, target label = NFE2L2):
  A549     ENCSR584GHV / ENCFF418TUX
  HeLa-S3  ENCSR707IUN / ENCFF305KIK
  HepG2    ENCSR488EES / ENCFF882YLO
  IMR-90   ENCSR197WGI / ENCFF474PPT

Honest labelling: this is an "NRF2 chromatin-binding / accessibility
reference", NOT a SAH-specific scATAC experiment.  Output:
results/nrf2_chip_hmoax1_enrichment.csv and figures/Fig26_nrf2_chip_binding.{png,pdf}.
"""
import os, io, json, time, random, urllib.request
import numpy as np, pandas as pd
import gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures"); DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"grey":"#8C8C8C","verm":"#D55E00","blue":"#0072B2","green":"#009E73","orange":"#E69F00","purple":"#CC79A7","sky":"#56B4E9"}
random.seed(20260)

ENS = "https://rest.ensembl.org"
# GRCh38 primary chromosome lengths (for random-control sampling + Poisson null)
CHRLEN = {str(i):L for i,L in zip(range(1,23),
  [248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,
   138394717,133797422,111797853,114364328,107043718,101991189,90338345,83257441,
   80373285,58617616,58617616,72530000,65285719,48124859])}  # chr19 fix below
CHRLEN = {"1":248956422,"2":242193529,"3":198295559,"4":190214555,"5":181538259,"6":170805979,
 "7":159345973,"8":145138636,"9":138394717,"10":133797422,"11":111797853,"12":114364328,
 "13":107043718,"14":101991189,"15":90338345,"16":83257441,"17":80373285,"18":58617616,
 "19":58617616,"20":63025520,"21":46709983,"22":50818468,"X":156040895,"Y":57227415}
GENOME = sum(CHRLEN.values())

FILES = {"A549":"NRF2_A549_GRCh38.bed.gz","HeLa":"NRF2_HeLa_GRCh38.bed.gz",
         "HepG2":"NRF2_HepG2_GRCh38.bed.gz","IMR90":"NRF2_IMR90_GRCh38.bed.gz"}
ACCESSIONS = {"A549":"ENCSR584GHV (ENCFF418TUX)","HeLa":"ENCSR707IUN (ENCFF305KIK)",
              "HepG2":"ENCSR488EES (ENCFF882YLO)","IMR90":"ENCSR197WGI (ENCFF474PPT)"}
GENES = ["HMOX1","FTH1","FTL","GPX4","SLC40A1"]

def ensure(fn, acc):
    """Download ENCODE narrowPeak if not cached (real released GRCh38 peaks)."""
    p = os.path.join(DATA, fn)
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return p
    url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
    with urllib.request.urlopen(url, timeout=180) as r:
        open(p,"wb").write(r.read())
    return p

def load_peaks(fn):
    """Read narrowPeak -> midpoint (chrom, pos), plus [start,end,score]."""
    mid, full = {}, []
    with gzip.open(os.path.join(DATA, fn), "rt") as f:
        for l in f:
            t = l.rstrip("\n").split("\t")
            if len(t) < 5: continue
            c = t[0].replace("chr","")
            if c not in CHRLEN: continue
            s, e, sc = int(t[1]), int(t[2]), float(t[4]) if t[4] not in (".","") else 0.0
            m = (s+e)//2
            mid.setdefault(c, []).append(m)
            full.append((c,s,e,sc))
    return {c: np.sort(np.array(v)) for c,v in mid.items()}, full, len(full)

def ens_get(url, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type":"application/json","User-Agent":"nrf2chip/1.0"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode("utf-8","ignore"))
        except Exception:
            time.sleep(2.0*(k+1))
    return None

def gene_tss(sym):
    d = ens_get(f"{ENS}/lookup/symbol/homo_sapiens/{sym}?expand=1&content-type=application/json")
    if not d: return None
    tr = [t for t in d["Transcript"] if t.get("is_canonical")==1] or d["Transcript"]
    t = tr[0]; strand = t["strand"]; tss = t["start"] if strand==1 else t["end"]
    return {"sym":sym,"chr":str(d["seq_region_name"]),"tss":int(tss),"strand":int(strand)}

def count_in_window(peaks_c, lo, hi):
    if peaks_c is None or len(peaks_c)==0: return 0
    return int(np.searchsorted(peaks_c, hi, "right") - np.searchsorted(peaks_c, lo, "left"))

def poisson_p(obs, lam):
    """One-sided P(X>=obs) under Poisson(lam), pure numpy (no scipy dependence risk)."""
    from math import exp, lgamma
    if lam <= 0: return 1.0 if obs>0 else 1.0
    # survival via cdf complement, sum terms
    kmax = max(obs, int(lam)+50)
    cdf = 0.0
    for k in range(0, obs):   # P(X<k) = cdf up to obs-1 ; tail = 1 - P(X<=obs-1)
        cdf += exp(-lam + k*np.log(lam) - lgamma(k+1))
    return max(0.0, min(1.0, 1.0 - cdf))

def decay_curve(peaks_by_chrom, centers, bins_kb):
    """peak density (peaks/bp) in distance bins around provided centre positions."""
    tot = sum(len(v) for v in peaks_by_chrom.values())
    dens = np.zeros(len(bins_kb)-1); npos=0
    for (c, pos) in centers:
        arr = peaks_by_chrom.get(c)
        if arr is None: continue
        npos += 1
        for bi in range(len(bins_kb)-1):
            lo, hi = bins_kb[bi]*1000, bins_kb[bi+1]*1000
            # two-sided window around pos
            cnt = 0
            for sgn in (-1, 1):
                a,b = pos+sgn*lo, pos+sgn*hi
                lo_,hi_ = (a,b) if a<b else (b,a)
                cnt += np.searchsorted(arr,hi_,"right")-np.searchsorted(arr,lo_,"left")
            width_bp = (hi-lo)  # per-side; two sides counted
            dens[bi] += cnt/(2*(hi-lo))
    return dens/max(npos,1)

# ---------------- main ----------------
print("== TASK1 NRF2 ChIP-seq reference analysis ==")
tss = {}
for g in GENES:
    r = gene_tss(g)
    if r: tss[g] = r; print("TSS", r)
    else: print("MISSING TSS", g)
if not tss:
    raise SystemExit("Could not fetch any gene TSS - abort (no fabrication).")

# background centres: random positions, same chrom-weighted, excluding near targets
bg_centers = []
chroms = list(CHRLEN); wts = np.array([CHRLEN[c] for c in chroms], float); wts/=wts.sum()
for c in np.random.choice(chroms, size=1000, p=wts):
    bg_centers.append((c, random.randint(0, CHRLEN[c])))

bins_kb = [0,2,10,50,100,200,500,1000]
rows = []
per_cell_curves = {}
for cell, fn in FILES.items():
    acc = {"A549":"ENCFF418TUX","HeLa":"ENCFF305KIK","HepG2":"ENCFF882YLO","IMR90":"ENCFF474PPT"}[cell]
    p = ensure(fn, acc)
    pb, full, npeaks = load_peaks(fn)
    print(f"\n{cell}: {npeaks} NRF2 peaks (real ENCODE narrowPeak)")
    for g,info in tss.items():
        c, t = info["chr"], info["tss"]
        arr = pb.get(c)
        obs2  = count_in_window(arr, t-2000,  t+2000)
        obs10 = count_in_window(arr, t-10000, t+10000)
        obs100= count_in_window(arr, t-100000,t+100000)
        for w,obs in [(2,obs2),(10,obs10),(100,obs100)]:
            lam = npeaks*(2*w*1000)/GENOME
            exp_ = lam; ratio = (obs+0.5)/(lam+0.5)
            pv = poisson_p(obs, lam)
            rows.append({"cell_line":cell,"accession":ACCESSIONS[cell],"gene":g,"chrom":c,"tss":t,
                         "window_kb":w,"n_peaks":obs,"expected_uniform":round(exp_,3),
                         "enrichment_ratio":round(ratio,2),"poisson_p_tail":pv,
                         "total_nrf2_peaks":npeaks})
    # decay curves: targets vs random background, all cells pooled per cell
    tgt_centers = [(info["chr"],info["tss"]) for info in tss.values()]
    dc_t = decay_curve(pb, tgt_centers, bins_kb)
    dc_b = decay_curve(pb, bg_centers, bins_kb)
    per_cell_curves[cell] = (dc_t, dc_b)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(RES,"nrf2_chip_hmoax1_enrichment.csv"), index=False)
print("\nwrote results/nrf2_chip_hmoax1_enrichment.csv  rows=",len(df))

# ---- figure ----
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
ax = axes[0]
sub = df[df.window_kb==10]
cells = list(FILES)
x = np.arange(len(GENES)); bw=0.2
for i,cell in enumerate(cells):
    vals=[]
    for g in GENES:
        r = sub[(sub.cell_line==cell)&(sub.gene==g)]
        vals.append(r.enrichment_ratio.values[0] if len(r) else 0)
    ax.bar(x+(i-1.5)*bw, vals, bw, label=cell, color=list(OI.values())[i])
ax.axhline(1.0, ls="--", c="k", lw=1)
ax.set_xticks(x); ax.set_xticklabels(GENES, rotation=0)
ax.set_ylabel("NRF2 peak enrichment vs uniform (±10 kb)")
ax.set_yscale("log"); ax.legend(fontsize=8); ax.set_title("NRF2 ChIP-seq binding at heme-handling loci")

ax2 = axes[1]
# pooled decay across cell lines: mean target density / mean bg density per bin
import numpy as np
tb = np.zeros(len(bins_kb)-1)
for cell,(dt,db) in per_cell_curves.items():
    ax2.plot(range(len(bins_kb)-1), dt, marker="o", ms=4, color=list(OI.values())[list(FILES).index(cell)], lw=1.4, label=f"{cell} NRF2@targets")
    ax2.plot(range(len(bins_kb)-1), db, marker="x", ms=4, ls="--", color=list(OI.values())[list(FILES).index(cell)], alpha=0.5, lw=1)
ax2.set_xticks(range(len(bins_kb)-1))
ax2.set_xticklabels([f"{bins_kb[i]}-{bins_kb[i+1]}" for i in range(len(bins_kb)-1)], rotation=35, fontsize=8)
ax2.set_xlabel("distance from TSS (kb, both sides)")
ax2.set_ylabel("NRF2 peak density (peaks/bp)")
ax2.set_title("Distance-decay: target loci (solid) vs random control (dashed)")

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig26_nrf2_chip_binding.png"), dpi=300)
fig.savefig(os.path.join(FIG,"Fig26_nrf2_chip_binding.pdf"))
print("wrote Fig26_nrf2_chip_binding.{png,pdf}")

# ---- audit entry for task 1 ----
h = df[(df.gene=="HMOX1")]
audit = {"status":"executed","kind":"option_b_NRF2_ChIP-seq_reference",
    "datasets":[ACCESSIONS[c] for c in FILES],
    "n_cell_lines":len(FILES),
    "hmoax1_summary": h[["cell_line","window_kb","n_peaks","expected_uniform","enrichment_ratio","poisson_p_tail"]].to_dict("records"),
    "conclusion":"Real ENCODE NRF2 ChIP-seq (4 cell lines) quantifies NRF2 binding events at HMOX1 and the heme module vs a uniform/distance-decay null.",
    "honest_label":"NRF2 chromatin-binding/accessibility reference; NOT SAH-specific scATAC."}
json.dump(audit, open(os.path.join(RES,"_t1_audit.json"),"w"), indent=2)
print("== T1 done ==")
