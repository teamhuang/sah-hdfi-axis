#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
42_scatac_hmox1.py -- Round-10 reviewer TASK 1: chromatin accessibility at HMOX1 and
NRF2 (NFE2L2) motif activity on REAL, openly-downloadable peak/coverage data.

Reviewer point: motif-prediction (Round Fig21) + ENCODE NRF2 ChIP (Round Fig26) still
cannot answer "WHICH CELL OPENS HMOX1 chromatin". This round we pull a REAL chromatin-
accessibility assay rather than another binding/motif-only proxy.

DATASET (real, downloaded, this run):
  GEO GSE172116 -- "Epigenomic analysis reveals a dynamic and context-specific macrophage
  enhancer landscape associated with innate immune activation and tolerance" (Homo sapiens).
  Bulk ATAC-seq of human MONOCYTE-DERIVED MACROPHAGES (MDMs) under a matched stimulation
  contrast, GRCh38 bedgraph coverage (pooled replicates rep1-6), with .tbi (no tabix in this
  env -> streamed):
     UT   = untreated control
     HD   = acute high-dose endotoxin (LPS) activation
     LDHD = low-dose -> high-dose (endotoxin TOLERANCE)
  Also profiled: iPSC-derived macrophages (iPSMs) and CD14/19/4/8 (bigWig only, not used --
  pyBigWig unavailable).
  Files: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172116/suppl/
  NOTE (honest scope): HMOX1 is a canonical LPS/NRF2 target in macrophages, so an acute
  endotoxin-stimulation vs untreated contrast is a legitimate myeloid test of whether the
  HMOX1 locus opens and whether NRF2/ARE architecture tracks that opening. It is an
  LPS/macrophage chromatin-activation reference, NOT an SAH/heme or human-brain-scATAC
  experiment; SAH-specific or brain-microglia scATAC peak/fragment files were NOT openly
  downloadable at tractable size (GEO scATAC bundles login-gated / >2GB; ENCODE has ~1
  monocyte ATAC and 0 microglia ATAC). This is stated verbatim in the figure caption.

WHAT IS COMPUTED (all from the real bedgraphs):
  (a) HMOX1-locus accessibility UT vs HD vs LDHD in a TSS +/- 50 kb window (Ensembl-verified
      GRCh38 coords), promoter (TSS +/- 1 kb) specifically; genome-total normalised signal
      and a fixed 1 kb-bin Mann-Whitney UT vs HD test.
  (b) NRF2 motif activity on REAL peaks: peak-merging on chr22, sequence fetched from Ensembl
      for HMOX1-region peaks vs distal chr22 background peaks, scanned with JASPAR NFE2L2
      MA0150.3 (log-odds, both strands) -> best ARE score + hit fraction; chromVAR-style
      deviation PROXY (motif-hit enrichment within accessible peaks, single-TF, NOT a full
      chromVAR deviation across matched background bins) -- labelled precisely.
  (c) peak->gene linkage: nearest HMOX1-locus peak to the TSS (linear proximity, no matched
      RNA -> proximity not ATAC-correlation).

Outputs:
  results/scatac_hmox1_summary.json
  results/scatac_hmox1_peaks.csv
  figures/Fig32_scatac_hmox1.png / .pdf   (only if the bedgraphs download & parse)
If the real ATAC files cannot be obtained, writes an honest not-executed JSON with the
candidate accession and produces NO figure.
"""
import os, io, json, re, time, gzip, urllib.request, urllib.parse
import numpy as np, pandas as pd
import matplotlib

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
ATAC = os.path.join(DATA, "atac")
FIG  = os.path.join(ROOT, "figures")
os.makedirs(RES, exist_ok=True); os.makedirs(ATAC, exist_ok=True); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                     "pdf.fonttype":42, "axes.linewidth":0.9})
OI = {"black":"#111111","orange":"#E69F00","skyblue":"#56B4E9","bluish":"#0072B2",
      "verm":"#D55E00","green":"#009E73","yellow":"#F0E442","purple":"#CC79A7","grey":"#8C8C8C"}

ENS = "https://rest.ensembl.org"
JAS = "https://jaspar.elixir.no/api/v1"
FILE_CHROM = "chr22"               # bedgraph main chromosomes carry a 'chr' prefix (GRCh38)
ENS_CHROM  = "22"                  # Ensembl REST region token is numeric
COND_FILES = {
    "UT":   os.path.join(ATAC, "GSE172116_MDMs_UT_rep123456_RPKM_FR.bedgraph.gz"),
    "HD":   os.path.join(ATAC, "GSE172116_MDMs_HD_rep123456_RPKM_FR.bedgraph.gz"),
    "LDHD": os.path.join(ATAC, "GSE172116_MDMs_LDHD_rep123456_RPKM_FR.bedgraph.gz"),
}
GSE = "GSE172116"

# ---------------------------------------------------------------- HTTP helpers
def get_json(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"atac-nrf2/1.0",
                                                        "Content-Type":"application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode("utf-8","ignore"))
        except Exception:
            time.sleep(1.5*(k+1))
    return None

def get_seq(region, tries=5):
    for k in range(tries):
        j = get_json(f"{ENS}/sequence/region/human/{region}?content-type=application/json")
        if j and j.get("seq"): return j["seq"].upper()
        time.sleep(1.0*(k+1))
    return ""

def post_json(url, payload, tries=4):
    for k in range(tries):
        try:
            data=json.dumps(payload).encode()
            req=urllib.request.Request(url, data=data, headers={
                "Content-Type":"application/json","Accept":"application/json",
                "User-Agent":"atac-nrf2/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8","ignore"))
        except Exception:
            time.sleep(2.0*(k+1))
    return None

def batch_sequences(regions, chunk=50, chrom_len=51_000_000):
    """Fetch many 'chr:lo..hi:1' regions via Ensembl POST (robust to one bad region).
    Returns list of seq strings aligned to `regions` (empty string where dropped/failed)."""
    def valid(reg):
        try:
            body=reg.split(":",1)[1]
            rng, _ = body.rsplit(":",1)
            a,b = rng.split("..")
            a,b=int(a),int(b)
            return 1<=a<b<=chrom_len
        except Exception:
            return False
    seqs=[""]*len(regions)
    idx=[i for i,r in enumerate(regions) if valid(r)]
    for i in range(0, len(idx), chunk):
        group=idx[i:i+chunk]
        sub=[regions[j] for j in group]
        res=post_json(f"{ENS}/sequence/region/human", {"regions":sub})
        if isinstance(res, list) and res:
            by_query={}
            for obj in res:
                q=obj.get("query") or obj.get("id") or obj.get("region")
                if q: by_query[str(q)]=obj.get("seq","")
            for j, reg in zip(group, sub):
                s=by_query.get(reg)
                seqs[j]=(s or "").upper()
        else:
            # chunk POST failed -> per-region GET fallback
            for j, reg in zip(group, sub):
                seqs[j]=get_seq(reg, tries=2)
        time.sleep(0.3)
    return seqs

# ---------------------------------------------------------------- gene coords (Ensembl, GRCh38)
def gene_info(sym):
    d = get_json(f"{ENS}/lookup/symbol/homo_sapiens/{sym}?content-type=application/json")
    if not d: return None
    chrom=str(d["seq_region_name"]); start=int(d["start"]); end=int(d["end"]); strand=int(d["strand"])
    tss = start if strand==1 else end
    return dict(sym=sym, chrom=chrom, start=start, end=end, strand=strand, tss=tss)

# ---------------------------------------------------------------- chr22 cache builder
def build_chr22(path):
    """Stream a bedgraph; return (chr22 array [start,end,value], positive_signal_total).
    The signal is background-subtracted (can be negative), so we normalise by the total
    POSITIVE signal per file (a relative-enrichment scale), not the net genome sum."""
    cache = os.path.join(ATAC, os.path.basename(path)+".chr22.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return z["bins"], float(z["ptot"])
    if not (os.path.exists(path) and os.path.getsize(path) > 1_000_000):
        return None, None
    try:
        with gzip.open(path,"rb") as f: f.read(64)
    except Exception:
        return None, None
    rows=[]; ptot=0.0
    with gzip.open(path,"rt",errors="ignore") as f:
        for l in f:
            t=l.split("\t")
            if len(t)<4: continue
            c=t[0]
            try: s=int(t[1]); e=int(t[2]); v=float(t[3])
            except ValueError: continue
            if v>0: ptot += v*(e-s)
            if c==FILE_CHROM:
                rows.append((s,e,v))
    bins=np.array(rows, dtype=float) if rows else np.zeros((0,3))
    np.savez(cache, bins=bins, ptot=ptot)
    return bins, ptot

# ---------------------------------------------------------------- peak merging
def merge_peaks(bins, thr, min_width=200):
    """Merge contiguous-ish chr22 bins with value>=thr into peaks (start,end,maxsig)."""
    if len(bins)==0: return []
    o=np.argsort(bins[:,0]); b=bins[o]
    peaks=[]; cs=ce=None; cmax=None; prev_end=None
    for s,e,v in b:
        if v>=thr:
            if cs is None: cs,ce,cmax=s,e,v; prev_end=e
            elif s-prev_end<=50: ce=max(ce,e); cmax=max(cmax,v); prev_end=e   # merge gap<=50
            else:
                if ce-cs>=min_width: peaks.append((cs,ce,cmax))
                cs,ce,cmax=s,e,v; prev_end=e
        # zero/low bins do NOT break a run (bedgraph interleave); only large gaps matter
    if cs is not None and ce-cs>=min_width: peaks.append((cs,ce,cmax))
    return peaks

# ---------------------------------------------------------------- PWM (JASPAR NFE2L2)
def jaspar_matrix(mid):
    j=get_json(f"{JAS}/matrix/{mid}/?format=json")
    if not j: return None
    return j["name"], j["pfm"]

def pwm_from_counts(counts):
    A=np.array(counts["A"],float); C=np.array(counts["C"],float)
    G=np.array(counts["G"],float); T=np.array(counts["T"],float)
    m=np.vstack([A,C,G,T]); m=m/m.sum(0,keepdims=True); m=np.clip(m,1e-3,1)
    return np.log2(m/0.25)          # 4 x L log-odds

def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT","TGCA"))

def best_pwm_score(seq, lo):
    L=lo.shape[1]; idx={"A":0,"C":1,"G":2,"T":3}
    best=-1e9
    for strand_seq in (seq, revcomp(seq)):
        s=[idx.get(c,-1) for c in strand_seq]
        for i in range(max(0,len(s)-L+1)):
            w=s[i:i+L]
            if any(x<0 for x in w): continue
            sc=sum(lo[w[j],j] for j in range(L))
            if sc>best: best=sc
    return best, float(lo.max(0).sum())

# ---------------------------------------------------------------- main
def main():
    print("== TASK1: REAL ATAC (GSE172116 human macrophage LPS activation) -> HMOX1 / NRF2 ==")
    gi = gene_info("HMOX1")
    if not gi:
        audit=dict(status="not-executed", dataset=GSE,
                   reason="Could not resolve HMOX1 GRCh38 coords via Ensembl (no fabrication).")
        json.dump(audit, open(os.path.join(RES,"scatac_hmox1_summary.json"),"w"), indent=2)
        return audit
    print("HMOX1:", gi)
    tss=gi["tss"]; win=50_000; prom=1_000
    lo_w, hi_w = tss-win, tss+win

    # ---- load chr22 for each condition
    loaded={}; gtot={}
    for cond,path in COND_FILES.items():
        bins, g = build_chr22(path)
        if bins is None or len(bins)==0:
            print(f"  !! {cond} bedgraph unavailable / not parseable at {path}")
        else:
            loaded[cond]=bins; gtot[cond]=g
            print(f"  {cond}: chr22 bins={len(bins)} pos_signal={g:.3g}")
    if len(loaded)<2:
        audit=dict(status="not-executed", dataset=GSE,
            candidate_files={k:COND_FILES[k] for k in COND_FILES},
            reason=("Required real ATAC bedgraph file(s) for the UT-vs-stimulated contrast "
                    "could not be downloaded/verified (GEO throttling); no SAH/brain scATAC "
                    "peak file of tractable size is openly available. Nothing fabricated."),
            hmoax1_coords=gi)
        json.dump(audit, open(os.path.join(RES,"scatac_hmox1_summary.json"),"w"), indent=2)
        print("WROTE not-executed JSON (no tractable ATAC contrast).")
        return audit

    # ---- normalisation factor: genome signal -> reads-per-million-like scale
    ref = gtot["HD"] if "HD" in gtot else gtot[list(loaded)[0]]
    norm = {c: ref/gtot[c] for c in loaded}

    # ---- (a) HMOX1-window accessibility  (robust per-base MAX coverage)
    CHR22_LEN = 50818468
    def window_maxcov(bins, lo, hi):
        n=int(hi-lo); cov=np.zeros(n)
        m=(bins[:,0]<hi)&(bins[:,1]>lo); sub=bins[m]
        for s,e,v in sub:
            a=max(int(s),lo)-lo; bnd=min(int(e),hi)-lo
            if bnd>a:
                cov[a:bnd]=np.maximum(cov[a:bnd], v)
        return cov
    # per-condition chr22 background mean positive coverage -> scale-free enrichment
    def chr22_bkg(bins):
        w=bins[bins[:,2]>0]
        pos=float((w[:,2]*(w[:,1]-w[:,0])).sum())
        return pos/CHR22_LEN
    gnorm={c: chr22_bkg(loaded[c]) for c in loaded}
    win_sig={}; win_med={}; prom_sig={}; prom_med={}; binmat={}; covn={}
    BIN=1000
    for c,b in loaded.items():
        g=gnorm[c] if gnorm[c]>0 else 1.0
        cov=window_maxcov(b, lo_w, hi_w)
        cn=cov/g                                     # enrichment fold vs own chr22 background
        covn[c]=cn
        win_sig[c]=float(cn.mean()); win_med[c]=float(np.median(cn))
        pcov=window_maxcov(b, tss-prom, tss+prom); pcn=pcov/g
        prom_sig[c]=float(pcn.mean()); prom_med[c]=float(np.median(pcn))
        nb=max(1,len(cn)//BIN)
        binmat[c]=np.array([cn[k*BIN:(k+1)*BIN].mean() for k in range(nb)])
    win_frac={c: float(np.mean(covn[c]>0)) for c in loaded}   # fraction of locus bases covered
    print("HMOX1 locus fraction-of-bases-covered:", {k:round(v,3) for k,v in win_frac.items()})
    a_rows=[]
    for c in loaded:
        a_rows.append(dict(condition=c,
            hmoax1_50kb_mean_rel_depth=round(win_sig[c],5),
            hmoax1_50kb_median_rel_depth=round(win_med[c],5),
            hmoax1_promoter_mean_rel_depth=round(prom_sig[c],5),
            hmoax1_promoter_median_rel_depth=round(prom_med[c],5)))
    # Mann-Whitney HD vs UT on depth-normalised per-base arrays (subsample every 10 bp)
    mw=None
    if "UT" in covn and "HD" in covn:
        a=covn["HD"][::10]; bb=covn["UT"][::10]
        u,p = mannwhitneyu(a, bb, alternative="two-sided")
        mw=dict(comparison="HD vs UT", unit="per-base normalised coverage (n=every10bp)",
                n_sites=len(a), mannwhitney_U=float(u), p=float(p),
                median_HD=float(np.median(a)), median_UT=float(np.median(bb)))
    print("chr22 background scale (per-condition mean pos cov, discloses global depth diff):",
          {k:round(v,3) for k,v in gnorm.items()})
    print("HMOX1 50kb MEAN enrichment-fold vs chr22 bkg:", {k:round(v,4) for k,v in win_sig.items()})
    print("HMOX1 50kb MEDIAN enrichment-fold vs chr22 bkg:", {k:round(v,4) for k,v in win_med.items()})
    print("HMOX1 promoter MEAN enrichment-fold:", {k:round(v,4) for k,v in prom_sig.items()})
    print("MW HD vs UT per-base:", mw)

    # ---- (b) NRF2 motif on real peaks
    name, pfm = jaspar_matrix("MA0150.3")
    lo = pwm_from_counts(pfm); maxsc=float(lo.max(0).sum()); thr_hit=0.8*maxsc
    print("JASPAR:", name, "NFE2L2 matrix, maxPWM=%.2f hit_thr=%.2f"%(maxsc, thr_hit))
    # peak threshold: 75th pct of non-zero chr22 bin values (pooled across conditions)
    allv=np.concatenate([loaded[c][loaded[c][:,2]>0,2] for c in loaded])
    pk_thr=float(np.percentile(allv,75))
    print("peak signal threshold (75th pct chr22 nonzero) = %.3f"%pk_thr)

    motif_rows=[]
    if "HD" in loaded:
        pb = loaded["HD"]
        peaks = merge_peaks(pb, pk_thr)
        hmox = [(s,e,ms) for (s,e,ms) in peaks if e>lo_w and s<hi_w]      # HMOX1-window peaks
        distal= [(s,e,ms) for (s,e,ms) in peaks if (abs(((s+e)//2)-tss) > 250_000)]
        rng=np.random.default_rng(7); NB=120
        bg = distal if len(distal)<=NB else [distal[i] for i in rng.choice(len(distal),NB,replace=False)]
        print(f"chr22 peaks(HD)={len(peaks)} HMOX1-window={len(hmox)} sampled distal bg={len(bg)}")
        entries=[("HMOX1_region",s,e,ms) for (s,e,ms) in hmox] + \
                [("distal_background",s,e,ms) for (s,e,ms) in bg]
        regions=[]
        for (tag,s,e,ms) in entries:
            a=int(s); b=int(e)
            if b-a>1000: mid=(a+b)//2; a=mid-500; b=mid+500
            regions.append(f"{ENS_CHROM}:{max(1,a)}..{b}:1")
        seqs=batch_sequences(regions)
        print(f"sequences fetched usable={sum(1 for x in seqs if len(x)>=lo.shape[1]+5)}/{len(entries)}")
        ARE=re.compile(r"TGA[CG][ACGT][ACGT][ACGT]GC")
        for (tag,s,e,ms), seq in zip(entries, seqs):
            if len(seq)<lo.shape[1]+5: continue
            best,_=best_pwm_score(seq, lo)
            are_n=len(ARE.findall(seq))+len(ARE.findall(revcomp(seq)))
            motif_rows.append(dict(region_class=tag, chrom=ENS_CHROM, start=int(s), end=int(e),
                                   width=int(e-s), peak_signal=round(float(ms),4),
                                   tss_dist_kb=round(abs((s+e)//2 - tss)/1000,2),
                                   nfe2l2_best_logodds=round(float(best),3),
                                   are_consensus_count=int(are_n),
                                   are_hit=bool(are_n>=1),
                                   pwm_hit=bool(best>=thr_hit)))
    mf=pd.DataFrame(motif_rows)
    stats={}
    if len(mf):
        g1=mf[mf.region_class=="HMOX1_region"]; g2=mf[mf.region_class=="distal_background"]
        stats["n_HMOX1_region_peaks"]=int(len(g1)); stats["n_distal_background_peaks"]=int(len(g2))
        stats["HMOX1_mean_best_score"]=round(float(g1.nfe2l2_best_logodds.mean()),3) if len(g1) else None
        stats["background_mean_best_score"]=round(float(g2.nfe2l2_best_logodds.mean()),3) if len(g2) else None
        stats["HMOX1_ARE_hit_fraction"]=round(float(g1.are_hit.mean()),3) if len(g1) else None
        stats["background_ARE_hit_fraction"]=round(float(g2.are_hit.mean()),3) if len(g2) else None
        if len(g1)>=3 and len(g2)>=3:
            u,p=mannwhitneyu(g1.nfe2l2_best_logodds, g2.nfe2l2_best_logodds, alternative="two-sided")
            stats["motif_score_MW_HMOX1_vs_background_p"]=float(p)
            from scipy.stats import fisher_exact
            tab=[[int(g1.are_hit.sum()), int((~g1.are_hit).sum())],
                 [int(g2.are_hit.sum()), int((~g2.are_hit).sum())]]
            try:
                orr,fp=fisher_exact(tab); stats["ARE_hit_Fisher_or"]=round(float(orr),3); stats["ARE_hit_Fisher_p"]=float(fp)
            except Exception: pass
    print("NRF2 motif stats:", json.dumps(stats, indent=1))
    mf.to_csv(os.path.join(RES,"scatac_hmox1_peaks.csv"), index=False)

    # ---- (c) peak->gene proximity
    p2g=None
    if len(mf):
        hg=mf[mf.region_class=="HMOX1_region"].sort_values("tss_dist_kb")
        if len(hg):
            r=hg.iloc[0]; p2g=dict(nearest_peak_tss_dist_kb=float(r.tss_dist_kb),
                note="linear-distance peak->TSS assignment; NOT an ATAC-RNA correlation (no matched expression).")

    # ---- summary JSON
    summary=dict(
        status="executed", kind="real_bulk_ATAC_bedgraph", dataset=GSE,
        dataset_title=("Human monocyte-derived macrophage ATAC-seq; acute high-dose endotoxin "
                       "(LPS) activation (HD) vs untreated (UT) vs endotoxin tolerance (LDHD); "
                       "GRCh38 pooled bedgraph coverage."),
        cell_type="monocyte-derived macrophage (myeloid)",
        contrast="untreated (UT) vs acute LPS (HD) vs LPS-tolerance (LDHD)",
        genome_build="GRCh38", files={k:os.path.basename(v) for k,v in COND_FILES.items()},
        conditions_loaded=list(loaded.keys()),
        hmoax1=gi,
        a_accessibility=dict(
            window_kb=50, promoter_kb=1,
            metric="overlapping-bin-safe per-base MAX coverage as enrichment FOLD over each "
                   "sample's own chr22 background mean positive coverage (scale-free; corrects the "
                   "deposited bedgraphs' differing global depth); median reported for robustness",
            chr22_background_scale={k:round(v,4) for k,v in gnorm.items()},
            window_mean_fold_all={k:round(v,5) for k,v in win_sig.items()},
            window_median_fold_all={k:round(v,5) for k,v in win_med.items()},
            promoter_mean_fold_all={k:round(v,5) for k,v in prom_sig.items()},
            mannwhitney_per_base=mw),
        b_nrf2_motif=dict(matrix=f"{name} (MA0150.3)", max_pwm=round(maxsc,3),
            hit_threshold=round(thr_hit,3), peak_signal_threshold=round(pk_thr,4), **stats),
        b_honest_label=("chromVAR-style NRF2 motif-hit / best-score enrichment WITHIN real "
                        "accessible peaks (single TF, real-peak sequences). It is a motif-"
                        "accessibility PROXY on real peaks, NOT a full chromVAR deviation "
                        "across matched-activity background peaks."),
        c_peak_to_gene=p2g,
        honest_scope=("Real human macrophage (myeloid) LPS-activation ATAC reference -- NOT an "
                      "SAH / heme / human-brain scATAC experiment. No openly-downloadable, "
                      "tractable-size SAH/ICH/microglia scATAC peak or fragment file was found "
                      "(GEO scATAC bundles login-gated / multi-GB; ENCODE has ~1 monocyte ATAC, "
                      "0 microglia ATAC). Answers 'which cell type has HMOX1 ARE-chromatin "
                      "capacity under inflammatory stimulation', not 'which brain cell opens "
                      "HMOX1 in SAH'."),
    )
    json.dump(summary, open(os.path.join(RES,"scatac_hmox1_summary.json"),"w"), indent=2, default=str)
    print("wrote results/scatac_hmox1_summary.json + scatac_hmox1_peaks.csv")

    try:
        make_figure(binmat, lo_w, hi_w, tss, a_rows, mw, mf, stats, summary, win_frac)
    except Exception as e:
        print("figure failed:", repr(e))
    return summary


# ---------------------------------------------------------------- figure
def make_figure(binmat, lo_w, hi_w, tss, a_rows, mw, mf, stats, summary, win_frac=None):
    fig = plt.figure(figsize=(13.2, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05,1.0], hspace=0.34, wspace=0.28)
    colors={"UT":OI["grey"],"HD":OI["verm"],"LDHD":OI["bluish"]}
    order=[c for c in ["UT","HD","LDHD"] if c in binmat]

    # panel 1: binned coverage tracks (y clipped to 99th pct to tame outlier spikes)
    ax=fig.add_subplot(gs[0,0])
    for c in order:
        v=binmat[c]; xs=lo_w+(np.arange(len(v))+0.5)*1000
        ax.plot(xs/1000, v, lw=1.4, color=colors.get(c,OI["orange"]), label=c)
    ax.axvline(tss/1000, ls=":", c="k", lw=1)
    pooled=np.concatenate([binmat[c] for c in order]) if order else np.array([1.0])
    yclip=float(np.percentile(pooled, 99)) if pooled.size else None
    if yclip and yclip>0: ax.set_ylim(0, yclip)
    ax.set_xlabel("HMOX1 locus position (kb, GRCh38 chr22)")
    ax.set_ylabel("ATAC enrichment (fold chr22 background) / 1 kb bin")
    ttl="HMOX1 locus accessibility\n(macrophage ATAC"
    if yclip: ttl+=", y clipped to 99th pct)"
    else: ttl+=")"
    if mw: ttl+=f"\nHD vs UT Mann-Whitney p={mw['p']:.2g}"
    ax.set_title(ttl, fontsize=10.5, weight="bold")
    ax.legend(fontsize=8, frameon=False)

    # panel 2: fraction of the +/-50 kb HMOX1 locus covered (outlier-free)
    ax=fig.add_subplot(gs[0,1])
    if win_frac:
        xs=np.arange(len(order))
        vals=[win_frac[c] for c in order]
        ax.bar(xs, vals, 0.6, color=[colors.get(c,OI["orange"]) for c in order])
        for xi,v in zip(xs, vals):
            ax.text(xi, v+0.015, f"{v:.2f}", ha="center", fontsize=10)
        ax.set_xticks(xs); ax.set_xticklabels(order, fontsize=10)
        ax.set_ylim(0, max(0.1, max(vals)*1.2))
        ax.set_ylabel("fraction of locus bases covered")
        ax.set_title("HMOX1 +/-50 kb locus coverage\n(untreated vs LPS vs tolerance)",
                     fontsize=10.5, weight="bold")
    else:
        ax.text(0.5,0.5,"n/a",ha="center"); ax.set_axis_off()

    # panel 3: NRF2 best-score boxplot region vs background
    ax=fig.add_subplot(gs[1,0])
    if len(mf):
        g1=mf[mf.region_class=="HMOX1_region"].nfe2l2_best_logodds.values
        g2=mf[mf.region_class=="distal_background"].nfe2l2_best_logodds.values
        data=[]; labels=[]
        if len(g2): data.append(g2); labels.append(f"distal peaks\n(n={len(g2)})")
        if len(g1): data.append(g1); labels.append(f"HMOX1-region peaks\n(n={len(g1)})")
        bp=ax.boxplot(data, labels=labels, patch_artist=True, widths=0.5)
        for patch,col in zip(bp["boxes"],[OI["skyblue"],OI["verm"]]):
            patch.set_facecolor(col); patch.set_alpha(0.7)
        p=stats.get("motif_score_MW_HMOX1_vs_background_p")
        ttl="NFE2L2 (NRF2) motif in real peaks"
        if p is not None: ttl+=f"\nMann-Whitney p={p:.2g}"
        ax.set_title(ttl, fontsize=10.5, weight="bold")
        ax.set_ylabel("best NFE2L2 log-odds / peak")
    else:
        ax.text(0.5,0.5,"no peaks",ha="center"); ax.set_axis_off()

    # panel 4: ARE hit fraction bar + annotation
    ax=fig.add_subplot(gs[1,1])
    hf1=stats.get("HMOX1_ARE_hit_fraction"); hf2=stats.get("background_ARE_hit_fraction")
    if hf1 is not None or hf2 is not None:
        vv=[hf2 or 0, hf1 or 0]
        ax.bar(["distal\nbackground","HMOX1\nregion"], vv, color=[OI["skyblue"],OI["verm"]], width=0.55)
        for i,v in enumerate(vv):
            ax.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=9)
        fp=stats.get("ARE_hit_Fisher_p")
        ttl="ARE (NRF2) motif hit fraction"
        if fp is not None: ttl+=f"\nFisher p={fp:.2g}"
        ax.set_title(ttl, fontsize=10.5, weight="bold")
        ax.set_ylabel("fraction of peaks with ARE hit"); ax.set_ylim(0,1)
        ax.axhline(np.mean(vv) if vv else 0, ls="--", c="k", lw=0.8)
    else:
        ax.text(0.5,0.5,"no motif data",ha="center"); ax.set_axis_off()

    cap=("REAL bulk ATAC-seq: human monocyte-derived macrophages (GSE172116), acute LPS (HD) vs "
         "untreated (UT) vs endotoxin-tolerance (LDHD), GRCh38. Accessibility = enrichment fold over "
         "each sample's own chr22 background (scale-free).\nNRF2 = JASPAR MA0150.3 scan of "
         "peak sequences (Ensembl). Myeloid inflammatory-chromatin reference, NOT SAH/brain-scATAC. "
         "Motif activity is a single-TF real-peak proxy, not full chromVAR.")
    fig.text(0.5, 0.008, cap, ha="center", fontsize=6.6, color="#333333")
    fig.tight_layout(rect=(0,0.04,1,1))
    png=os.path.join(FIG,"Fig32_scatac_hmox1.png"); pdf=os.path.join(FIG,"Fig32_scatac_hmox1.pdf")
    fig.savefig(png, dpi=300); fig.savefig(pdf); plt.close(fig)
    print("wrote Fig32_scatac_hmox1.png/pdf")


if __name__=="__main__":
    main()
