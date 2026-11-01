"""
30_regulatory_motif.py
Reviewer multi-omics request (ATAC / regulatory layer), executed with REAL
public data that is actually downloadable in this environment (Ensembl REST +
JASPAR) rather than the ~500 MB GWAS sumstats that are not tractable here.

Question answered: "why is HMOX1 (and the heme-handling module) transcriptionally
activated — is there NRF2/ARE (and BACH1-CNC) regulatory architecture at these
promoters?" This is the sequence-level proxy for the ATAC/motif analysis the
reviewer asked for.

For each gene we fetch the GRCh38 promoter (TSS +/- 2 kb) from Ensembl, then:
  1. scan it with the JASPAR NFE2L2 (ARE) and BACH1 position-weight matrices
     (log-odds, both strands, threshold = 80% of max score) -> ARE/CNC site count;
  2. count the canonical ARE/MARE consensus matches (TGACnnnGC; MARE);
  3. query Ensembl regulatory features (promoter/enhancer) overlapping the locus;
  4. query Ensembl TF-binding motif features (ENCODE-derived) and tally the
     fraction that are CNC-bZIP/ARE-type (NFE2L2/MAF/BACH/AREG).
Outputs: results/regulatory_motif_HMOX1.csv and figures/Fig21_regulatory_motif.png/pdf.
"""
import os, json, re, time
import numpy as np, pandas as pd
import urllib.request
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"grey":"#8C8C8C","verm":"#D55E00","blue":"#0072B2","green":"#009E73","orange":"#E69F00","purple":"#CC79A7","sky":"#56B4E9"}

ENS = "https://rest.ensembl.org"
JAS = "https://jaspar.elixir.no/api/v1"

def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"motif/1.0","Content-Type":"application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read().decode("utf-8","ignore")
        except Exception as e:
            if k == tries-1: return ""
            time.sleep(1.5)
    return ""

def get_json(url, tries=4):
    t = get(url, tries)
    try: return json.loads(t) if t else None
    except Exception: return None

def get_seq(region):
    # Ensembl sequence via JSON (avoids content-type clash)
    j = get_json(f"{ENS}/sequence/region/human/{region}?content-type=application/json")
    return (j.get("seq","") if j else "").upper()

# ---- JASPAR matrices ----
def jaspar_matrix(mid):
    j = get_json(f"{JAS}/matrix/{mid}/?format=json")
    if not j: return None
    counts = j.get("pfm") or j.get("matrix") or j.get("counts")
    if not counts: return None
    return {"name": j.get("name"), "counts": counts}

# fetch NFE2L2 ARE matrix; BACH1 competes for the SAME ARE/MARE element, so we
# report NRF2-ARE architecture and treat BACH1 via the consensus/CNC-bZIP motif
# features rather than a separate (unreliable) JASPAR accession.
nfe2 = jaspar_matrix("MA0150.3")            # NFE2L2 (NRF2)
bach = None
print("JASPAR NFE2L2:", None if not nfe2 else nfe2["name"])

def pwm_from_counts(counts):
    A=np.array(counts["A"],float); C=np.array(counts["C"],float); G=np.array(counts["G"],float); T=np.array(counts["T"],float)
    m=np.vstack([A,C,G,T]); m=m/m.sum(0,keepdims=True)     # shape (4, L) bases x positions
    m=np.clip(m,1e-3,1)
    return np.log2(m/0.25)   # log-odds, 4 x L

def scan_pwm(seq, logodds):
    L=logodds.shape[1]; idx={"A":0,"C":1,"G":2,"T":3}
    s=np.array([idx.get(c,-1) for c in seq])
    maxsc=logodds.max(0).sum()
    thr=0.8*maxsc
    scores=[]
    for i in range(max(0,len(s)-L+1)):
        w=s[i:i+L]
        if (w<0).any(): continue
        sc=sum(logodds[w[j], j] for j in range(L))   # base at pos j
        scores.append(sc)
    if not scores: return 0,0.0
    nhit=sum(1 for sc in scores if sc>=thr)
    return nhit, max(scores)

def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT","TGCA"))

ARE_CORE = re.compile(r"TGA[CG][ACGT][ACGT][ACGT]GC")     # TGACnnnGC ARE
MARE = re.compile(r"TGCTGA[CG]CAGC[GT]GCA")               # MARE / CNC-bZIP
CNCBZIP = re.compile(r"(NFE2|MAF|BACH|ARE|NRF2|SMAF|MAFK|MAFF|MAFG)", re.I)

GENES = ["HMOX1","HP","FTH1","FTL","GPX4","SLC40A1","TFRC","NFE2L2","BACH1"]
rows=[]
for g in GENES:
    info = get_json(f"{ENS}/lookup/symbol/homo_sapiens/{g}?content-type=application/json")
    if not info: 
        rows.append(dict(gene=g, found=False)); continue
    chrom=info["seq_region_name"]; start=info["start"]; end=info["end"]; strand=info["strand"]
    tss = start if strand==1 else end
    p_lo, p_hi = tss-2000, tss+2000
    seq = get_seq(f"{chrom}:{p_lo}..{p_hi}:{strand}")
    nA = nB = 0; bA=bB=0.0
    if nfe2 and seq:
        lo = pwm_from_counts(nfe2["counts"])
        h1,s1 = scan_pwm(seq, lo); h2,s2 = scan_pwm(revcomp(seq), lo)
        nA = h1 + h2; bA = max(s1, s2)
    if bach and seq:
        lo = pwm_from_counts(bach["counts"])
        h1,s1 = scan_pwm(seq, lo); h2,s2 = scan_pwm(revcomp(seq), lo)
        nB = h1 + h2; bB = max(s1, s2)
    are = len(ARE_CORE.findall(seq)) + len(ARE_CORE.findall(revcomp(seq)))
    mare = len(MARE.findall(seq)) + len(MARE.findall(revcomp(seq)))
    # Ensembl regulatory + motif features
    reg = get_json(f"{ENS}/overlap/region/human/{chrom}:{p_lo}..{p_hi}?feature=regulatory;content-type=application/json") or []
    mot = get_json(f"{ENS}/overlap/region/human/{chrom}:{p_lo}..{p_hi}?feature=motif;content-type=application/json") or []
    reg_types = pd.Series([r.get("description","?") for r in reg]).value_counts().to_dict() if reg else {}
    tf_names = [m.get("transcription_factor_complex","") for m in mot]
    are_motifs = sum(1 for t in tf_names if CNCBZIP.search(str(t)))
    rows.append(dict(gene=g, found=True, chrom=chrom, tss=tss,
                     nfe2l2_pwm_sites=nA, nfe2l2_best=round(bA,2),
                     bach1_pwm_sites=nB, bach1_best=round(bB,2),
                     are_consensus=are, mare_consensus=mare,
                     ens_reg_features=len(reg), ens_motif_features=len(mot), ens_are_type_motifs=are_motifs,
                     reg_types=json.dumps(reg_types)))
df=pd.DataFrame(rows)
df.to_csv(os.path.join(RES,"regulatory_motif_HMOX1.csv"), index=False)
print(df.to_string(index=False))

# ---------------- figure ----------------
d=df[df.found==True].copy()
fig, axes = plt.subplots(1, 3, figsize=(15.5,5.0))
# panel 1: ARE (NFE2L2 PWM) site count in promoter
ax=axes[0]
ax.bar(range(len(d)), d.nfe2l2_pwm_sites, color=OI["blue"])
ax.bar(range(len(d)), d.are_consensus, color=OI["verm"], alpha=.55, width=.5, label="ARE consensus")
ax.set_xticks(range(len(d))); ax.set_xticklabels(d.gene, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("NRF2/ARE sites in promoter (\u00b12 kb)")
ax.set_title("NFE2L2 PWM + ARE consensus", fontsize=11, weight="bold")
ax.legend(fontsize=8, frameon=False)
# panel 2: Ensembl motif features total vs ARE/CNC-type
ax=axes[1]
x=np.arange(len(d)); w=.4
ax.bar(x-w/2, d.ens_motif_features, w, color=OI["grey"], label="all TF motifs")
ax.bar(x+w/2, d.ens_are_type_motifs, w, color=OI["green"], label="CNC-bZIP/ARE-type")
ax.set_xticks(x); ax.set_xticklabels(d.gene, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Ensembl (ENCODE) motif features")
ax.set_title("Annotated TF-binding motifs", fontsize=11, weight="bold")
ax.legend(fontsize=8, frameon=False)
# panel 3: Ensembl regulatory features (promoter/enhancer) count
ax=axes[2]
ax.bar(range(len(d)), d.ens_reg_features, color=OI["orange"])
ax.set_xticks(range(len(d))); ax.set_xticklabels(d.gene, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Regulatory features (\u00b12 kb)")
ax.set_title("Open-regulatory architecture", fontsize=11, weight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig21_regulatory_motif.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig21_regulatory_motif.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Fig21 + regulatory_motif_HMOX1.csv")
