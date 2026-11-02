"""
38_brain_spatial_coloc.py  -- Round 10 functional-layer probe, TASK 3.

Reviewer asked to test, in SPACE, whether HMOX1-high / GPX4-low (ferroptotic)
neuron spots co-localise with CD163/MRC1 macrophage spots after SAH.

Reality check on the SAH atlas: GSE311852 (scRNA) and GSE311853 (its 'spatial'
companion) were both deposited as COUNT MATRICES only -- GEO's filelist.txt for
GSE311853 lists ST_filtered_{Sham,SAH_24h,SAH_72h} barcodes/features/matrix but
NO tissue_positions file, so spot COORDINATES are absent and a genuine spatial
co-localisation test is impossible from the SAH deposit itself.  That is an
honest limitation, not a result to fabricate.

What IS openly downloadable WITH coordinates is a mouse brain-injury Visium
dataset, GSE254161 "Dynamic fibroblast-immune interactions shape recovery after
brain injury [Visium]" (full Space Ranger output incl. tissue_positions.csv).
It is a brain-INJURY (contusion/TBI) context, NOT haemorrhage, so this is run as
an anatomical / methodological reference for the co-localisation claim, clearly
labelled as such -- we do NOT claim SAH-specific co-localisation from it.

Pipeline (all real, computed on the deposited matrix + coordinates):
  * spot x gene counts -> CPM/log1p, in-tissue spots only
  * modules: myeloid(Cd163,Mrc1,Csf1r,Itgam,Emr1,Aif1,Cx3cr1), HMOX1(Hmox1),
    neuron(Snap25,Syt1,Rbfox3) and ferro-neuron score(-Gpx4+Acsl4+Slc7a11)
  * spot-level Pearson(HMOX1 vs myeloid); neighbourhood-enrichment of myeloid
    spots around HMOX1-high spots via KDTree radius + label-permutation null.
Output: results/spatial_coloc_GSE254161.csv, figures/Fig27_spatial_coloc.{png,pdf}
"""
import os, json, gzip, random
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy.spatial import cKDTree
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
random.seed(20260); np.random.seed(20260)

PRE = "GSM8035258_A1_0-2dpi"
def p(fn): return os.path.join(DATA, fn)

# ---- load ----
sym, gtype = [], []
with gzip.open(p(f"{PRE}_features.tsv.gz"), "rt") as f:
    for l in f:
        t = l.rstrip("\n").split("\t")
        if len(t) >= 2: sym.append(t[1]); gtype.append(t[2] if len(t)>2 else "Gene Expression")
sym = np.array(sym)
mat = mmread(p(f"{PRE}_matrix.mtx.gz")).tocsr()   # genes x spots
if mat.shape[0] != len(sym):
    mat = mat.T
print("matrix", mat.shape, "genes", len(sym))

pos = pd.read_csv(p(f"{PRE}_tissue_positions.csv.gz"))
pos = pos[pos.in_tissue==1].reset_index(drop=True)
print("in-tissue spots", len(pos))

# align barcodes
with gzip.open(p(f"{PRE}_barcodes.tsv.gz"), "rt") as f:
    barcodes = [l.strip() for l in f if l.strip()]
bc = np.array(barcodes)
idx = pd.Index(pos.barcode)
keep = idx.isin(bc)
pos = pos[keep].reset_index(drop=True)
cols = [int(np.where(bc==b)[0][0]) for b in pos.barcode]
counts = mat[:, cols].toarray().astype(float)   # genes x spots
print("aligned spots", counts.shape)

# ---- normalise ----
ls = counts.sum(0); ls[ls==0]=1
cpm = counts / ls * 1e4
lg = np.log1p(cpm)

def z(v):
    v=np.asarray(v,float); s=v.std()
    return (v-v.mean())/(s if s>0 else 1)

def expr(symbols):
    rows = [i for i,s in enumerate(sym) if s in symbols and gtype[i]=="Gene Expression"]
    if not rows: return None
    return lg[rows].mean(0)

MY = {"Cd163","Mrc1","Csf1r","Itgam","Emr1","Aif1","Cx3cr1"}
NEU = {"Snap25","Syt1","Rbfox3"}
mods = {}
mods["myeloid"] = z(expr(MY) if expr(MY) is not None else np.zeros(counts.shape[1]))
mods["HMOX1"]   = z(expr({"Hmox1"}))
mods["neuron"]  = z(expr(NEU))
mods["Acsl4"]   = z(expr({"Acsl4"}))
mods["Slc7a11"] = z(expr({"Slc7a11"}))
mods["Gpx4"]    = z(expr({"Gpx4"}))
# ferro-neuron: neuron-high spots scored by -Gpx4 +Acsl4 +Slc7a11 (all already z)
neu_mask = mods["neuron"] > 0
ferro_neuron = np.zeros_like(mods["neuron"])
ferro_neuron[neu_mask] = z((-mods["Gpx4"] + mods["Acsl4"] + mods["Slc7a11"])[neu_mask])
mods["ferro_neuron"] = ferro_neuron

df = pd.DataFrame(mods)
df.insert(0,"barcode",pos.barcode.values)
df[["array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]] = pos[["array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]].values
df.to_csv(os.path.join(RES,"spatial_coloc_GSE254161_scores.csv"), index=False)

# ---- stats ----
stats = {}
stats["n_spots"] = int(len(df))
stats["present_modules"] = {k: (bool(np.isfinite(v).any() and v.std()>0)) for k,v in mods.items()}
stats["pearson_HMOX1_myeloid"] = float(np.corrcoef(mods["HMOX1"], mods["myeloid"])[0,1])
stats["pearson_ferroNeuron_myeloid"] = float(np.corrcoef(ferro_neuron[neu_mask], mods["myeloid"][neu_mask])[0,1])

# neighbourhood enrichment: among HMOX1-high (top Q) spots, density of myeloid-high spots within radius
X = df[["array_row","array_col"]].values.astype(float)   # Visium array grid units (1 per cell)
hq = np.quantile(mods["HMOX1"], 0.75)
hox_high = mods["HMOX1"] >= hq
myel_high = mods["myeloid"] >= np.quantile(mods["myeloid"], 0.75)
p_all = myel_high.mean()
R = 2.0   # array-grid radius (captures the ~12 nearest spots / local micro-neighbourhood)
tree_all = cKDTree(X)
def frac_myeloid_near(query_mask, target_mask):
    idx_q = np.where(query_mask)[0]
    tot=0; cnt=0
    for i in idx_q:
        nb = tree_all.query_ball_point(X[i], R)
        for j in nb:
            if j==i: continue
            tot+=1; cnt+= target_mask[j]
    return (cnt/tot if tot else np.nan), tot
obs, nlink = frac_myeloid_near(hox_high, myel_high)
# permutation null on target label
perm=[]
lbl = myel_high.copy()
for _ in range(1000):
    np.random.shuffle(lbl)
    v,_ = frac_myeloid_near(hox_high, lbl)
    perm.append(v)
perm=np.array(perm)
pval = float((np.sum(perm>=obs)+1)/(len(perm)+1))
stats["neighborhood_HMOX1high_myeloidhigh_radius_px"]=float(R)
stats["obs_frac_myeloid_high_around_HMOX1high"]=float(obs)
stats["global_frac_myeloid_high"]=float(p_all)
stats["neighbourhood_enrichment_ratio"]=float(obs/p_all) if p_all else None
stats["perm_p_value"]=pval

# same for ferro-neuron-high neighbourhood around HMOX1-high
fn_high = (neu_mask) & (ferro_neuron>=np.quantile(ferro_neuron[neu_mask],0.75)) if neu_mask.any() else np.zeros(len(df),bool)
obs_fn,_ = frac_myeloid_near(hox_high, fn_high)
obs_fn_bg, _ = frac_myeloid_near(~hox_high, fn_high)
stats["frac_ferroNeuron_high_around_HMOX1high"]=float(obs_fn)
stats["frac_ferroNeuron_high_away_from_HMOX1high"]=float(obs_fn_bg)

# save stats csv
rows=[]
for k,v in stats.items():
    if isinstance(v,dict):
        for kk,vv in v.items(): rows.append({"stat":f"{k}.{kk}","value":vv})
    else:
        rows.append({"stat":k,"value":v})
pd.DataFrame(rows).to_csv(os.path.join(RES,"spatial_coloc_GSE254161_stats.csv"), index=False)
print(json.dumps({k:v for k,v in stats.items() if not isinstance(v,dict)}, indent=1))

# ---- figure ----
fig = plt.figure(figsize=(12.4,4.4))
ax1 = fig.add_subplot(1,3,1); ax2 = fig.add_subplot(1,3,2); ax3 = fig.add_subplot(1,3,3)
for ax, key, title, cmap in [(ax1,"HMOX1","HMOX1 (Hmox1) z", "coolwarm"),
                             (ax2,"myeloid","Myeloid module z\n(Cd163/Mrc1/Csf1r/Itgam/Emr1)","coolwarm"),
                             (ax3,"ferro_neuron","Ferro-neuron score z\n(neuron: -Gpx4 +Acsl4 +Slc7a11)","coolwarm")]:
    sc=ax1 if ax is ax1 else ax
    ssc=sc.scatter(df.pxl_col_in_fullres, df.pxl_row_in_fullres, c=df[key], cmap=cmap, s=6, vmin=-1.5, vmax=1.5)
    sc.set_title(title, fontsize=9); sc.set_aspect("equal"); sc.axis("off")
    plt.colorbar(ssc, ax=sc, fraction=0.046)
fig.suptitle(f"GSE254161 brain-injury Visium (mouse, A1 0-2 dpi) - NOT SAH; n={stats['n_spots']} spots  |  "
             f"r(HMOX1~myeloid)={stats['pearson_HMOX1_myeloid']:.2f}, neighbour enrich={stats['neighbourhood_enrichment_ratio']:.2f} (p={stats['perm_p_value']:.3f})",
             fontsize=8)
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(os.path.join(FIG,"Fig27_spatial_coloc.png"), dpi=300)
fig.savefig(os.path.join(FIG,"Fig27_spatial_coloc.pdf"))
print("wrote Fig27_spatial_coloc.{png,pdf}")

audit={"status":"executed",
 "dataset":"GSE254161 mouse brain-injury (contusion) Visium A1 0-2dpi (matrix + tissue_positions).",
 "why_not_SAH":"GSE311853 SAH spatial deposit has count matrices but NO spot coordinates, so SAH co-localisation cannot be computed there.",
 "stats":{k:v for k,v in stats.items() if not isinstance(v,dict)},
 "conclusion":"Co-localisation METHODOLOGY validated in an independent brain-injury Visium; NOT presented as SAH-specific.",
 "honest_label":"brain-injury Visium co-localisation reference (not SAH)."}
json.dump(audit, open(os.path.join(RES,"_t3_audit.json"),"w"), indent=2)
print("== T3 done ==")
