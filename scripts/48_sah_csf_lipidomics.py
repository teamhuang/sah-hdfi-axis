"""
48_sah_csf_lipidomics.py -- Round-11 reviewer TASK: LIPID-LEVEL corroboration of the
HDFI -> ferroptosis arm using REAL human SAH cerebrospinal-fluid lipidomics.

The reviewer's PhD plan asks for a "ferroptosis lipid index" (PUFA-PL / oxidised PE /
arachidonate) measured in human SAH, cross-linked to HDFI.

HARD-HONESTY RULE: only real, openly obtainable data. Nothing is fabricated, imputed
or invented. Where a required quantitative object cannot be obtained we say so.

DATA-ACCESS REALITY ESTABLISHED IN THIS RUN
-------------------------------------------
The one openly described human SAH-CSF lipidomics study is
  PMC10525923 / Biomedicines 2023 (11,2360), "Distinct Cerebrospinal Fluid Lipid
  Signature in Patients with SAH-Induced Hydrocephalus" (SAH n=13 vs control
  n=11, cisternal CSF, LC-MS/MS, 244 classified lipids).
Its deposited code repository (the Data-Availability link,
  https://github.com/Sorennorge/MacAulayLab-SAH-Metabolomics) contains ONLY the
  analysis scripts; the per-sample normalised abundance matrices are read from
  local *.xlsx inputs that are NOT published in the repository, and the journal
  supplement is a PDF of figures with NO numeric table.
  => The per-sample human-CSF abundance matrix is therefore NOT openly downloadable
     and we CANNOT rebuild the volcano / recompute per-sample fold-changes. We do not
     fake them.

What we CAN do honestly:
  (A) transcribe, verbatim, the authors' reported differential-lipid findings
      (species, class, direction, and whether they cleared multiple-testing
      correction) -- read directly from the open-access main text (_tmp XML we
      fetched); these are the authors' numbers, never re-derived here; and
  (B) anchor the *quantitative* ferroptosis-execution lipid evidence on the REAL,
      numeric, openly downloadable datasets used in Round 10
      (Metabolomics Workbench ST003337 valinomycin-ferroptosis MEFs and
      ST004928 RSL3/GPX4-inhibition MDSCs, results/lipidomics_ferroptosis.csv),
      which supply the actual ox-PE fold-changes with CIs.

So the "ferroptosis lipid" claim rests on (B) real execution-layer numbers, while
(A) is a directional human-CSF corroboration, explicitly labelled verbatim.

Outputs:
  results/tables/sah_csf_lipidomics.csv   (verbatim human findings + numeric mouse anchor)
  results/sah_csf_lipidomics_summary.json
  figures/Fig38_sah_csf_lipidomics.png + .pdf  (Okabe-Ito, Arial, 300 dpi, no panel letters)
"""
import os, re, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB  = os.path.join(ROOT, "results", "tables")
FIG  = os.path.join(ROOT, "figures")
for d in (TAB, FIG): os.makedirs(d, exist_ok=True)
CSV10 = os.path.join(ROOT, "results", "lipidomics_ferroptosis.csv")   # Round-10 numeric
XMLPATH = os.path.join(ROOT, "_tmp_probe", "pmc10525923.xml")         # open-access XML

OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", grey="#555555", yellow="#F0E442")
rcParams.update({"font.family":"Arial","font.size":9,"axes.linewidth":0.9,
                 "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,
                 "axes.spines.right":False,"pdf.fonttype":42})

COHORT = "PMC10525923 (Biomedicines 2023;11:2360): human cisternal CSF, " \
         "SAH n=13 vs control n=11, LC-MS/MS, 244 classified lipids."

# -----------------------------------------------------------------------------
# (A) VERBATIM authors' differential-lipid findings (transcribed, not re-derived)
# -----------------------------------------------------------------------------
def verbatim_human_findings():
    # The species/values below are READ from the published open-access main text.
    # We verify the key tokens are actually present in the XML so nothing is invented.
    xml = open(XMLPATH, encoding="utf-8").read() if os.path.exists(XMLPATH) else ""
    txt = re.sub("<[^>]+>", " ", xml); txt = re.sub("  +", " ", txt)
    def present(tok): return bool(re.search(re.escape(tok), txt))
    # FDR-significantly up in SAH (the paper's Figure-4C corrected set, n=7 species / 4 classes)
    fdr = [
        ("Arachidonoyl amide", "Amide",            "up", True,  "PUFA (20:4)-derived; ferroptosis-relevant"),
        ("Amide C18",          "Amide",            "up", True,  "fatty-acide amide"),
        ("PC 30:0",            "Phosphatidylcholine","up",True, ""),
        ("PC 31:0",            "Phosphatidylcholine","up",True, ""),
        ("PC 32:0",            "Phosphatidylcholine","up",True, ""),
        ("PC (16:0/16:0)",     "Phosphatidylcholine","up",True, ""),
        ("Plasmenyl-PE 43:6",  "plasmenyl-PE (ether PE)","up",True, "ether-PE; PUFA-PL class"),
        ("Plasmenyl-PE 36:1 (18:0/18:1)","plasmenyl-PE (ether PE)","up",True,"ether-PE"),
    ]
    # up but significant ONLY WITHOUT multiple-testing correction (>2.5 log2FC; 5-16x)
    uncorr = [
        ("Fatty acids (class cluster)","Fatty acid","up (uncorrected)",True,">2.5 log2FC (5-16x)"),
        ("Lysophosphatidylcholines (LPC)","LPC","up (uncorrected)",True,">2.5 log2FC (5-16x)"),
        ("Monoacylglycerols (MAG)","MAG","up (uncorrected)",True,">2.5 log2FC (5-16x)"),
        ("Phosphatidic acids (PA22:0,24:0,26:0,28:0,45:4)","Phosphatidic acid","up (uncorrected)",True,"5 PA species; >2.5 log2FC (5-16x)"),
    ]
    down = [
        ("Some amide-class lipids","Amide","down",True,"authors: amides partly downregulated"),
    ]
    rows = []
    for name, cls, direction, checked, note in fdr + uncorr + down:
        rows.append(dict(layer="Human SAH-CSF lipidomics (VERBATIM published values)",
                         cohort=COHORT, analyte=name, lipid_class=cls, direction=direction,
                         significance=("FDR-corrected" if "uncorrected" not in direction and direction=="up"
                                       else ("nominal only (no FDR)" if "uncorrected" in direction else "as reported")),
                         fold_range="", ferroptosis_note=note,
                         provenance="transcribed from open-access main text (authors' Figure 4C / results); NOT re-derived here",
                         token_present_in_xml=present(name.split(" (")[0]) or present(name)))
    # cohort-level statement + data-availability note
    rows.append(dict(layer="Human SAH-CSF lipidomics (data availability)",
                     cohort=COHORT, analyte="__per_sample_abundance_matrix__",
                     lipid_class="-", direction="-",
                     significance="NOT downloadable",
                     fold_range="", ferroptosis_note="MacAulayLab GitHub repo contains CODE ONLY; "
                     "normalised per-sample .xlsx inputs not published; journal supplement is figures-only",
                     provenance="checked this run: repo file listing = scripts only; supp = PDF figures",
                     token_present_in_xml=False))
    return pd.DataFrame(rows)

# -----------------------------------------------------------------------------
# (B) REAL numeric ferroptosis-execution lipids (Round-10, openly downloadable)
# -----------------------------------------------------------------------------
def load_mouse_anchor():
    if not os.path.exists(CSV10):
        return pd.DataFrame()
    df = pd.read_csv(CSV10)
    df = df[df["status"].astype(str).str.contains("measured", na=False)].copy()
    # keep the execution-relevant classes: oxidised PE and free arachidonate
    sel = df[df["marker_class"].str.contains("ox-PE|free arachidonate", case=False, na=False)]
    return sel[["analyte","dataset","marker_class","contrast","tissue","fold_change",
                "log2FC","fold_change_CI_lo","fold_change_CI_hi","welch_p","cliffs_delta"]]

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
def make_figure(hum_fdr, hum_uncorr, mouse, path_png, path_pdf):
    fig = plt.figure(figsize=(10.4, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.02, 1.12, 0.96], wspace=0.56, left=0.145, right=0.975, top=0.80, bottom=0.20)

    # --- Panel A: human CSF verbatim differential lipids (class-level direction) ---
    axA = fig.add_subplot(gs[0, 0])
    cats = ["Arachidonoyl\namide  (FDR)", "Amide C18\n(FDR)",
            "PC 30:0/31:0/\n32:0  (FDR)", "Ether (plasmenyl)\nPE 43:6/36:1 (FDR)",
            "Phosphatidic\nacids  (nominal)", "Fatty acids /\nLPC / MAG (nom.)",
            "Some amides\ndown"]
    # x position encodes direction; log2FC shown on secondary for the nominal 5-16x band
    ypos = np.arange(len(cats))[::-1]
    for i,(yy,cat) in enumerate(zip(ypos, cats)):
        if "FDR" in cat:
            axA.plot([1, 3.2], [yy, yy], color=OI["verm"], lw=3, solid_capstyle="round")
            axA.scatter([3.2], [yy], color=OI["verm"], s=48, zorder=5, marker=">")
        elif "nominal" in cat or "nom." in cat:
            axA.add_patch(Rectangle((2.5, yy-0.16), 3.0, 0.32, color=OI["orange"], alpha=0.7))
            axA.text(5.6, yy, "5-16x", va="center", fontsize=7, color=OI["grey"])
        elif "down" in cat:
            axA.plot([-2.0, -0.3], [yy, yy], color=OI["blue"], lw=3, solid_capstyle="round")
            axA.scatter([-2.0], [yy], color=OI["blue"], s=48, zorder=5, marker="<")
    axA.axvline(0, color="#999999", lw=0.8)
    axA.set_yticks(ypos); axA.set_yticklabels(cats, fontsize=7)
    axA.set_xlim(-3.2, 6.6); axA.set_xlabel("SAH vs control direction", fontsize=8)
    axA.set_title("Human SAH-CSF lipidome\n(authors' reported findings — VERBATIM)", fontsize=9, fontweight="bold")

    # --- Panel B: mouse execution lipids real log2FC with CI whiskers ---
    axB = fig.add_subplot(gs[0, 1])
    if len(mouse):
        def shortlab(a):
            mm = re.match(r"(.*?)(?:T(\d)_(\d)I)?$", str(a))
            base, g1, g2 = mm.group(1), mm.group(2), mm.group(3)
            return f"{base} #{g1}{g2}" if (g1 and g2) else base
        m = mouse.sort_values("log2FC", ascending=True).reset_index(drop=True)
        y = np.arange(len(m))
        fc = m["log2FC"].values
        lo = np.log2(m["fold_change_CI_lo"].values)
        hi = np.log2(m["fold_change_CI_hi"].values)
        col = [OI["verm"] if "ox-PE" in c else OI["orange"] for c in m["marker_class"]]
        for yy,(v,l,h,cl) in enumerate(zip(fc,lo,hi,col)):
            axB.plot([l,h],[yy,yy], color=cl, lw=1.4)
            axB.scatter([v],[yy], color=cl, s=42, zorder=5, edgecolor="black", linewidth=0.3)
        axB.set_yticks(y); axB.set_yticklabels(
            [shortlab(a) for a in m["analyte"]], fontsize=6.2)
        axB.set_xlabel("log2 fold-change (ferroptosis vs control)", fontsize=8)
        axB.axvline(0, color="#999999", lw=0.8)
        axB.set_title("Ferroptosis EXECUTION lipids\n(real LC-MS/MS, mouse models)", fontsize=9, fontweight="bold")

    # --- Panel C: evidence map for the ferroptosis-lipid index ---
    axC = fig.add_subplot(gs[0, 2]); axC.axis("off")
    lines = [
        ("Ferroptosis-lipid evidence chain", True),
        ("", False),
        ("RNA layer  :  PUFA-PL machinery (ACSL4) +", False),
        ("              ferritin/ferritin-loss program", False),
        ("              in HDFI axis  (Fig 18/24/36)", False),
        ("", False),
        ("Execution   :  ox-PE 17-28x, free AA 2.44x,", False),
        ("  lipids     PE-20:4 substrates consumed", True),
        ("              (mouse ferroptosis models)", False),
        ("", False),
        ("Human CSF   :  ether/plasmenyl-PE +", False),
        ("  direction   arachidonoyl + PUFA-containing", False),
        ("              species UP in SAH (verbatim)", False),
        ("", False),
        ("NOT computable:", True),
        ("  a single per-patient human", False),
        ("  ferroptosis-lipid index, because", False),
        ("  the human CSF per-sample matrix", False),
        ("  is not openly deposited", False),
    ]
    ytop = 0.98; step = 0.045
    for i,(tx,bold) in enumerate(lines):
        axC.text(0.02, ytop-i*step, tx, fontsize=7.4,
                 fontweight=("bold" if bold else "normal"),
                 color=(OI["verm"] if bold else "#222222"))
    axC.set_title("Integration & limits", fontsize=9, fontweight="bold")

    fig.suptitle("Ferroptosis lipid layer of the HDFI axis — human SAH-CSF lipidomics (verbatim) + "
                 "real execution-lipid anchor",
                 fontsize=11, fontweight="bold", y=0.985)
    note = ("Human CSF per-sample abundance is NOT openly downloadable (MacAulayLab GitHub = code only; journal supplement is figures-only),\n"
            "so panel A reproduces the authors' published differential-lipid findings VERBATIM (13 SAH vs 11 control) and is never re-derived.\n"
            "Panel B numbers are real but from mouse cell-culture ferroptosis (valinomycin MEFs / RSL3 MDSCs), not human SAH; no paired within-patient HDFI-lipid correlation is claimed.")
    fig.text(0.5, 0.012, note, ha="center", va="bottom", fontsize=6.6, color=OI["grey"])
    fig.savefig(path_png, dpi=300)
    fig.savefig(path_pdf)
    plt.close(fig)

# -----------------------------------------------------------------------------
def main():
    hum = verbatim_human_findings()
    mouse = load_mouse_anchor()
    mouse = mouse.assign(layer="Mouse ferroptosis EXECUTION lipids (Round 10, real numeric anchor)",
                         significance="Welch p on log2FC")

    out = pd.concat([hum, mouse], axis=0, ignore_index=True)
    # reorder columns present
    cols = ["layer","cohort","analyte","lipid_class","direction","significance",
            "fold_range","ferroptosis_note","provenance","token_present_in_xml",
            "dataset","contrast","tissue","fold_change","log2FC","fold_change_CI_lo",
            "fold_change_CI_hi","welch_p","cliffs_delta"]
    cols = [c for c in cols if c in out.columns]
    out = out[cols]
    csv_path = os.path.join(TAB, "sah_csf_lipidomics.csv")
    out.to_csv(csv_path, index=False)

    # summary
    fdr_n = int((hum["significance"].str.contains("FDR-corrected")).sum())
    unc_n = int((hum["significance"].str.contains("nominal")).sum())
    oxpe = mouse[mouse["marker_class"].str.contains("ox-PE", na=False)] if len(mouse) else pd.DataFrame()
    summary = dict(
        human_cohort=COHORT,
        human_fdr_significant_up_classes=fdr_n,
        human_nominal_only_up_clusters=unc_n,
        human_per_sample_matrix_downloadable=False,
        human_provenance="VERBATIM published differential-lipid findings; per-sample abundance NOT open (code-only repo)",
        execution_anchor="Round-10 mouse ferroptosis LC-MS/MS (ST003337 valinomycin MEFs, ST004928 RSL3 MDSCs)",
        oxpe_n=int(len(oxpe)),
        oxpe_fold_range=(None if oxpe.empty else [float(oxpe["fold_change"].min()), float(oxpe["fold_change"].max())]),
        oxpe_best_p=(None if oxpe.empty else float(oxpe["welch_p"].min())),
        honest_limits=["4-HNE/MDA/GSH absent in panels and not imputed",
                       "execution numbers are mouse cell-culture, not human SAH",
                       "no within-patient HDFI <-> lipid correlation possible"],
        conclusion=("Ferroptosis-execution lipids are quantitatively supported by REAL public LC-MS/MS; "
                    "the human SAH-CSF lipidome directionally corroborates elevated PUFA-/ether-PE species but is "
                    "reported VERBATIM (its per-sample matrix is not openly downloadable). A single per-patient human "
                    "ferroptosis-lipid index therefore remains an author-must-produce measurement."),
    )
    with open(os.path.join(ROOT, "results", "sah_csf_lipidomics_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    make_figure(hum[hum["significance"].str.contains("FDR")],
                hum[hum["significance"].str.contains("nominal")], mouse,
                os.path.join(FIG, "Fig38_sah_csf_lipidomics.png"),
                os.path.join(FIG, "Fig38_sah_csf_lipidomics.pdf"))
    print("wrote", csv_path)
    print(json.dumps(summary, indent=2, default=str)[:800])
    print("Fig38 png/pdf written.")

if __name__ == "__main__":
    main()
