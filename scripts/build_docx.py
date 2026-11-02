"""build_docx.py -- render a submission markdown to a Brain/Stroke-style .docx.

Rules honoured:
  * Bold + black ONLY for: document title, section/subsection headings, table
    captions, figure captions. Body text is never bold.
  * Tables rendered as three-line (top + header + bottom rules, no verticals).
  * Bracketed numbered citations [1], [2,3], [25-27] rendered as superscript.
  * Figure images inserted at the END (under their legend), sized to a 5 cm tall /
    <=6.5 cm wide footprint; NO panel letters drawn on the image (letters are in the
    caption). Multi-panel source figures are inserted as their split panels.

Usage: python build_docx.py <in.md> <out.docx> <main|supp>
"""
import os, re, sys, glob
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

FIGDIR = "figures"
BODY_FONT = "Times New Roman"
HEAD_FONT = "Arial"

def _set_font(run, name=BODY_FONT, size=11, bold=False, italic=False, color=(0, 0, 0)):
    run.font.name = name
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    run.font.color.rgb = RGBColor(*color)
    rpr = run._element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts"); rpr.append(rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rf.set(qn(a), name)

CITE = re.compile(r"(\[[0-9]+(?:[\s,\-–][0-9]+)*\])")

def add_body(doc, text, size=11, italic=False, align=None, first_bold_caption=False):
    p = doc.add_paragraph()
    if align:
        p.alignment = align
    # caption form: "**Figure 2.** rest" / "**Table S1.** rest" -> leading bold, rest normal
    m = re.match(r"^\*\*(.+?)\.\*\*\s*(.*)$", text)
    if m:
        r = p.add_run(m.group(1) + "."); _set_font(r, HEAD_FONT, size, bold=True)
        rest = m.group(2)
        _add_cited(p, rest, size)
        return p
    # strip any other inline bold -> body must not be bold
    text = text.replace("**", "").replace("__", "")
    _add_cited(p, text, size, italic=italic)
    return p

def _add_cited(p, text, size=11, italic=False):
    if not text:
        return
    for part in CITE.split(text):
        if not part:
            continue
        if CITE.fullmatch(part):
            r = p.add_run(part[1:-1])          # drop brackets
            _set_font(r, BODY_FONT, size, italic=italic)
            r.font.superscript = True
        else:
            r = p.add_run(part)
            _set_font(r, BODY_FONT, size, italic=italic)

def add_heading(doc, text, level):
    text = text.replace("**", "").strip()
    p = doc.add_paragraph()
    sizes = {0: 15, 1: 13, 2: 12, 3: 11}
    r = p.add_run(text)
    _set_font(r, HEAD_FONT, sizes.get(level, 11), bold=True)
    p.space_before = Pt(8); p.space_after = Pt(4)
    return p

def _border(el, sz=8, val="single"):
    b = OxmlElement(el)
    b.set(qn("w:val"), val); b.set(qn("w:sz"), str(sz)); b.set(qn("w:color"), "000000")
    return b

def three_line(table, n_header=1):
    tbl = table._tbl
    tblPr = tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge, sz, on in (("top", 12, True), ("left", 0, False), ("bottom", 12, True),
                         ("right", 0, False), ("insideH", 0, False), ("insideV", 0, False)):
        e = _border(f"w:{edge}", sz, "single" if on else "none")
        borders.append(e)
    tblPr.append(borders)
    # header bottom rule via first row
    for i, row in enumerate(table.rows):
        if i < n_header:
            for cell in row.cells:
                tcPr = cell._tc.get_or_add_tcPr()
                b = OxmlElement("w:tcBorders")
                b.append(_border("w:bottom", 8, "single"))
                tcPr.append(b)

def render_md_table(doc, rows):
    header = [c.strip() for c in rows[0].strip().strip("|").split("|")]
    body = []
    for r in rows[2:]:
        body.append([c.strip() for c in r.strip().strip("|").split("|")])
    t = doc.add_table(rows=1 + len(body), cols=len(header))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for j, h in enumerate(header):
        cell = t.rows[0].cells[j]
        cell.text = ""
        _set_font(cell.paragraphs[0].add_run(h.replace("**", "")), HEAD_FONT, 9, bold=True)
    for i, row in enumerate(body):
        for j in range(len(header)):
            val = row[j] if j < len(row) else ""
            cell = t.rows[i + 1].cells[j]
            cell.text = ""
            _add_cited(cell.paragraphs[0], val.replace("**", ""), size=9)
    three_line(t)
    return t

def insert_figure_panels(doc, base_names):
    """Insert every panel image for the given source-figure basenames, sized to footprint."""
    files = []
    for b in base_names:
        hits = sorted(glob.glob(os.path.join(FIGDIR, b + "*.png")))
        # prefer split panels if present, else the whole figure
        panels = [h for h in hits if re.search(r"_[a-z]\.png$", h)]
        files.extend(panels if panels else [h for h in hits if not re.search(r"_[a-z]\.png$", h)])
    for f in files:
        try:
            from PIL import Image
            w, h = Image.open(f).size
        except Exception:
            w, h = (1500, 1200)
        aspect = w / h if h else 1.0
        width_cm = 5.0 if aspect < 1.2 else 6.5
        height_cm = width_cm / aspect
        if height_cm > 5.0:
            height_cm = 5.0; width_cm = height_cm * aspect
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(f, width=Cm(width_cm), height=Cm(height_cm))
    return len(files)

def main():
    in_md, out_docx, kind = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else "main")
    # figure manifest: main figure number -> source basenames (order = panel order)
    MANIFEST = {
        1: ["Fig10_axis_schematic"],
        2: ["Fig1_volcano_GSE36791", "Fig3_heatmap_candidates", "Fig4_enrichment",
            "Fig14_bulk_deconvolution", "Fig19_composition_adjusted"],
        3: ["Fig8_scrna_summary", "Fig12_scrna_deep", "Fig13_scrna_temporal", "Fig20_cellstate"],
        4: ["Fig16_cellchat_style", "Fig25_nichenet_style", "Fig34_cite_seq_cd163"],
        5: ["Fig15_regulon_bach1_nrf2", "Fig18_hdfi_centrepiece", "Fig36_hdfi_latent_axis",
            "Fig41_hdfi_robustness"],
        6: ["Fig21_regulatory_motif", "Fig26_nrf2_chip_binding", "Fig32_scatac_hmox1"],
        7: ["Fig33_lipidomics_ferroptosis", "Fig37_proteomics_hdfi",
            "Fig40_csf_hdfi_neuronal_injury", "Fig23_hdfi_external"],
        8: ["Fig0_graphical_abstract"],
    }
    SUPP_MANIFEST = {
        1: ["Fig7_wgcna"], 2: ["Fig5_ml_ranking"], 3: ["Fig6_score_and_validation"],
        4: ["Fig9_xgb"], 5: ["umap_celltype", "umap_heme", "dotplot__markers"],
        6: ["Fig39_mouse_sah_atlas"], 7: ["Fig11_scrna_trajectory"],
        8: ["Fig17_ferroptosis_gsva"], 9: ["Fig29_lipidomics_stroke"],
        10: ["Fig38_sah_csf_lipidomics"], 11: ["Fig22_mr_forest"], 12: ["Fig24_mr_panel"],
        13: ["Fig27_spatial_coloc"], 14: ["Fig28_drugrepurposing_targets"],
        15: ["Fig35_drug_reversal_lincs"],
    }
    man = MANIFEST if kind == "main" else SUPP_MANIFEST

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = BODY_FONT; st.font.size = Pt(11)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Cm(2.2); s.left_margin = s.right_margin = Cm(2.4)

    lines = open(in_md, encoding="utf-8").read().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if not s or s in ("---", "***"):
            i += 1; continue
        if s.startswith("|") and i + 1 < len(lines) and set(lines[i + 1].replace("|", "").strip()) <= set("-: "):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(lines[i]); i += 1
            render_md_table(doc, rows)
            doc.add_paragraph()
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            level = len(m.group(1))
            add_heading(doc, m.group(2), level - 1)
            i += 1; continue
        # figure legend caption -> insert panels right after the legend paragraph
        fc = re.match(r"^\*\*(Figure\s+)(S?)(\d+)\b", s)
        if fc:
            add_body(doc, s, size=10)
            is_supp, num = bool(fc.group(2)), int(fc.group(3))
            want_supp = (kind == "supp")
            if is_supp == want_supp and num in man:
                insert_figure_panels(doc, man[num])
            i += 1; continue
        # table caption
        tc = re.match(r"^\*\*(Table\s+S?)(\d+)\.?\*\*", s)
        if tc:
            add_body(doc, s, size=10)
            i += 1; continue
        add_body(doc, s, size=11)
        i += 1

    doc.save(out_docx)
    print("wrote", out_docx)

if __name__ == "__main__":
    main()
