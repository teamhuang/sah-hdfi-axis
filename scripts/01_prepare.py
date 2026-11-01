import os, gzip, csv, re
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

# ---- 1. probe -> gene symbol map from GPL10558 annotation ----
def load_probe_symbol():
    m = {}
    with gzip.open(os.path.join(DATA, "GPL10558.annot.gz"), "rt", encoding="utf-8", errors="replace") as f:
        rdr = csv.reader(f, delimiter="\t")
        header = None
        for row in rdr:
            if not row:
                continue
            if row[0] == "ID":
                header = row
                continue
            if header and len(row) == len(header) and row[0].startswith("ILMN_"):
                pid = row[0]
                sym = row[2].strip()
                gid = row[3].strip()
                if sym and sym not in ("---", "NO_ID", "NA"):
                    sym = sym.split("///")[0].split()[0]
                    if re.fullmatch(r"(HGN\w+|LOC\d+|MIR\d+\w*|[A-Z0-9.\-]+)", sym) and not sym.startswith("---"):
                        m[pid] = sym  # keep named gene symbols only (drop ENTREZID-only probes)
    return m

# ---- 2. parse a GEO series matrix ----
def parse_matrix(fn):
    with gzip.open(os.path.join(DATA, fn), "rt", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    chars = []
    title = geoid = source = None
    begin = None
    for i, l in enumerate(lines):
        ls = l.rstrip("\n")
        if ls.startswith("!Sample_geo_accession"):
            geoid = [x.strip('"') for x in ls.split("\t")[1:]]
        elif ls.startswith("!Sample_title"):
            title = [x.strip('"') for x in ls.split("\t")[1:]]
        elif ls.startswith("!Sample_source_name_ch1"):
            source = [x.strip('"') for x in ls.split("\t")[1:]]
        elif ls.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in ls.split("\t")[1:]]
            chars.append(vals)
        elif ls.startswith("!series_matrix_table_begin"):
            begin = i + 1
            break
    # data table
    rows = []
    gsm_cols = None
    for ls in lines[begin:]:
        if ls.startswith("!series_matrix_table_end"):
            break
        parts = next(csv.reader([ls], delimiter="\t"))
        if parts[0] == "ID_REF":
            gsm_cols = [x.strip('"') for x in parts[1:]]
            continue
        rows.append(parts)
    expr = pd.DataFrame(
        {g: pd.to_numeric([r[i + 1] for r in rows], errors="coerce") for i, g in enumerate(gsm_cols)},
        index=[r[0] for r in rows],
    )
    expr.index.name = "probe"
    meta = pd.DataFrame({"gsm": geoid, "title": title, "source": source})
    # extract group + sex from characteristics
    def pick(key):
        for vals in chars:
            if key.lower() in vals[0].lower():
                return [v.split(":", 1)[1].strip() for v in vals]
        return None
    grp = pick("group") or pick("disease status")
    sex = pick("Sex")
    meta["group_raw"] = grp
    meta["sex"] = sex
    return expr, meta

def collapse(expr, pmap, min_sym=1):
    expr = expr.copy()
    expr["symbol"] = [pmap.get(p, None) for p in expr.index]
    expr = expr.dropna(subset=["symbol"])
    # keep probes with finite values, collapse by mean across probes per symbol (linear),
    # then log2 later in R. Use mean of available values.
    g = expr.groupby("symbol").mean()
    return g

# ---- 3. harmonise group labels ----
def label_group(s):
    s2 = s.lower()
    if "ruptured" in s2 or "aneurysm" in s2 or s2 == "case":
        return "SAH"
    if "control" in s2:
        return "Control"
    return None

if __name__ == "__main__":
    pmap = load_probe_symbol()
    print(f"probe->symbol map: {len(pmap)} probes")
    for fn, tag in [("GSE36791_series_matrix.txt.gz", "GSE36791"),
                    ("GSE73378-GPL10558_series_matrix.txt.gz", "GSE73378")]:
        expr, meta = parse_matrix(fn)
        print(f"\n{tag}: {expr.shape[0]} probes x {expr.shape[1]} samples; expr NA frac={expr.isna().mean().mean():.4f}")
        sym = collapse(expr, pmap)
        sym.to_csv(os.path.join(DATA, f"expr_linear_{tag}.tsv"), sep="\t")
        meta["gsm"] = meta["gsm"].str.strip()
        meta = meta[meta.gsm.isin(sym.columns)]
        meta["group"] = meta["group_raw"].map(label_group)
        meta.to_csv(os.path.join(DATA, f"pheno_{tag}.tsv"), sep="\t", index=False)
        print(f"{tag}: collapsed to {sym.shape[0]} gene symbols x {sym.shape[1]} samples")
        print(f"{tag}: {sym.columns.difference(meta.gsm).tolist()} unmatched cols")
        print(meta.groupby("group").size())
        print(meta.groupby(["group", "sex"]).size())
