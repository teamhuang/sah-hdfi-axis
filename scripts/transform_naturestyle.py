"""transform_naturestyle.py -- one-time injector.

Adds `import naturestyle as _ns; _ns.savefig_patch()` just after a figure script's
matplotlib import so every fig.savefig() in that script is auto-converted:
in-figure notes stripped, multi-panel composites split into independent <name>_a/_b
panels, single panels rendered at the 5 cm footprint. Idempotent.
"""
import os, re, glob

HERE = os.path.dirname(os.path.abspath(__file__))
INJECT = ("\n# --- naturestyle: strip in-figure notes; split composites to independent "
          "panels; 5 cm footprint ---\nimport naturestyle as _ns\n")

# basenames to keep as one whole image (schematics / non-panel layouts)
WHOLE = {
    "20_graphical_abstract.py": ["Fig0_graphical_abstract"],
    "15_newfigures.py": ["Fig10_axis_schematic"],
    "48_sah_csf_lipidomics.py": ["Fig38_sah_csf_lipidomics"],
}

changed = []
for path in sorted(glob.glob(os.path.join(HERE, "*.py"))):
    name = os.path.basename(path)
    if name in ("naturestyle.py", "transform_naturestyle.py", "build_supp_tables.py"):
        continue
    src = open(path, encoding="utf-8").read()
    if "savefig" not in src:
        continue
    if "import naturestyle" in src:
        continue
    # injection point: after the first matplotlib-related import line, else after top imports
    lines = src.splitlines(keepends=True)
    anchor = None
    for i, ln in enumerate(lines[:60]):
        if re.search(r"^\s*(import matplotlib|matplotlib\.use|from matplotlib import|import matplotlib\.pyplot)", ln):
            anchor = i
            break
    if anchor is None:
        for i, ln in enumerate(lines[:40]):
            if ln.startswith(("import ", "from ")) and "os.path.dirname" not in ln:
                anchor = i
    if anchor is None:
        anchor = 0
    block = INJECT
    if name in WHOLE:
        block += "_ns.mark_whole(" + ", ".join(f'"{w}"' for w in WHOLE[name]) + ")\n"
    block += "_ns.savefig_patch()\n"
    lines.insert(anchor + 1, block)
    open(path, "w", encoding="utf-8").write("".join(lines))
    changed.append(name)

print("patched", len(changed), "scripts:")
for c in changed:
    print("  ", c)
