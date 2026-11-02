"""renumber_refs.py -- reorder the reference list by first appearance in the body
and rewrite every in-text citation to the new numbering (Brain/Stroke: numbered,
superscript, cited in order, none uncited)."""
import re, sys

path = "manuscript/MAIN_submission.md"
text = open(path, encoding="utf-8").read()
cut = text.find("## References")
body, refblock = text[:cut], text[cut:]

# --- parse reference list ---
refpat = re.compile(r"^(\d+)\.\s+(.*)$", re.M)
refs = {}
for m in refpat.finditer(refblock):
    refs[int(m.group(1))] = m.group(2).strip()

# --- expand a citation group into ints ---
def expand(grp):
    out = []
    for part in re.split(r",", grp):
        part = part.strip()
        rng = re.match(r"^(\d+)\s*[\u2013-]\s*(\d+)$", part)
        if rng:
            out += list(range(int(rng.group(1)), int(rng.group(2)) + 1))
        elif part.isdigit():
            out.append(int(part))
    return out

# --- assign new numbers by first appearance ---
order = []
def cite_repl(m):
    nums = expand(m.group(1))
    for n in nums:
        if n not in order:
            order.append(n)
    new = sorted({order.index(n) + 1 for n in nums})
    return "[" + ",".join(str(x) for x in new) + "]"

body2 = re.sub(r"\[([0-9][0-9,\s\u2013-]*)\]", cite_repl, body)

# --- rebuild reference list in new order ---
lines = []
for i, old in enumerate(order, start=1):
    lines.append(f"{i}. {refs[old]}")
newrefs = "## References\n\n" + "\n".join(lines) + "\n"

open(path, "w", encoding="utf-8").write(body2 + newrefs)
print(f"renumbered {len(order)} references; dropped uncited:",
      sorted(set(refs) - set(order)))
