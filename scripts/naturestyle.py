"""naturestyle.py -- single shared publication-figure style.

Enforces, for every figure produced by the pipeline:
  * Arial, Okabe-Ito (colour-blind-safe) palette, light spines, no legend frame;
  * removal of in-figure "notes" (axes titles + suptitles) -- captions live in the
    manuscript, not on the image;
  * a fixed footprint: height 5 cm, width either 5 cm (square panels) or 6.5 cm
    (landscape panels), chosen from the panel aspect ratio.

Import it near the top of a figure script (right after `import matplotlib`) with:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import naturestyle as ns
    ns.here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

then call ns.save(fig, "FigX", width_cm=5 or 6.5) instead of fig.savefig.
ns.savefig_patch() can also monkeypatch plt.Figure.savefig for scripts that call
fig.savefig directly, so a single import converts an old script with no other edits.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.transforms import Bbox

_cm = 1.0 / 2.54
here = os.getcwd()
FIGDIR = os.path.join(here, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# Okabe-Ito colour-blind-safe qualitative palette
OKABE = {
    "blue": "#0072B2", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "vermil": "#D55E00", "purple": "#CC79A7", "black": "#000000",
    "grey": "#7F7F7F", "lightgrey": "#B9BEC4",
}
# sequential ramps (perceptually uniform, colourblind-friendly)
SEQ = "cividis"
DIV = "RdBu_r"

def apply():
    """Apply Nature-family rcParams (idempotent)."""
    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "axes.prop_cycle": plt.cycler(color=[
            OKABE["blue"], OKABE["vermil"], OKABE["green"], OKABE["orange"],
            OKABE["purple"], OKABE["sky"], OKABE["black"]]),
        "pdf.fonttype": 42,   # editable text in the vector PDF
        "ps.fonttype": 42,
    })

apply()

# colour accessors for scripts
BLUE, ORANGE, SKY, GREEN = OKABE["blue"], OKABE["orange"], OKABE["sky"], OKABE["green"]
YELLOW, VERMIL, PURPLE, BLACK = OKABE["yellow"], OKABE["vermil"], OKABE["purple"], OKABE["black"]
GREY, LIGHTGREY = OKABE["grey"], OKABE["lightgrey"]


def _strip_notes(fig):
    """Remove suptitles and per-axes titles (these are the on-image 'notes')."""
    try:
        if fig._suptitle is not None:
            fig.suptitle("")
    except Exception:
        pass
    for ax in fig.axes:
        try:
            if ax.get_title():
                ax.set_title("")
        except Exception:
            pass


def save(fig, name, width_cm=6.5, height_cm=5.0, strip_notes=True, outdir=None,
         tight=True, dpi=300):
    """Save a panel as PNG + editable PDF at the fixed Nature footprint.

    width_cm must be 5 or 6.5; height_cm defaults to 5. In-figure notes stripped.
    """
    if width_cm not in (5, 6.5):
        width_cm = 6.5 if width_cm > 5.5 else 5
    if strip_notes:
        _strip_notes(fig)
    outdir = outdir or FIGDIR
    fig.set_size_inches(width_cm * _cm, height_cm * _cm, forward=True)
    try:
        fig.tight_layout(pad=0.4)
    except Exception:
        pass
    base = os.path.join(outdir, name)
    kw = dict(dpi=dpi, facecolor="white")
    if tight:
        kw["bbox_inches"] = "tight"
    fig.savefig(base + ".png", **kw)
    fig.savefig(base + ".pdf", format="pdf", **kw)
    plt.close(fig)
    return base


_ORIG = matplotlib.figure.Figure.savefig

# basenames (no extension) that must stay a single whole image (schematics,
# graphical abstract); multi-axes figures not in this set are split into panels.
WHOLE = set()

# Split multi-panel source figures into independent panels (readable at 5 cm);
# colourbars are merged into their parent panel. Set False to keep whole composites.
DO_SPLIT = True


def mark_whole(*names):
    """Register figure basenames that should never be split into panels."""
    for n in names:
        WHOLE.add(os.path.splitext(os.path.basename(n))[0])


def _panel_label(i):
    return chr(ord("a") + i) if i < 26 else str(i + 1)


def _is_thin(ax):
    p = ax.get_position()
    return p.width < 0.035 or p.height < 0.035


def _is_colorbar(ax):
    """Robust colourbar test: matplotlib tags the colourbar axes with _colorbar_info.

    A width/height threshold is unreliable (a colourbar can be ~0.05 wide after
    tight_layout), so we key off the marker the Colorbar machinery sets and keep
    the thin-axes test only as a fallback for annotation strips.
    """
    return hasattr(ax, "_colorbar_info") or hasattr(ax, "_colorbar") or _is_thin(ax)


def _panels_with_colorbars(self):
    """Return [(data_axes, [its colourbar axes, ...])] in visual order.

    A colourbar is attached to its TRUE parent axes (recorded by matplotlib in
    ``_colorbar_info['parents']``), not to the geometrically nearest axes -- the
    latter mis-assigns colourbars in tight gridspecs and bleeds their tick labels
    into a neighbouring panel's crop.
    """
    data = [ax for ax in self.axes if not _is_colorbar(ax)]
    cbs = [ax for ax in self.axes if _is_colorbar(ax)]
    def key(ax):
        p = ax.get_position()
        return (-round(p.y0, 2), round(p.x0, 2))
    data = sorted(data, key=key)
    dataset = set(data)
    attach = {ax: [] for ax in data}
    def center(ax):
        p = ax.get_position(); return (p.x0 + p.width / 2, p.y0 + p.height / 2)
    for cb in cbs:
        parents = getattr(cb, "_colorbar_info", {}).get("parents", []) if hasattr(cb, "_colorbar_info") else []
        parent = next((p for p in parents if p in dataset), None)
        if parent is None:                       # fallback: nearest data axes
            cx, cy = center(cb)
            parent = min(data, key=lambda d: ((center(d)[0] - cx) ** 2 + (center(d)[1] - cy) ** 2) ** 0.5)
        attach[parent].append(cb)
    return [(ax, attach[ax]) for ax in data]


def _data_axes(self):
    return [ax for ax, _ in _panels_with_colorbars(self)]


def savefig_patch():
    """Monkeypatch Figure.savefig so legacy scripts are converted by one import.

    Any fig.savefig(...) then strips in-figure notes and either (a) splits a
    genuine multi-panel figure into independent <name>_a/_b/... PNG+PDF panels, or
    (b) renders a single-panel figure at the 5 cm footprint (width 5 cm square,
    6.5 cm landscape). Figures registered via mark_whole() are never split."""
    def _patched(self, fname, *a, **k):
        try:
            return _restyled(self, fname, k)
        except Exception as exc:                       # never break an analysis run
            print("naturestyle fallback:", str(exc)[:120])
            k.setdefault("dpi", 300)
            k.setdefault("facecolor", "white")
            return _ORIG(self, fname, *a, **k)

    matplotlib.figure.Figure.savefig = _patched


def _whole(self, fname, dpi):
    """Strip notes, render a single image (or registered-whole figure) at 5 cm footprint."""
    _strip_notes(self)
    self.set_dpi(150)
    try:
        bb = self.get_tightbbox(self.canvas.get_renderer())
        ratio = (bb.width / bb.height) if bb and bb.height else 1.0
    except Exception:
        ratio = 1.0
    w_cm = 5.0 if ratio < 1.2 else 6.5
    self.set_size_inches(w_cm * _cm, 5.0 * _cm, forward=False)
    try:
        self.tight_layout(pad=0.4)
    except Exception:
        pass
    return _ORIG(self, fname, dpi=dpi, facecolor="white", bbox_inches="tight")


def _split(self, fname, base, panels, dpi):
    """Save each content axes (with its colourbar) as an independent <base>_a/_b ... PNG+PDF."""
    self.set_dpi(150)
    rpix = self.canvas.get_renderer().dpi          # pixels per inch of the tight bbox
    outs = []
    seen_tb = []                                   # guard against duplicate crops
    panel_i = 0
    for i, (ax, cbs) in enumerate(panels):
        try:
            boxes = [ax.get_tightbbox(self.canvas.get_renderer())]
            for cb in cbs:
                boxes.append(cb.get_tightbbox(self.canvas.get_renderer()))
            tb = Bbox.union(boxes)
            # get_tightbbox returns a Bbox in PIXELS; savefig's bbox_inches wants INCHES
            tb = Bbox([[tb.x0 / rpix, tb.y0 / rpix], [tb.x1 / rpix, tb.y1 / rpix]])
        except Exception:
            tb = None
        if tb is None or tb.width < 0.5 or tb.height < 0.5:
            continue
        sig = (round(tb.x0, 2), round(tb.y0, 2), round(tb.x1, 2), round(tb.y1, 2))
        if sig in seen_tb:                          # same region already written -> skip
            continue
        seen_tb.append(sig)
        pbase = f"{base}_{_panel_label(panel_i)}"
        panel_i += 1
        _ORIG(self, pbase + ".png", bbox_inches=tb, dpi=dpi, facecolor="white")
        try:
            _ORIG(self, pbase + ".pdf", format="pdf", bbox_inches=tb, dpi=dpi, facecolor="white")
        except Exception:
            pass
        outs.append(pbase)
    return outs


def _restyled(self, fname, k):
    dpi = k.get("dpi", 300)
    _strip_notes(self)
    stem = os.path.splitext(os.path.basename(str(fname)))[0]
    panels = _panels_with_colorbars(self)
    base = os.path.splitext(str(fname))[0]
    if DO_SPLIT and len(panels) > 1 and stem not in WHOLE:
        outs = _split(self, fname, base, panels, dpi)
        if outs:
            return outs
    return _whole(self, fname, dpi)

