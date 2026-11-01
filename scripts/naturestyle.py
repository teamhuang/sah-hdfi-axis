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


def savefig_patch():
    """Monkeypatch Figure.savefig so legacy scripts are converted by one import.

    Any fig.savefig(...) then strips notes and renders at the 5 cm footprint with
    width chosen by aspect (square -> 5 cm, landscape -> 6.5 cm)."""
    _orig = matplotlib.figure.Figure.savefig

    def _patched(self, fname, *a, **k):
        _strip_notes(self)
        try:
            self.canvas.draw()
            bb = self.get_tightbbox(self.canvas.get_renderer())
            ratio = (bb.width / bb.height) if bb and bb.height else 1.0
        except Exception:
            ratio = 1.0
        h_cm = 5.0
        w_cm = 5.0 if ratio < 1.15 else 6.5
        self.set_size_inches(w_cm * _cm, h_cm * _cm, forward=True)
        try:
            self.tight_layout(pad=0.4)
        except Exception:
            pass
        k.setdefault("dpi", 300)
        k.setdefault("facecolor", "white")
        k["bbox_inches"] = "tight"
        return _orig(self, fname, *a, **k)

    matplotlib.figure.Figure.savefig = _patched
