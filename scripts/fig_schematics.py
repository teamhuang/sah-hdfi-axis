"""fig_schematics.py -- clean redraw of the two framing schematics.

Fig10_axis_schematic.png  (Main Figure 1: overview of the axis + study design)
Fig0_graphical_abstract.png (Main Figure 8: integrated model)

Minimal box-and-arrow flows, short labels only, no on-image notes/titles (explanation
lives in the caption), sized to a 6.5 cm landscape footprint.
"""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import naturestyle as ns

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
C = ns.OKABE


def box(ax, x, y, w, h, txt, fc, fs=6.5, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.02",
                                fc=fc, ec="none"))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
            fontsize=fs, color=tc, weight="bold", linespacing=1.1)


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10,
                                 lw=1.0, color="#444444", shrinkA=1, shrinkB=1))


def canvas():
    fig = plt.figure(figsize=(6.5 / 2.54, 5.0 / 2.54))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    return fig, ax


def fig10():
    fig, ax = canvas()
    w, h = 0.18, 0.20
    xs = [0.012, 0.209, 0.406, 0.603, 0.80]
    r1, r2 = 0.66, 0.16
    fs = 5.4
    # Row 1: heme release -> clearance -> catabolism -> labile iron / ROS (left -> right)
    box(ax, xs[0], r1, w, h, "SAH\nRBC lysis", C["vermil"], fs=fs)
    box(ax, xs[1], r1, w, h, "free\nHb", C["orange"], fs=fs)
    box(ax, xs[2], r1, w, h, "HP\nCD163/MRC1", C["green"], fs=fs)
    box(ax, xs[3], r1, w, h, "HMOX1\nheme", C["blue"], fs=fs)
    box(ax, xs[4], r1, w, h, "Fe2+\nROS", C["purple"], fs=fs)
    for i in range(4):
        arrow(ax, xs[i] + w, r1 + h / 2, xs[i + 1], r1 + h / 2)
    # Connector: straight drop from Fe2+/ROS down into the buffering arm (serpentine)
    arrow(ax, xs[4] + w / 2, r1, xs[4] + w / 2, r2 + h)
    # Row 2: buffer overwhelmed -> failure state -> ferroptosis -> injury (right -> left)
    box(ax, xs[4], r2, w, h, "ferritin\n+ GPX4", C["sky"], tc="#123", fs=fs)
    box(ax, xs[3], r2, w, h, "HDFI\nstate", C["vermil"], fs=fs)
    box(ax, xs[2], r2, w, h, "neuronal\nferroptosis", C["blue"], fs=fs)
    box(ax, xs[1], r2, w, h, "secondary\ninjury", C["green"], fs=fs)
    for i in (4, 3, 2):
        arrow(ax, xs[i], r2 + h / 2, xs[i - 1] + w, r2 + h / 2)
    fig.savefig(os.path.join(FIG, "Fig10_axis_schematic.png"), dpi=300)
    fig.savefig(os.path.join(FIG, "Fig10_axis_schematic.pdf"))
    plt.close(fig)


def fig0():
    fig, ax = canvas()
    w, h = 0.18, 0.22
    xs = [0.012, 0.209, 0.406, 0.603, 0.80]
    y = 0.39
    fs = 5.4
    box(ax, xs[0], y, w, h, "SAH\nred-cell\nheme stress", C["vermil"], fs=fs)
    box(ax, xs[1], y, w, h, "myeloid\nheme\nhandling", C["green"], fs=fs)
    box(ax, xs[2], y, w, h, "HDFI\nstate", C["orange"], fs=fs)
    box(ax, xs[3], y, w, h, "neuronal\nferroptosis", C["blue"], fs=fs)
    box(ax, xs[4], y, w, h, "secondary\nbrain injury", C["purple"], fs=fs)
    for i in range(4):
        arrow(ax, xs[i] + w, y + h / 2, xs[i + 1], y + h / 2)
    fig.savefig(os.path.join(FIG, "Fig0_graphical_abstract.png"), dpi=300)
    fig.savefig(os.path.join(FIG, "Fig0_graphical_abstract.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig10(); fig0()
    print("redrew Fig10 + Fig0 schematics")
