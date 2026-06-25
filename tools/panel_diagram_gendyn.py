#!/usr/bin/env python3
"""GENDYN 12 HP panel layout preview + footprint check (top polygon display).

Mirrors the control coordinates planned for src/GENDYN.cpp GENDYWidget after the
widen-to-12HP + display rework, in the shared Axon/Soma/Haptik visual style.
Flags overlaps / out-of-bounds, writes a PNG:

    python3 tools/panel_diagram.py     # writes gendyn_panel.png, prints check
"""
import matplotlib, math, os, random
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

W, H = 60.96, 128.5
D_KNOB, D_TRIM, D_JACK, SCREW = 10.0, 7.0, 8.4, 5.0
fig, ax = plt.subplots(figsize=(W/10, H/10), dpi=130)
ax.add_patch(Rectangle((0, 0), W, H, facecolor="#16162a", edgecolor="black"))

overlaps, placed = [], []
def circ(name, cx, cy, d, color):
    r = d/2
    if cx-r < 0 or cx+r > W or cy-r < 0 or cy+r > H: overlaps.append(f"OOB {name}")
    for n2, x2, y2, r2 in placed:
        if ((cx-x2)**2+(cy-y2)**2)**0.5 < (r+r2)-0.3: overlaps.append(f"OVERLAP {name}<->{n2}")
    placed.append((name, cx, cy, r))
    ax.add_patch(Circle((cx, cy), r, facecolor=color, edgecolor="white", lw=0.5, alpha=0.9))
    ax.text(cx, cy, name, ha="center", va="center", fontsize=3.2, color="white")
def sw(name, cx, cy):
    placed.append((name, cx, cy, 2.5))
    ax.add_patch(FancyBboxPatch((cx-2.5, cy-4.5), 5, 9, boxstyle="round,pad=0.2",
                                facecolor="#444455", edgecolor="white", lw=0.5))
    ax.text(cx, cy, name, ha="center", va="center", fontsize=3.0, color="white", rotation=90)

for x, y in [(1, 1), (54.96, 1), (1, 122), (54.96, 122)]:
    ax.add_patch(Circle((x+SCREW/2, y+SCREW/2), SCREW/2, facecolor="#888888", edgecolor="black", lw=0.4))
    placed.append(("screw", x+SCREW/2, y+SCREW/2, SCREW/2))

# ── morphing-polygon display across the top (shared screen style) ──
dx, dy, dw, dh = 5.5, 8, 50, 36
ax.add_patch(FancyBboxPatch((dx, dy), dw, dh, boxstyle="round,pad=0.4",
                            facecolor="#070712", edgecolor="#2b2b4d", lw=1.2))
random.seed(3)
N = 13
xs = [dx + dw*(i+0.5)/N for i in range(N)]
ys = [dy + dh/2 + (random.random()*2-1)*dh*0.32 for _ in range(N)]
ax.plot(xs+[xs[0]], ys+[ys[0]], color="#55e0a0", lw=1.2)        # GENDYN accent: green
for x, y in zip(xs, ys): ax.add_patch(Circle((x, y), 0.5, facecolor="#c0ffd8", edgecolor="none"))
ax.text(dx+3, dy+4, "GENDYN", ha="left", va="center", fontsize=6, color="#9affb8", weight="bold")

# ── top control row (non-CV): N | FREQ | LOCK | DIST | PERSIST ──
circ("N", 8.0, 54, D_KNOB, "#333344")
circ("FREQ", 20.5, 54, D_KNOB, "#333344")
sw("LOCK", 30.5, 54)
circ("DIST", 41.0, 54, D_KNOB, "#333344")
circ("PERS", 53.0, 54, D_KNOB, "#333344")

# ── 4 CV channel strips: knob (74) / attenuverter (86) / jack (96) ──
cols = [(9.0, "S.AMP"), (24.32, "S.DUR"), (39.64, "B.AMP"), (54.96, "B.WID")]
for x, nm in cols:
    circ(nm, x, 74, D_KNOB, "#333344")
    circ(nm+".a", x, 86, D_TRIM, "#555533")
    circ(nm+".cv", x, 96, D_JACK, "#224444")

# ── outputs ──
circ("OUT", 15.0, 112, D_JACK, "#442222")
circ("TRIG", 30.5, 112, D_JACK, "#442222")
circ("FREQ", 46.0, 112, D_JACK, "#442222")

ax.set_xlim(-2, W+2); ax.set_ylim(H+2, -2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("GENDYN 12HP — morphing-polygon display", fontsize=8)
out = os.path.join(os.getcwd(), "gendyn_panel.png")
plt.tight_layout(); plt.savefig(out, dpi=140, bbox_inches="tight")
print("wrote", out)
print("Layout check:", "  ".join(overlaps) if overlaps else "no overlaps / out-of-bounds")
