#!/usr/bin/env python3
"""Axon panel layout preview + footprint check (12 HP, top phase-portrait display).

Mirrors the control coordinates in src/Axon.cpp AxonWidget. Draws them at
realistic sizes, flags overlaps / out-of-bounds, and writes a preview PNG:

    python3 tools/panel_diagram.py     # writes axon_panel.png, prints check
"""
import matplotlib, math, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

W, H = 60.96, 128.5
D_KNOB, D_TRIM, D_JACK, SCREW = 10.0, 7.0, 8.4, 5.0
fig, ax = plt.subplots(figsize=(W/10, H/10), dpi=130)
ax.add_patch(Rectangle((0, 0), W, H, facecolor="#16162a", edgecolor="black"))

overlaps, placed = [], []
def circ(name, cx, cy, d, color):
    r = d / 2
    if cx - r < 0 or cx + r > W or cy - r < 0 or cy + r > H:
        overlaps.append(f"OOB {name}")
    for n2, x2, y2, r2 in placed:
        if ((cx - x2)**2 + (cy - y2)**2)**0.5 < (r + r2) - 0.3:
            overlaps.append(f"OVERLAP {name}<->{n2}")
    placed.append((name, cx, cy, r))
    ax.add_patch(Circle((cx, cy), r, facecolor=color, edgecolor="white", lw=0.5, alpha=0.9))
    ax.text(cx, cy, name, ha="center", va="center", fontsize=3.4, color="white")

for x, y in [(1, 1), (54.96, 1), (1, 122), (54.96, 122)]:
    ax.add_patch(Circle((x + SCREW/2, y + SCREW/2), SCREW/2, facecolor="#888888", edgecolor="black", lw=0.4))
    placed.append(("screw", x + SCREW/2, y + SCREW/2, SCREW/2))

# ── phase-portrait display across the top ──
dx, dy, dw, dh = 5.5, 8, 50, 34
ax.add_patch(FancyBboxPatch((dx, dy), dw, dh, boxstyle="round,pad=0.4",
                            facecolor="#070712", edgecolor="#2b2b4d", lw=1.2))
# illustrative limit-cycle orbit
cx, cy = dx + dw/2, dy + dh/2
pts = [(cx + 0.40*dw*math.cos(t), cy - 0.30*dh*math.sin(t))
       for t in [2*math.pi*i/80 for i in range(81)]]
ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#55c8ff", lw=1.2)
ax.text(dx + 3, dy + 4, "AXON", ha="left", va="center", fontsize=6, color="#9ab0ff", weight="bold")
ax.text(dx + dw - 3, dy + dh - 3, "v–w", ha="right", va="center", fontsize=4, color="#6080b0")

# ── controls below ──
# knob row (y=54): PITCH | CURRENT | EPS | SHAPE
circ("PITCH", 9.0, 54, D_KNOB, "#333344")
circ("CURR", 24.32, 54, D_KNOB, "#333344")
circ("EPS", 39.64, 54, D_KNOB, "#333344")
circ("SHAPE", 54.96, 54, D_KNOB, "#333344")
# CV strips under CURRENT and EPS: attenuverter (y=72) + input jack (y=84)
circ("I.a", 24.32, 72, D_TRIM, "#555533")
circ("e.a", 39.64, 72, D_TRIM, "#555533")
circ("I.cv", 24.32, 84, D_JACK, "#224444")
circ("e.cv", 39.64, 84, D_JACK, "#224444")
# I/O row (y=112): three inputs | three outputs
circ("VOCT", 6.5, 112, D_JACK, "#224444")
circ("TRIG", 15.5, 112, D_JACK, "#224444")
circ("SYNC", 24.5, 112, D_JACK, "#224444")
circ("OUT", 36.5, 112, D_JACK, "#442222")
circ("SPK", 45.5, 112, D_JACK, "#442222")
circ("W", 54.5, 112, D_JACK, "#442222")

ax.set_xlim(-2, W + 2); ax.set_ylim(H + 2, -2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Axon 12HP — FHN phase-portrait display (poly; +SYNC)", fontsize=8)
out = os.path.join(os.getcwd(), "axon_panel.png")
plt.tight_layout(); plt.savefig(out, dpi=140, bbox_inches="tight")
print("wrote", out)
print("Layout check:", "  ".join(overlaps) if overlaps else "no overlaps / out-of-bounds")
