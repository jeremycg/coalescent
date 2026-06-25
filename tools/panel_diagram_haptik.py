#!/usr/bin/env python3
"""Haptik panel layout preview + footprint check (18 HP, top display).

Mirrors the control coordinates in src/Haptik.cpp HaptikWidget and the ring
display. Draws them at realistic sizes, flags overlaps / out-of-bounds, and
writes a preview PNG:

    python3 tools/panel_diagram.py     # writes haptik_panel.png, prints check
"""
import matplotlib, math, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

W, H = 91.44, 128.5
D_KNOB, D_TRIM, D_JACK, SCREW = 10.0, 7.0, 8.4, 5.0
fig, ax = plt.subplots(figsize=(W/10, H/10), dpi=130)
ax.add_patch(Rectangle((0,0), W, H, facecolor="#16162a", edgecolor="black"))

overlaps, placed = [], []
def circ(name, cx, cy, d, color):
    r=d/2
    if cx-r<0 or cx+r>W or cy-r<0 or cy+r>H: overlaps.append(f"OOB {name}")
    for n2,x2,y2,r2 in placed:
        if ((cx-x2)**2+(cy-y2)**2)**0.5 < (r+r2)-0.3: overlaps.append(f"OVERLAP {name}<->{n2}")
    placed.append((name,cx,cy,r))
    ax.add_patch(Circle((cx,cy), r, facecolor=color, edgecolor="white", lw=0.5, alpha=0.9))
    ax.text(cx,cy,name,ha="center",va="center",fontsize=3.6,color="white")
def sw(name,cx,cy):
    placed.append((name,cx,cy,4.5))
    ax.add_patch(FancyBboxPatch((cx-2.5,cy-4.5),5,9,boxstyle="round,pad=0.2",facecolor="#444455",edgecolor="white",lw=0.5))
    ax.text(cx,cy,name,ha="center",va="center",fontsize=3.2,color="white",rotation=90)

for x,y in [(1,1),(85.42,1),(1,122),(85.42,122)]:
    ax.add_patch(Circle((x+SCREW/2,y+SCREW/2),SCREW/2,facecolor="#888888",edgecolor="black",lw=0.4))
    placed.append(("screw",x+SCREW/2,y+SCREW/2,SCREW/2))

# ── big ring display across the top ──
dx,dy,dw,dh = 6.5, 8, 78, 46
ax.add_patch(FancyBboxPatch((dx,dy),dw,dh,boxstyle="round,pad=0.4",facecolor="#070712",edgecolor="#2b2b4d",lw=1.2))
cx,cy = dx+dw/2, dy+dh/2
halfmin=min(dw/2,dh/2); R=halfmin*0.62; amp=halfmin*0.33
# concentric guides + spokes
for f in (0.4,0.7,1.0):
    ax.add_patch(Circle((cx,cy),R*f,facecolor="none",edgecolor="#40708f",lw=0.5,alpha=0.25))
for k in range(12):
    th=2*math.pi*k/12
    ax.plot([cx,cx+R*math.cos(th)],[cy,cy+R*math.sin(th)],color="#305070",lw=0.4,alpha=0.2)
# ring shape
N=28
pts=[]
for i in range(N):
    th=2*math.pi*i/N - math.pi/2
    r=R+0.6*math.sin(2*math.pi*3*i/N)*amp
    pts.append((cx+r*math.cos(th), cy+r*math.sin(th)))
poly=plt.Polygon(pts, closed=True, facecolor="#40c0ff", edgecolor="#8ae0ff", lw=1.3, alpha=0.18)
ax.add_patch(poly)
ax.plot([p[0] for p in pts]+[pts[0][0]],[p[1] for p in pts]+[pts[0][1]],color="#8ae0ff",lw=1.2)
for p in pts: ax.add_patch(Circle(p,0.45,facecolor="#bfeeff",edgecolor="none"))
# driver + scan comet
dth=2*math.pi*7/N-math.pi/2; dr=R+0.6*math.sin(2*math.pi*3*7/N)*amp
ax.add_patch(Circle((cx+dr*math.cos(dth),cy+dr*math.sin(dth)),1.3,facecolor="#ff9b3a",edgecolor="none"))
# single scan dot (no trail — true position is audio-rate / unviewable)
sth=2*math.pi*0.32-math.pi/2
ax.plot([cx,cx+R*math.cos(sth)],[cy,cy+R*math.sin(sth)],color="#fff09a",lw=0.6,alpha=0.4)
ax.add_patch(Circle((cx+R*math.cos(sth),cy+R*math.sin(sth)),2.4,facecolor="#fff09a",edgecolor="none",alpha=0.45))
ax.add_patch(Circle((cx+R*math.cos(sth),cy+R*math.sin(sth)),1.3,facecolor="#fff2a0",edgecolor="none"))
ax.text(dx+4,dy+5,"HAPTIK",ha="left",va="center",fontsize=6,color="#9ab0ff",weight="bold")
ax.text(dx+dw-4,dy+5,"FREEZE",ha="right",va="center",fontsize=4.5,color="#6ab0ff")

# ── controls below ──
# Zone A: voice knobs
circ("N",11,64,D_KNOB,"#333344"); circ("PITCH",26,64,D_KNOB,"#333344")
circ("EXC",41,64,D_KNOB,"#333344"); circ("DRV",56,64,D_KNOB,"#333344")
sw("FRZ",70,64); sw("MODE",81,64)
# Zone B: CV channel strips (knob / att / jack)
for x,nm in [(13,"RATE"),(32,"COUP"),(51,"DAMP"),(70,"INJ")]:
    circ(nm,x,82,D_KNOB,"#333344"); circ(nm+".a",x,93,D_TRIM,"#555533"); circ(nm+".cv",x,102,D_JACK,"#224444")
# Zone C: I/O
circ("VOCT",13,113,D_JACK,"#224444"); circ("TRIG",27,113,D_JACK,"#224444"); circ("EXT",41,113,D_JACK,"#224444")
circ("OUT",64,113,D_JACK,"#442222"); circ("MOT",78,113,D_JACK,"#442222")

ax.set_xlim(-2,W+2); ax.set_ylim(H+2,-2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Haptik 18HP — top display redesign",fontsize=8)
out=os.path.join(os.getcwd(),"haptik_panel.png")
plt.tight_layout(); plt.savefig(out,dpi=140,bbox_inches="tight")
print("wrote",out)
print("Layout check:", "  ".join(overlaps) if overlaps else "no overlaps / out-of-bounds")
