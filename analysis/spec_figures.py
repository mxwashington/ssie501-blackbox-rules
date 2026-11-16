"""Draw the two figures of MODEL_SPECIFICATION.md (not a claim script; run_all.py does not run it).

    python analysis/spec_figures.py   ->  figures/11_region_map.png, figures/11_region1_moves.png
"""
from pathlib import Path
FIG = Path(__file__).resolve().parent.parent / "figures"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
plt.rcParams["font.family"]="DejaVu Sans"
# ---------- Figure 1: region map ----------
COL={"r1":"#b9cdf2","r2min":"#f4b9a3","r2mul":"#f9dccf","r3":"#b5e3cf","r4":"#a9afcc"}
EDGE={"r1":"#6f93d8","r2min":"#e08a6b","r2mul":"#e9ae95","r3":"#62c09a","r4":"#6a7299"}
def part(r,c):
    if r<5 and c<15 or 5<=r<10 and c<5: return "r1"
    if r<10 and c>=15: return "r2min"
    if 5<=r<10: return "r2mul"
    if 10<=r<15 and c<15 or r>=15 and c<5: return "r3"
    return "r4"
inputs={(8,8),(8,13),(8,18),(13,3),(13,8),(13,13),(18,8),(18,13),(18,18),(3,18),(18,3),(13,18)}
fig=plt.figure(figsize=(8.6,6.2),dpi=200)
ax=fig.add_axes([0.04,0.08,0.58,0.78])
for r in range(20):
    for c in range(20):
        p=part(r,c)
        ax.add_patch(FancyBboxPatch((c+0.08,r+0.08),0.84,0.84,boxstyle="round,pad=0,rounding_size=0.18",fc=COL[p],ec=EDGE[p],lw=0.6))
        if (r+1,c+1) in inputs:
            ax.add_patch(FancyBboxPatch((c+0.08,r+0.08),0.84,0.84,boxstyle="round,pad=0,rounding_size=0.18",fc="none",ec="#1d1d1b",lw=1.8))
for x in (5,10,15):
    ax.plot([x,x],[0,10] if x!=5 else [0,10],color="white",lw=0)  # spacer
# block letters for Region 1
for lab,(r,c) in {"A":(2.5,2.5),"B":(2.5,7.5),"C":(2.5,12.5),"D":(7.5,2.5)}.items():
    ax.text(c,r,lab,ha="center",va="center",fontsize=13,fontweight="bold",color="#1f3f7a")
for lab,(r,c) in {"2":(7.5,10.5),"3":(12.5,10.5),"4":(17.5,10.5)}.items():
    ax.text(c,r,lab,ha="center",va="center",fontsize=12,fontweight="bold",color="#1d1d1b")
ax.text(17.5,5,"2",ha="center",va="center",fontsize=12,fontweight="bold",color="#1d1d1b")
ax.set_xlim(0,20); ax.set_ylim(20,0); ax.set_aspect("equal")
ax.set_xticks([0.5,4.5,9.5,14.5,19.5]); ax.set_xticklabels([1,5,10,15,20],fontsize=8)
ax.set_yticks([0.5,4.5,9.5,14.5,19.5]); ax.set_yticklabels([1,5,10,15,20],fontsize=8)
ax.xaxis.tick_top(); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel("column",fontsize=9); ax.set_ylabel("row",fontsize=9)
fig.text(0.04,0.955,"Four L-shaped regions of 100 cells each",fontsize=13,fontweight="bold")
fig.text(0.04,0.915,"Region 1 is drawn as its four 5 × 5 blocks A-D; Region 2 is split into its two parts",fontsize=8.5,color="#555")
leg=[("r1","Region 1 (blocks A-D)","counters steered by logic gates;\nfreezes as 4s and 5s"),
     ("r2min","Region 2, minimum part","keeps the smaller of itself and a\nrandom cell of rows 11-20; ends as 0s"),
     ("r2mul","Region 2, multiply part","multiplies by a random cell of its\npart, keeps the last digit; ends as 0s"),
     ("r3","Region 3","redrawn from a fixed table;\nnever freezes"),
     ("r4","Region 4","swaps pairs only if neither cell\nloses matching neighbours")]
y=0.84
for k,t,d in leg:
    fig.patches.append(FancyBboxPatch((0.655,y-0.012),0.025,0.033,boxstyle="round,pad=0,rounding_size=0.004",fc=COL[k],ec=EDGE[k],lw=0.8,transform=fig.transFigure))
    fig.text(0.69,y+0.004,t,fontsize=9,fontweight="bold",va="center")
    fig.text(0.69,y-0.045,d,fontsize=7.8,color="#444",va="center")
    y-=0.14
fig.patches.append(FancyBboxPatch((0.655,y-0.012),0.025,0.033,boxstyle="round,pad=0,rounding_size=0.004",fc="white",ec="#1d1d1b",lw=1.8,transform=fig.transFigure))
fig.text(0.69,y+0.004,"Outlined cells",fontsize=9,fontweight="bold",va="center")
fig.text(0.69,y-0.045,"the 12 cells outside Region 1\nthat its gates read (p, q, c)",fontsize=7.8,color="#444",va="center")
fig.savefig(FIG / "11_region_map.png",facecolor="white")
# ---------- Figure 2: Region 1 moves ----------
fig,ax=plt.subplots(figsize=(8.6,3.6),dpi=200); ax.set_xlim(0,10); ax.set_ylim(-0.2,4.4); ax.axis("off")
def lane(y,label):
    ax.add_patch(FancyBboxPatch((1.0,y-0.45),5.9,0.9,boxstyle="round,pad=0,rounding_size=0.15",fc="#f2f2f0",ec="#bbb",lw=0.8))
    ax.text(1.05,(y+0.6) if y>2 else (y-0.75),label,fontsize=8.5,color="#555")
def node(x,y,t,end=False):
    ax.add_patch(plt.Circle((x,y),0.3,fc="#dfe8f8" if end else "white",ec="#4b74c9" if end else "#888",lw=1.6 if end else 1))
    ax.text(x,y,t,ha="center",va="center",fontsize=11,fontweight="bold" if end else None)
def arr(x1,y1,x2,y2,dash=False,rad=0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=10,lw=1,color="#666",linestyle=(0,(3,3)) if dash else "-",connectionstyle=f"arc3,rad={rad}"))
yl,yh=3.2,1.1
lane(yl,"low branch (0-4)"); lane(yh,"high branch (5-9)")
xs=[1.6,2.95,4.3,5.65]
for x,t in zip(xs,"0123"): node(x,yl,t)
for x,t in zip(xs,"9876"): node(x,yh,t)
node(7.6,yl,"4",True); node(7.6,yh,"5",True)
for a,b in zip(xs,xs[1:]+[7.6]):
    arr(a+0.3,yl,b-0.3,yl); arr(a+0.3,yh,b-0.3,yh)
arr(2.95,yl-0.32,1.75,yh+0.32,dash=True)
ax.text(2.55,2.15,"jump to 9 when\nthe gate says HIGH",fontsize=8,ha="left",va="center",color="#333")
arr(1.45,yh+0.32,1.45,yl-0.32,dash=True)
ax.text(0.05,2.15,"jump to 0\nwhen the gate\nsays LOW",fontsize=8,ha="left",va="center",color="#333")
ax.text(8.05,yl,"stays while its gate says LOW;\njumps to 9 when it says HIGH",fontsize=8,va="center")
ax.text(8.05,yh,"stays while its gate says HIGH;\njumps to 0 when it says LOW",fontsize=8,va="center")
ax.text(0.05,4.25,"Region 1: advance along the branch, or jump when the gate points to the other half",fontsize=11,fontweight="bold")
ax.text(0.05,-0.05,"Solid arrow: advance one step (gate agrees with the current half). Dashed arrow: jump (gate points to the other half). Which cell is picked is random; what it does is not.",fontsize=7.6,color="#555")
fig.savefig(FIG / "11_region1_moves.png",facecolor="white",bbox_inches="tight",pad_inches=0.1)
