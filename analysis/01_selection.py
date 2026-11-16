"""01 - Areas, starting grids, and HOW cells are picked each step.

Claims proven here: C01-C05 (see REPORT.md).
Run:  python analysis/01_selection.py
"""
import sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from scipy.stats import chisquare
from blackbox.data import load_discovery, load_heldout, load_group2
from blackbox.regions import AREA, MASK, VN
from blackbox.events import changes_per_step, single_events
from blackbox.report import Log, FIGURES
from blackbox.plotting import plt, AREA_COLOR, INK

log = Log("01_selection")
D, H, G2 = load_discovery(), load_heldout(), load_group2()
log(f"datasets: discovery={len(D)} runs x {D[0].shape[0]} frames, heldout={len(H)} runs x {H[0].shape[0]} frames, group2={len(G2)} run x {G2[0].shape[0] if G2 else 0} frames")

# ---------------------------------------------------------------- C01 area shapes
sizes = {a: int(MASK[a].sum()) for a in (1, 2, 3, 4)}
cnt = np.zeros((20, 20))
for a in D:
    cnt += (a[1:] != a[:-1]).sum(0)
per_area = {a: float(cnt[MASK[a]].mean() / len(D)) for a in (1, 2, 3, 4)}
# how sharp is each boundary? mean change count of cells just inside vs the neighbouring area
log.claim("C01", f"Four L-shaped areas of 100 cells each {sizes}; mean changes per cell per 2,000-step run "
          f"(discovery) = " + ", ".join(f"Area {a}: {v:.1f}" for a, v in per_area.items()),
          sizes=sizes, mean_changes_per_cell=per_area)

fig, ax = plt.subplots(figsize=(5.6, 5))
im = ax.imshow(cnt / len(D), cmap="Blues")
for a in (1, 2, 3, 4):
    ys, xs = np.nonzero(MASK[a])
    ax.text(xs.mean(), ys.mean(), f"Area {a}", ha="center", va="center", color=INK, fontsize=11, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8))
ax.contour(AREA, levels=[1.5, 2.5, 3.5], colors="white", linewidths=2)
ax.set_title("How often each cell changed (per run, 2,000 steps)")
ax.set_xticks(range(0, 20, 5)); ax.set_yticks(range(0, 20, 5)); ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.8, label="changes per run")
fig.tight_layout(); fig.savefig(FIGURES / "01_change_map.png"); plt.close(fig)

# ---------------------------------------------------------------- C02 cells changed per step
res = {}
for name, runs in [("discovery", D), ("heldout", H), ("group2", G2)]:
    if not runs:
        continue
    res[name] = {a: dict(sorted(collections.Counter(changes_per_step(runs, a)).items())) for a in (1, 2, 3, 4)}
    for a in (1, 2, 3, 4):
        log(f"   {name:9s} Area {a}: cells changed per step -> number of steps {res[name][a]}")
max_change = {name: {a: max(v[a]) for a in v} for name, v in res.items()}
area4_values = {name: sorted(v[4]) for name, v in res.items()}
log.claim("C02", "In every dataset, Areas 1-3 never change more than 1 cell in a step, and Area 4 changes either 0 or exactly 2 cells.",
          max_cells_changed=max_change, area4_counts_seen=area4_values, table=res)

# ---------------------------------------------------------------- C03 random (not raster) order
# If cells were updated in a top-to-bottom / left-to-right sweep, consecutive changes would hit
# consecutive cells. Measure how often the next change in Area 3 is the next cell in raster order.
idx = {tuple(rc): i for i, rc in enumerate(np.argwhere(MASK[3]))}
seq_next = total = 0
counts_by_cell = np.zeros(100)
for a in D:
    prev = None
    for t in range(len(a) - 1):
        ch = np.argwhere((a[t + 1] != a[t]) & MASK[3])
        if len(ch) != 1:
            continue
        i = idx[tuple(ch[0])]
        counts_by_cell[i] += 1
        if prev is not None:
            total += 1
            seq_next += (i == (prev + 1) % 100)
        prev = i
frac_next = seq_next / total
p_uniform = chisquare(counts_by_cell).pvalue
log.claim("C03", f"Cells are picked in random order, not a sweep: in Area 3 the next changed cell is the next cell in "
          f"reading order only {frac_next:.2%} of the time (a sweep would give ~100%, random picking ~1%). "
          f"Changes are spread evenly over the 100 cells (chi-square vs uniform p = {p_uniform:.2f}).",
          frac_next_in_raster=frac_next, chi2_uniform_p=p_uniform, n_pairs=total)

# ---------------------------------------------------------------- C04 starting grids
starts = [r[0] for r in D] + [h[0] for h in H]
distinct = len({s.tobytes() for s in starts})
share = np.bincount(np.concatenate([s.ravel() for s in starts]), minlength=10)[:10]
share = share / share.sum()
log.claim("C04", f"Reset gives a new random grid every time: {distinct} distinct starting grids out of {len(starts)} resets; "
          f"digits 0-9 each fill {share.min():.1%}-{share.max():.1%} of a starting grid.",
          distinct_starts=distinct, n_starts=len(starts), start_digit_share=share)

# ---------------------------------------------------------------- C05 change does not depend on matching neighbours
# P(change | cell picked) ~ 100 x P(cell changes in a step), because one of 100 cells is picked per step (C02).
table = {}
for area in (1, 2, 3):
    m = MASK[area]
    tot = collections.Counter(); ch = collections.Counter()
    for a in D:
        for t in range(len(a) - 1):
            g = a[t]
            eq = np.zeros_like(g)
            for dr, dc in VN:
                s = np.full(g.shape, -1)
                r0, r1 = max(0, -dr), min(20, 20 - dr); c0, c1 = max(0, -dc), min(20, 20 - dc)
                s[r0:r1, c0:c1] = g[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
                eq += (s == g)
            v = g[m]; e = eq[m]; c = (a[t + 1] != g)[m]
            for vv, ee, cc in zip(v.tolist(), e.tolist(), c.tolist()):
                tot[(vv, ee)] += 1; ch[(vv, ee)] += cc
    table[area] = {f"{v}|{e}": round(100 * ch[(v, e)] / tot[(v, e)], 3)
                   for (v, e) in sorted(tot) if tot[(v, e)] >= 5000 and v <= 9}
    log(f"   Area {area}: P(change | picked) by 'value|#side-neighbours equal to it' (n>=5000): {table[area]}")
# summary: for Area 3 and Area 1 non-terminal values, spread across neighbour counts
def spread(area, values):
    out = {}
    for v in values:
        xs = [p for k, p in table[area].items() if int(k.split('|')[0]) == v]
        if len(xs) >= 2:
            out[v] = max(xs) - min(xs)
    return out
sp3 = spread(3, range(10)); sp1 = spread(1, [0, 1, 2, 3, 6, 7, 8, 9])
log.claim("C05", "Whether a picked cell changes does NOT depend on how many neighbours match it: in Area 3 and for "
          "Area 1's non-4/5 values the change probability is flat across 0-4 matching side-neighbours "
          f"(largest spread {100*max(list(sp3.values()) + list(sp1.values())):.0f} percentage points; Area 1 non-4/5 values change ~100% of the time when picked).",
          table=table, spread_area3=sp3, spread_area1=sp1,
          area1_nonterminal_min=min(p for k, p in table[1].items() if int(k.split('|')[0]) not in (4, 5)),
          area1_nonterminal_max=max(p for k, p in table[1].items() if int(k.split('|')[0]) not in (4, 5)))
log.save()
