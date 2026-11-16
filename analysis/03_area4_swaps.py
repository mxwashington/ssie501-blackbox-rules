"""03 - Area 4: two random cells swap, and only when the swap increases matching neighbours.

Claims proven here: C11-C15 (see REPORT.md).
Run:  python analysis/03_area4_swaps.py      (about 5-8 minutes)
"""
import sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_discovery, load_heldout, load_group2
from blackbox.regions import AREA, MASK, MOORE, VN
from blackbox.events import swap_events
from blackbox.report import Log, FIGURES
from blackbox.plotting import plt, AREA_COLOR, INK, MUTED

log = Log("03_area4_swaps")
D, H, G2 = load_discovery(), load_heldout(), load_group2()
CELLS = np.argwhere(MASK[4]); K = len(CELLS)
IDX = {tuple(x): i for i, x in enumerate(CELLS)}

# ---------------------------------------------------------------- C11 swaps conserve counts
exact = total = 0; counts_change = 0
for runs in (D, H, G2):
    for ri, t, p1, p2 in swap_events(runs):
        a = runs[ri]; total += 1
        exact += (a[t][p1] == a[t + 1][p2]) and (a[t][p2] == a[t + 1][p1])
    for a in runs:
        v = a[:, MASK[4]]
        hist = np.stack([np.bincount(row, minlength=14) for row in v])
        counts_change += int((hist != hist[0]).any(1).sum())
log.claim("C11", f"Every Area 4 change is two cells trading values ({exact} of {total} two-cell changes are exact swaps, all datasets), "
          f"so Area 4's count of each value never changes (frames whose value counts differ from the run's first frame: {counts_change}).",
          exact=exact, total=total, frames_with_changed_counts=counts_change)


# ---------------------------------------------------------------- energy definitions
def like_count(g, r, c, v, nb, wrap, kind):
    s = 0
    for dr, dc in nb:
        rr, cc = r + dr, c + dc
        if wrap:
            rr %= 20; cc %= 20
        elif not (0 <= rr < 20 and 0 <= cc < 20):
            continue
        x = int(g[rr, cc])
        s += (x == v) if kind == "equal" else (-abs(x - v) if kind == "absdiff" else -(x - v) ** 2)
    return s


def gain(g, p1, p2, nb=MOORE, wrap=False, kind="equal"):
    """Change in (score of cell 1 + score of cell 2) if the two values were swapped."""
    a, b = int(g[p1]), int(g[p2])
    before = like_count(g, *p1, a, nb, wrap, kind) + like_count(g, *p2, b, nb, wrap, kind)
    g2 = g.copy(); g2[p1], g2[p2] = b, a
    after = like_count(g2, *p1, b, nb, wrap, kind) + like_count(g2, *p2, a, nb, wrap, kind)
    return after - before


# ---------------------------------------------------------------- C12 every swap strictly improves (8 neighbours, no wrap)
res = {}
for name, runs in [("discovery", D), ("heldout", H), ("group2", G2)]:
    n = strict = fixed = viol = 0; boundary = 0
    for ri, t, p1, p2 in swap_events(runs):
        a = runs[ri]; n += 1
        if gain(a[t], p1, p2) > 0:
            strict += 1
        else:
            # Other areas also update in the same step; if the box applies them first, use their new values.
            g2 = a[t].copy(); out = AREA != 4; g2[out] = a[t + 1][out]
            if gain(g2, p1, p2) > 0:
                fixed += 1
                # does either swapped cell touch a cell of another area?
                boundary += any(0 <= q[0] + dr < 20 and 0 <= q[1] + dc < 20 and AREA[q[0] + dr, q[1] + dc] != 4
                                for q in (p1, p2) for dr, dc in MOORE)
            else:
                viol += 1
    res[name] = dict(swaps=n, improve_at_t=strict, improve_after_other_areas=fixed, of_which_touch_boundary=boundary, violations=viol)
    log(f"   {name}: {res[name]}")
log.claim("C12", "Every Area 4 swap strictly increases the number of same-value neighbours (8 surrounding cells, grid edges not wrapped, "
          "cells of other areas included): "
          + "; ".join(f"{k} {v['improve_at_t'] + v['improve_after_other_areas']}/{v['swaps']}" for k, v in res.items())
          + f". {res['discovery']['improve_after_other_areas']} discovery swaps only count as improving once the same-step update of a "
          f"neighbouring area is applied; {res['discovery']['of_which_touch_boundary']} of them touch an area boundary.", **res)

# ---------------------------------------------------------------- C13 alternative definitions fail
alts = {}
for nbn, nb in [("8-neighbour", MOORE), ("4-neighbour", VN)]:
    for wrap in (False, True):
        for kind in ("equal", "absdiff", "squared"):
            key = f"{nbn}, {'wrap' if wrap else 'no wrap'}, {kind}"
            bad = sum(gain(H[ri][t], p1, p2, nb, wrap, kind) <= 0 for ri, t, p1, p2 in swap_events(H))
            alts[key] = bad
for k, v in alts.items():
    log(f"   held-out swaps NOT improving under '{k}': {v}")
ok = [k for k, v in alts.items() if v == 0]
log.claim("C13", f"Of 12 candidate definitions of 'better' (4 or 8 neighbours, wrapped or not, matching vs value-difference), only "
          f"{ok} has zero non-improving swaps in the held-out runs; the others have {min(v for v in alts.values() if v > 0)}-{max(alts.values())}.",
          alternatives=alts, zero_exception_definitions=ok)

# ---------------------------------------------------------------- vectorised gain for ALL pairs (8-neighbour, no wrap, equal)
NB = [[(r + dr, c + dc) for dr, dc in MOORE if 0 <= r + dr < 20 and 0 <= c + dc < 20] for r, c in CELLS]
ADJ = np.zeros((K, K), int)
for i, (r, c) in enumerate(CELLS):
    for j, (r2, c2) in enumerate(CELLS):
        ADJ[i, j] = int(i != j and max(abs(r - r2), abs(c - c2)) == 1)
IU = np.triu_indices(K, 1)


def all_pair_gains(g):
    vals = g[CELLS[:, 0], CELLS[:, 1]]
    N = np.zeros((K, 16), int)
    for i, nb in enumerate(NB):
        for rr, cc in nb:
            N[i, g[rr, cc]] += 1
    ar = np.arange(K)
    Nib = N[ar[:, None], vals[None, :]]          # cell i's neighbours equal to j's value
    Nja = N[ar[None, :], vals[:, None]]          # cell j's neighbours equal to i's value
    own = N[ar, vals]
    G = Nib + Nja - own[:, None] - own[None, :] - 2 * ADJ
    G[vals[:, None] == vals[None, :]] = 0
    return G


def exposure_and_swaps(runs, every=1):
    """Per step: share of all pairs with gain k (k=1..6+), and the gain of the swap that happened (if any)."""
    rows = []
    for ri, a in enumerate(runs):
        for t in range(0, len(a) - 1, every):
            G = all_pair_gains(a[t]); g = np.clip(G[IU], -1, 6)
            share = np.array([(g == k).mean() for k in range(1, 7)])
            ch = np.argwhere((a[t + 1] != a[t]) & MASK[4])
            k = 0
            if len(ch) == 2:
                k = int(np.clip(G[IDX[tuple(ch[0])], IDX[tuple(ch[1])]], -1, 6))
            rows.append((ri, t, share, k))
    return rows


log("   computing all-pair gains (every step of every run) ...")
RD = exposure_and_swaps(D, every=1)
RH = exposure_and_swaps(H, every=1)

# ---------------------------------------------------------------- C14 swap probability rises with the gain
def acceptance(rows):
    E = np.zeros(6); O = np.zeros(6)
    for _, _, share, k in rows:
        E += share
        if k >= 1:
            O[k - 1] += 1
    return O / E, O, E
accD, OD, ED = acceptance(RD)
accH, OH, EH = acceptance(RH)
log.claim("C14", "Not every improving pair is swapped: the chance that an improving pair is swapped (if one random pair is tried per step) "
          "rises with the size of the improvement. Held-out: "
          + ", ".join(f"+{k + 1}: {accH[k]:.2f}" for k in range(5))
          + "; discovery: " + ", ".join(f"+{k + 1}: {accD[k]:.2f}" for k in range(5)) + ".",
          heldout=accH[:5], discovery=accD[:5], heldout_swaps_by_gain=OH, heldout_exposure=EH)

# ---------------------------------------------------------------- C15 the slow-down comes from Area 4's own state
a_fit = accD.copy(); a_fit[5] = min(1.0, a_fit[5])
windows = [(0, 200), (200, 600), (600, 1200), (1200, 2000)]
cmp = {}
t_arr = np.array([t for _, t, _, _ in RH]); sh = np.stack([s for _, _, s, _ in RH]); sw = np.array([k >= 1 for *_, k in RH])
pred_step = sh @ a_fit
for lo, hi in windows:
    sel = (t_arr >= lo) & (t_arr < hi)
    cmp[f"{lo}-{hi}"] = (float(pred_step[sel].mean()), float(sw[sel].mean()))
log.claim("C15", "Area 4's slow-down can be predicted from Area 4's own grid alone (improving swaps run out as clusters form): using the "
          "discovery acceptance rates and nothing from Areas 1-2, the predicted held-out swap rate per step is "
          + "; ".join(f"steps {k}: predicted {p:.3f} vs observed {o:.3f}" for k, (p, o) in cmp.items()) + ".",
          windows=cmp)

# ---------------------------------------------------------------- figures
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
ax = axs[0]
ks = np.arange(1, 6)
ax.bar(ks - 0.2, accD[:5], width=0.38, color=AREA_COLOR[4], label="discovery runs")
ax.bar(ks + 0.2, accH[:5], width=0.38, color="#f6cf7a", label="held-out runs")
ax.set_xlabel("improvement in matching neighbours if swapped"); ax.set_ylabel("share of such pairs that swap")
ax.set_title("Bigger improvements are swapped more often"); ax.set_ylim(0, 1.05); ax.legend(loc="upper left", fontsize=8)
ax = axs[1]
bins = np.arange(0, 2000, 50)
centers = bins + 25
pb = [pred_step[(t_arr >= b) & (t_arr < b + 50)].mean() for b in bins]
ob = [sw[(t_arr >= b) & (t_arr < b + 50)].mean() for b in bins]
ax.plot(centers, ob, color=AREA_COLOR[4], label="observed (held-out)")
ax.plot(centers, pb, color=INK, linestyle="--", linewidth=1.5, label="predicted from Area 4's own grid")
ax.set_xlabel("step"); ax.set_ylabel("swaps per step"); ax.set_title("Area 4's slow-down follows from its own grid")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIGURES / "03_area4_swaps.png"); plt.close(fig)
log.save()
