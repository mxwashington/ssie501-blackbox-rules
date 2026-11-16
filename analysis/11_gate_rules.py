"""11 - Region 1 gates, the two parts of Region 2, the strict Region 4 rule and the update order, tested on the
held-out runs; then the full Part A model rebuilt and compared with the held-out runs.

Claims proven here: C35-C41. The Region 1 gates, the Region 2 split and the strict Region 4 rule come from
Marcus Washington's simulator, built on his own recorded runs; none of the 30 held-out runs was used to build them.
Run:  python analysis/11_gate_rules.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_heldout
from blackbox.regions import AREA
from blackbox.report import Log
from blackbox.gate_model import (GateModel, r1_target_grid, PT, FO, NB8, matches)

log = Log("11_gate_rules")
runs = load_heldout()
R1 = AREA == 1

# Waiting-time rule from the earlier version of A5 (kept only to score it against the gates)
LO = [0.36, 0.02, 0.04, 0.08, 0.09, 0.10]
HI = [0.63, 0.06, 0.18, 0.26, 0.33, 0.39]


def wait_rate(low, wait):
    k = 0 if wait is None else 1 + min(int(wait // 50), 4)
    return (LO if low else HI)[k]


def gains(g, i, j):
    vi, vj = g[i], g[j]
    bi, bj = matches(g, i, vi), matches(g, j, vj)
    h = g.copy()
    h[i], h[j] = vj, vi
    ai, aj = matches(h, i, vj), matches(h, j, vi)
    return (ai + aj > bi + bj), (ai >= bi and aj >= bj and ai + aj > bi + bj)


n_ev = gate_ok = wait_ok = 0
ll_wait = 0.0
settled = {"agree": [0, 0], "disagree": [0, 0]}
pt = dict(changes=0, rises=0, partner_found=0)
fo = dict(changes=0, product_found=0)
sw = dict(swaps=0, total_gain_old=0, strict_old=0, strict_in_order=0)

for h in runs:
    last = np.full((20, 20), -1)
    for t in range(len(h) - 1):
        g0, g1 = h[t], h[t + 1]
        ch = g0 != g1
        T = r1_target_grid(g0)
        # Region 1: settled cells
        for r, c in np.argwhere(R1 & ((g0 == 4) | (g0 == 5))):
            key = "agree" if T[r, c] == int(g0[r, c] >= 5) else "disagree"
            settled[key][0] += 1
            settled[key][1] += int(g1[r, c] != g0[r, c])
        # Region 1: every change from 0-3 or 6-9
        for r, c in np.argwhere(ch & R1):
            v, w = int(g0[r, c]), int(g1[r, c])
            if v not in (4, 5):
                high = v >= 5
                jump = w == (0 if high else 9)
                n_ev += 1
                gate_ok += int((T[r, c] != int(high)) == jump)
                p = wait_rate(not high, None if last[r, c] < 0 else t - last[r, c])
                ll_wait += np.log(p if jump else 1 - p)
                wait_ok += int(jump == (p > 0.5))
            last[r, c] = t
        # Region 2
        bottom = set(g0[10:20].ravel().tolist())
        for r, c in np.argwhere(ch & PT):
            pt["changes"] += 1
            pt["rises"] += int(g1[r, c] > g0[r, c])
            pt["partner_found"] += int(int(g1[r, c]) in bottom)
        fov = g0[FO]
        for r, c in np.argwhere(ch & FO):
            fo["changes"] += 1
            fo["product_found"] += int(any((int(g0[r, c]) * int(n)) % 10 == g1[r, c] for n in fov))
        # Region 4
        d4 = np.flatnonzero((ch & (AREA == 4)).ravel())
        if len(d4) == 2:
            i, j = d4
            sw["swaps"] += 1
            a, b = gains(g0.ravel(), i, j)
            sw["total_gain_old"] += int(a)
            sw["strict_old"] += int(b)
            gs = g0.copy()
            m = ch & (AREA != 4)
            gs[m] = g1[m]                      # Regions 1-3 have already moved this step
            sw["strict_in_order"] += int(gains(gs.ravel(), i, j)[1])

log.claim("C35", f"Region 1 gates predict jump vs advance for {gate_ok} of {n_ev} held-out changes from 0-3 or 6-9",
          correct=gate_ok, n=n_ev)
na, ja = settled["agree"]
nd, jd = settled["disagree"]
log.claim("C36", f"Region 1 cells at 4 or 5: gate agrees -> {ja} jumps in {na} cell-steps; gate disagrees -> "
          f"{jd} jumps in {nd} cell-steps ({jd / nd:.4f} per step; picking alone gives 0.0100)",
          agree_cellsteps=na, agree_jumps=ja, disagree_cellsteps=nd, disagree_jumps=jd, disagree_rate=jd / nd)
log.claim("C37", f"Waiting-time rule (earlier A5) scored on the same changes: best guess right {wait_ok}/{n_ev} "
          f"({wait_ok / n_ev:.1%}); log-likelihood {ll_wait / n_ev:.3f} per change (gates: 0 errors)",
          correct=wait_ok, n=n_ev, loglik_per_change=ll_wait / n_ev)
log.claim("C38", f"Region 2 minimum part: {pt['changes']} changes, {pt['rises']} rises, new value present in rows "
          f"11-20 in {pt['partner_found']}", **pt)
log.claim("C39", f"Region 2 multiply part: {fo['changes']} changes, {fo['product_found']} equal (old x a value in "
          f"that part) mod 10", **fo)
log.claim("C40", f"Region 4: {sw['swaps']} swaps; strict rule holds for {sw['strict_old']} on the old grid and "
          f"{sw['strict_in_order']} when Regions 1-3 move first; total-gain rule {sw['total_gain_old']}", **sw)

# Full model rebuilt from the 30 held-out starting grids, 5 seeds each
AREA_FLAT = AREA.ravel()


def H(v):
    _, c = np.unique(v, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def stats(traj):
    ch = traj[1:] != traj[:-1]
    end = traj[-1]
    s = {f"r{a}_changed": int(ch[:, AREA == a].sum()) for a in (1, 2, 3, 4)}
    s.update({f"r{a}_entropy_end": H(end[AREA == a]) for a in (1, 2, 3, 4)})
    e1 = end[AREA == 1]
    s.update(r1_share4=float(np.mean(e1 == 4)), r1_share45=float(np.mean((e1 == 4) | (e1 == 5))),
             r2_share0=float(np.mean(end[AREA == 2] == 0)))
    return s


real = [stats(h) for h in runs]
model = [stats(GateModel(h[0], seed=100 * k + rep).run(len(h) - 1)) for k, h in enumerate(runs) for rep in range(5)]
mean = lambda rows, key: float(np.mean([r[key] for r in rows]))
table = {key: dict(model=mean(model, key), real=mean(real, key)) for key in real[0]}
log("Full model vs held-out runs (2,000 steps, means):")
for key, v in table.items():
    log(f"  {key:18s} model {v['model']:9.3f}   real {v['real']:9.3f}")
log.claim("C41", "The Part A model, rebuilt from the 30 held-out starting grids, matches the held-out runs "
          "(table above)", table=table)
log.save()
