"""Score the Part A rules on every recorded change of a list of runs (used by analysis/11 and 12)."""
import numpy as np

from blackbox.regions import AREA
from blackbox.gate_model import r1_target_grid, PT, FO, matches

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


def evaluate(runs):
    """Score every rule on every recorded change of a list of runs."""
    n_ev = gate_ok = wait_ok = 0
    ll_wait = 0.0
    settled = {"agree": [0, 0], "disagree": [0, 0]}
    pt = dict(changes=0, rises=0, partner_found=0)
    fo = dict(changes=0, product_found=0)
    sw = dict(swaps=0, total_gain_old=0, strict_old=0, strict_in_order=0)
    first = dict(transitions=0, r1_changes=0, r1_gate_ok=0)
    for h in runs:
        last = np.full((20, 20), -1)
        for t in range(len(h) - 1):
            g0, g1 = h[t], h[t + 1]
            ch = g0 != g1
            T = r1_target_grid(g0)
            for r, c in np.argwhere(R1 & ((g0 == 4) | (g0 == 5))):
                key = "agree" if T[r, c] == int(g0[r, c] >= 5) else "disagree"
                settled[key][0] += 1
                settled[key][1] += int(g1[r, c] != g0[r, c])
            if t == 0:
                first["transitions"] += 1
            for r, c in np.argwhere(ch & R1):
                v, w = int(g0[r, c]), int(g1[r, c])
                if v not in (4, 5):
                    high = v >= 5
                    jump = w == (0 if high else 9)
                    ok = (T[r, c] != int(high)) == jump
                    n_ev += 1
                    gate_ok += int(ok)
                    if t == 0:
                        first["r1_changes"] += 1
                        first["r1_gate_ok"] += int(ok)
                    p = wait_rate(not high, None if last[r, c] < 0 else t - last[r, c])
                    ll_wait += np.log(p if jump else 1 - p)
                    wait_ok += int(jump == (p > 0.5))
                last[r, c] = t
            bottom = set(g0[10:20].ravel().tolist())
            for r, c in np.argwhere(ch & PT):
                pt["changes"] += 1
                pt["rises"] += int(g1[r, c] > g0[r, c])
                pt["partner_found"] += int(int(g1[r, c]) in bottom)
            fov = g0[FO]
            for r, c in np.argwhere(ch & FO):
                fo["changes"] += 1
                fo["product_found"] += int(any((int(g0[r, c]) * int(n)) % 10 == g1[r, c] for n in fov))
            d4 = np.flatnonzero((ch & (AREA == 4)).ravel())
            if len(d4) == 2:
                i, j = d4
                sw["swaps"] += 1
                a, b = gains(g0.ravel(), i, j)
                sw["total_gain_old"] += int(a)
                sw["strict_old"] += int(b)
                gs = g0.copy()
                m = ch & (AREA != 4)
                gs[m] = g1[m]
                sw["strict_in_order"] += int(gains(gs.ravel(), i, j)[1])
    return dict(n_ev=n_ev, gate_ok=gate_ok, wait_ok=wait_ok, ll_wait=ll_wait, settled=settled, pt=pt, fo=fo, sw=sw, first=first)


