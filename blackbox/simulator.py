"""A re-implementation of the Black Box built ONLY from the rules in REPORT.md.

Exact parts:  Area 3 (redraw from p), Area 4 (swap only if it increases matching 8-neighbours),
              one-cell / one-pair picking per step, the allowed-move sets of Areas 1 and 2.
Fitted parts: Area 1 jump probability (by branch and the cell's history, C18-C19), 4/5 exit rates (C17), Area 2 change and target
              probabilities (by old value only, C24), Area 4 acceptance by gain (C14).
All fitted numbers come from the DISCOVERY runs; the simulator is then compared with HELD-OUT runs.
"""
import numpy as np
from .regions import AREA, MASK, MOORE

CELLS = {a: np.argwhere(MASK[a]) for a in (1, 2, 3, 4)}


def fit_parameters(discovery_runs, p_area3, area4_accept, split_branches=True):
    """Estimate the Area 1 / Area 2 probabilities from discovery runs."""
    m1, m2 = MASK[1], MASK[2]
    low = lambda v: v <= 4
    # Area 1: jump probability for first move, and by time since last move (25-step bins, capped at 300)
    # branch index: 0 = low (0-3), 1 = high (6-9)
    first = np.zeros((2, 2)); byage = np.zeros((2, 14, 2))
    for a in discovery_runs:
        last = {}
        for t in range(len(a) - 1):
            for r, c in np.argwhere((a[t + 1] != a[t]) & m1):
                o, n = int(a[t, r, c]), int(a[t + 1, r, c])
                if o not in (4, 5):
                    j = low(o) != low(n); br = 0 if low(o) else 1
                    if (r, c) not in last:
                        first[br, 0] += 1; first[br, 1] += j
                    else:
                        b = min((t - last[(r, c)]) // 25, 13); byage[br, b, 0] += 1; byage[br, b, 1] += j
                last[(r, c)] = t
    if not split_branches:          # ablation: same jump rates on both branches
        first = np.repeat(first.sum(0, keepdims=True), 2, 0); byage = np.repeat(byage.sum(0, keepdims=True), 2, 0)
    jump_first = first[:, 1] / first[:, 0]
    jump_age = byage[:, :, 1] / np.maximum(byage[:, :, 0], 1)
    # Area 1: exit rate of 4 and of 5 when picked
    picked = {4: 0.0, 5: 0.0}; moved = {4: 0, 5: 0}
    for a in discovery_runs:
        v = a[:-1][:, m1]; ch = (a[1:] != a[:-1])[:, m1]
        for k in (4, 5):
            picked[k] += (v == k).sum() / 100; moved[k] += (ch & (v == k)).sum()
    exit45 = {k: moved[k] / picked[k] for k in (4, 5)}
    if not split_branches:
        exit45 = {k: (moved[4] + moved[5]) / (picked[4] + picked[5]) for k in (4, 5)}
    # Area 2: P(change | picked, old) and P(new | old, changed)
    occ = np.zeros(10); T = np.zeros((10, 10))
    for a in discovery_runs:
        v, n = a[:-1][:, m2], a[1:][:, m2]
        occ += np.bincount(v.ravel(), minlength=10)[:10]
        c = v != n
        np.add.at(T, (v[c], n[c]), 1)
    p_change2 = np.clip(100 * T.sum(1) / np.maximum(occ, 1), 0, 1)
    target2 = T / np.maximum(T.sum(1, keepdims=True), 1)
    return dict(jump_first=jump_first, jump_age=jump_age, exit45=exit45,
                p_change2=p_change2, target2=target2, p3=np.asarray(p_area3), a4=np.asarray(area4_accept))


def _gain(g, p1, p2):
    a, b = g[p1], g[p2]
    if a == b:
        return 0
    def cnt(r, c, v, skip):
        s = 0
        for dr, dc in MOORE:
            rr, cc = r + dr, c + dc
            if 0 <= rr < 20 and 0 <= cc < 20:
                x = g[rr, cc] if (rr, cc) != skip else None
                s += (x == v)
        return s
    # after the swap, the other swapped cell holds the opposite value: handle adjacency by skipping it and adding back
    adj = max(abs(p1[0] - p2[0]), abs(p1[1] - p2[1])) == 1
    before = cnt(*p1, a, None) + cnt(*p2, b, None)
    after = cnt(*p1, b, p2) + cnt(*p2, a, p1) + (2 * (a == b) if adj else 0)
    return after - before


def simulate(start, steps, P, rng):
    g = start.copy()
    frames = [g.copy()]
    last = {}             # Area 1: time of last move per cell
    for t in range(steps):
        new = g.copy()
        # Area 1
        r, c = CELLS[1][rng.integers(100)]
        v = g[r, c]
        if v in (4, 5):
            if rng.random() < P["exit45"][int(v)]:
                new[r, c] = 9 if v == 4 else 0; last[(r, c)] = t
        else:
            br = 0 if v <= 4 else 1
            pj = P["jump_first"][br] if (r, c) not in last else P["jump_age"][br, min((t - last[(r, c)]) // 25, 13)]
            if rng.random() < pj:
                new[r, c] = 9 if v <= 4 else 0
            else:
                new[r, c] = v + 1 if v <= 4 else v - 1
            last[(r, c)] = t
        # Area 2
        r, c = CELLS[2][rng.integers(100)]
        v = g[r, c]
        if v != 0 and rng.random() < P["p_change2"][v]:
            new[r, c] = rng.choice(10, p=P["target2"][v])
        # Area 3
        r, c = CELLS[3][rng.integers(100)]
        new[r, c] = rng.choice(len(P["p3"]), p=P["p3"])
        # Area 4
        i, j = rng.choice(100, 2, replace=False)
        p1, p2 = tuple(CELLS[4][i]), tuple(CELLS[4][j])
        k = _gain(g, p1, p2)
        if k > 0 and rng.random() < P["a4"][min(k, len(P["a4"])) - 1]:
            new[p1], new[p2] = g[p2], g[p1]
        g = new
        frames.append(g.copy())
    return np.array(frames)
