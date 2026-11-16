"""10 - Real worked examples, one per region, taken from held-out run 1 (shared dataset trajectory_001).

Claim proven here: C34. Prints, for one real change in each region, the before/after values and how the
written rule accounts for it (for Region 4: the two 3x3 neighbourhoods and the neighbour counts).
Run:  python analysis/10_worked_example.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_heldout
from blackbox.regions import MASK, MOORE
from blackbox.report import Log

log = Log("10_worked_example")
h = load_heldout()[0]          # trajectory_001; frame f is the grid after step f+1 (recorded_steps start at 1)
START = 200


def first_change(area, cond=lambda o, n: True, k=1, start=START):
    for t in range(start, len(h) - 1):
        ch = np.argwhere((h[t + 1] != h[t]) & MASK[area])
        if len(ch) == k:
            vals = [(tuple(map(int, p)), int(h[t][tuple(p)]), int(h[t + 1][tuple(p)])) for p in ch]
            if k == 2 or cond(vals[0][1], vals[0][2]):
                return t, vals
    raise RuntimeError


def nbhd(g, r, c):
    rows = []
    for dr in (-1, 0, 1):
        row = []
        for dc in (-1, 0, 1):
            rr, cc = r + dr, c + dc
            row.append("." if not (0 <= rr < 20 and 0 <= cc < 20) else ("[" + str(g[rr, cc]) + "]" if (dr, dc) == (0, 0) else str(g[rr, cc])))
        rows.append(" ".join(f"{x:>3}" for x in row))
    return rows


def count(g, r, c, v):
    return sum(1 for dr, dc in MOORE if 0 <= r + dr < 20 and 0 <= c + dc < 20 and g[r + dr, c + dc] == v)


ex = {}
t, [(p, o, n)] = first_change(1, lambda o, n: o not in (4, 5) and (o <= 4) == (n <= 4))
ex["region1_advance"] = dict(frame=t, cell=p, old=o, new=n)
log(f"Region 1 advance: frame {t} -> {t + 1}, cell (row {p[0]}, col {p[1]}): {o} -> {n}")
t, [(p, o, n)] = first_change(1, lambda o, n: o not in (4, 5) and (o <= 4) != (n <= 4))
ex["region1_jump"] = dict(frame=t, cell=p, old=o, new=n)
log(f"Region 1 jump:    frame {t} -> {t + 1}, cell (row {p[0]}, col {p[1]}): {o} -> {n}")
t, [(p, o, n)] = first_change(2, lambda o, n: n > o, start=0)   # rises are rare after the first ~200 steps
ex["region2_rise"] = dict(frame=t, cell=p, old=o, new=n)
log(f"Region 2 rise:    frame {t} -> {t + 1}, cell (row {p[0]}, col {p[1]}): {o} -> {n}  (allowed rises from {o}: {sorted({(o * k) % 10 for k in range(10)} - set(range(o + 1)))})")
t, [(p, o, n)] = first_change(2, lambda o, n: n < o)
ex["region2_drop"] = dict(frame=t, cell=p, old=o, new=n)
log(f"Region 2 drop:    frame {t} -> {t + 1}, cell (row {p[0]}, col {p[1]}): {o} -> {n}")
t, [(p, o, n)] = first_change(3)
ex["region3_redraw"] = dict(frame=t, cell=p, old=o, new=n)
log(f"Region 3 redraw:  frame {t} -> {t + 1}, cell (row {p[0]}, col {p[1]}): {o} -> {n}")
t, vals = first_change(4, k=2)
(p1, a, _), (p2, b, _) = vals
g = h[t]; g2 = g.copy(); g2[p1], g2[p2] = b, a
b1, b2 = count(g, *p1, a), count(g, *p2, b)
a1, a2 = count(g2, *p1, b), count(g2, *p2, a)
ex["region4_swap"] = dict(frame=t, cell1=p1, cell2=p2, value1=a, value2=b, before=(b1, b2), after=(a1, a2), gain=a1 + a2 - b1 - b2,
                          nbhd1_before=nbhd(g, *p1), nbhd2_before=nbhd(g, *p2), nbhd1_after=nbhd(g2, *p1), nbhd2_after=nbhd(g2, *p2))
log(f"Region 4 swap:    frame {t} -> {t + 1}, cells (row {p1[0]}, col {p1[1]}) = {a} and (row {p2[0]}, col {p2[1]}) = {b}")
for lab, rows in [("cell 1 before", ex["region4_swap"]["nbhd1_before"]), ("cell 2 before", ex["region4_swap"]["nbhd2_before"]),
                  ("cell 1 after", ex["region4_swap"]["nbhd1_after"]), ("cell 2 after", ex["region4_swap"]["nbhd2_after"])]:
    log(f"   {lab}:"); [log("      " + r) for r in rows]
log.claim("C34", f"Worked examples from held-out run 1 (frames counted from 0): Region 1 advance {ex['region1_advance']['old']}->{ex['region1_advance']['new']}, "
          f"jump {ex['region1_jump']['old']}->{ex['region1_jump']['new']}; Region 2 rise {ex['region2_rise']['old']}->{ex['region2_rise']['new']}, "
          f"drop {ex['region2_drop']['old']}->{ex['region2_drop']['new']}; Region 3 redraw {ex['region3_redraw']['old']}->{ex['region3_redraw']['new']}; "
          f"Region 4 swap of {a} and {b}: matching neighbours {b1}+{b2}={b1 + b2} before, {a1}+{a2}={a1 + a2} after, gain {a1 + a2 - b1 - b2}.",
          **ex)
log.save()
