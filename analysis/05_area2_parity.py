"""05 - Area 2: 0 is a dead end, drops can go anywhere lower, rises obey an odd/even (mod 10) rule.

Claims proven here: C21-C25 (see REPORT.md).
Run:  python analysis/05_area2_parity.py
"""
import sys, collections, math, random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_discovery, load_heldout, load_group2
from blackbox.regions import MASK, MOORE, neighbour_values
from blackbox.events import single_events
from blackbox.report import Log, FIGURES
from blackbox.plotting import plt, AREA_COLOR, INK, MUTED

log = Log("05_area2_parity")
D, H, G2 = load_discovery(), load_heldout(), load_group2()
m = MASK[2]
REACH = {o: {(o * k) % 10 for k in range(10)} for o in range(10)}      # values o*k mod 10 = multiples of gcd(o,10)
ALLOWED = {(o, t) for o in range(10) for t in range(10) if t != o and (t < o or t in REACH[o])}

# ---------------------------------------------------------------- C21 / C22 / C23 allowed moves
seen = {}
for name, runs in [("discovery", D), ("heldout", H), ("group2", G2)]:
    seen[name] = collections.Counter((e["old"], e["new"]) for e in single_events(runs, 2))
from0 = {k: sum(v for (o, n), v in c.items() if o == 0) for k, c in seen.items()}
log.claim("C21", f"0 is a dead end in Area 2: changes away from 0 = {from0}.", **from0)

rises_pred = {(o, t) for (o, t) in ALLOWED if t > o}
rises_seen = {k for k in seen["discovery"] if k[1] > k[0]}
viol = {k: sum(v for pair, v in c.items() if pair not in ALLOWED) for k, c in seen.items()}
log.claim("C22", "A value can only rise to a value it reaches by multiplying by a whole number mod 10 (= a multiple of gcd(value, 10)): "
          "even values only rise to even values, 5 never rises, 1/3/7 can rise to anything. On discovery data the set of rises seen "
          f"equals the predicted set exactly (missing: {sorted(rises_pred - rises_seen)}, extra: {sorted(rises_seen - rises_pred)}); "
          f"changes breaking the rule: {viol}.",
          predicted_rises=sorted(rises_pred), seen_rises=sorted(rises_seen), violations=viol)
for o in range(1, 10):
    log(f"   from {o}: rises allowed to {sorted(t for t in REACH[o] if t > o)}; rises seen {sorted(t for (a, t) in rises_seen if a == o)}")

drops_missing = sorted((o, t) for o in range(1, 10) for t in range(o) if seen["discovery"][(o, t)] == 0)
log.claim("C23", f"Drops are unrestricted: every drop to a lower value occurs (drops never seen on discovery data: {drops_missing}).",
          drops_never_seen=drops_missing)

# ---------------------------------------------------------------- C24 how often each value changes
def change_prob(runs):
    occ = collections.Counter(); ch = collections.Counter(); z = collections.Counter()
    for a in runs:
        v = a[:-1][:, m]; n = a[1:][:, m]; c = v != n
        for k in range(10):
            occ[k] += int((v == k).sum()); ch[k] += int((c & (v == k)).sum()); z[k] += int((c & (v == k) & (n == 0)).sum())
    return {k: 100 * ch[k] / occ[k] for k in range(1, 10)}, {k: 100 * z[k] / occ[k] for k in range(1, 10)}
pcD, p0D = change_prob(D); pcH, p0H = change_prob(H)
log.claim("C24", "When picked, bigger values generally change more often (estimated as 100 x per-step change rate, since one of 100 cells is picked per step - C02): P(change | picked) for values 1..9 = "
          + ", ".join(f"{pcD[k]:.2f}" for k in range(1, 10)) + " (discovery) and "
          + ", ".join(f"{pcH[k]:.2f}" for k in range(1, 10)) + " (held-out).",
          discovery=pcD, heldout=pcH, drop_to_zero_discovery=p0D, drop_to_zero_heldout=p0H)

# ---------------------------------------------------------------- C25 neighbours: more 0s nearby -> more drops to 0, but no copying
c = collections.Counter(); z = collections.Counter()
for e in single_events(D, 2, MOORE, wrap=False):
    k = sum(x == 0 for x in e["nbrs"]); c[k] += 1; z[k] += e["new"] == 0
tab = {k: z[k] / c[k] for k in sorted(c) if c[k] >= 300}
random.seed(0)
cells = np.argwhere(m); hit = ctl = n = 0
for e in single_events(D, 2, MOORE, wrap=False):
    g = D[e["run"]][e["t"]]
    hit += e["new"] in e["nbrs"]
    rr, cc = cells[random.randrange(len(cells))]
    ctl += e["new"] in neighbour_values(g, rr, cc); n += 1
log.claim("C25", "Neighbours are associated with the outcome but do not decide it: a changing cell drops to 0 more often when more of its 8 neighbours are 0 ("
          + ", ".join(f"{k} zero nbrs: {v:.0%}" for k, v in tab.items())
          + f"), but it is not copying a neighbour: its new value appears among its own neighbours {hit / n:.0%} of the time vs "
          f"{ctl / n:.0%} among a random Area 2 cell's neighbours.",
          p_drop_to_0_by_zero_neighbours=tab, new_in_own_nbrs=hit / n, new_in_random_cell_nbrs=ctl / n)

# ---------------------------------------------------------------- figure: allowed-move matrix with counts
C = seen["discovery"]
M = np.full((10, 10), np.nan)
for o in range(10):
    for t in range(10):
        if (o, t) in ALLOWED:
            M[o, t] = C[(o, t)]
fig, ax = plt.subplots(figsize=(5.6, 5))
allowed = np.zeros((10, 10))
for (o, t) in ALLOWED:
    allowed[o, t] = 1
ax.imshow(allowed, cmap="Oranges", vmin=0, vmax=3)
for o in range(10):
    for t in range(10):
        if o == t:
            continue
        txt = str(C[(o, t)]) if (o, t) in ALLOWED else "x"
        ax.text(t, o, txt, ha="center", va="center", fontsize=7, color=INK if (o, t) in ALLOWED else MUTED)
ax.set_xticks(range(10)); ax.set_yticks(range(10)); ax.grid(False)
ax.set_xlabel("new value"); ax.set_ylabel("old value")
ax.set_title("Area 2: shaded = allowed by the rule\nnumbers = changes seen (discovery); x = never seen")
fig.tight_layout(); fig.savefig(FIGURES / "05_area2_moves.png"); plt.close(fig)
log.save()
