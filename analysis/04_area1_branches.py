"""04 - Area 1: two one-way branches (0->1->2->3->4 and 9->8->7->6->5) with jumps between them.

Claims proven here: C16-C20 (see REPORT.md).
Run:  python analysis/04_area1_branches.py
"""
import sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_discovery, load_heldout, load_group2
from blackbox.regions import MASK, VN, MOORE
from blackbox.events import single_events
from blackbox.report import Log, FIGURES
from blackbox.plotting import plt, AREA_COLOR, INK, MUTED

log = Log("04_area1_branches")
D, H, G2 = load_discovery(), load_heldout(), load_group2()
low = lambda v: v <= 4
ADVANCE = {(v, v + 1) for v in range(4)} | {(v, v - 1) for v in range(6, 10)}
JUMP = {(v, 9) for v in range(5)} | {(v, 0) for v in range(5, 10)}

# ---------------------------------------------------------------- C16 only two kinds of move
res = {}
for name, runs in [("discovery", D), ("heldout", H), ("group2", G2)]:
    c = collections.Counter()
    for e in single_events(runs, 1):
        k = (e["old"], e["new"])
        c["advance" if k in ADVANCE else "jump" if k in JUMP else "other"] += 1
    res[name] = dict(c)
log.claim("C16", "Every Area 1 change is either one step along its branch (0->1->2->3->4 or 9->8->7->6->5) or a jump to the START "
          "of the other branch (0-4 -> 9, 5-9 -> 0). Counts (advance / jump / anything else): "
          + "; ".join(f"{k} {v.get('advance', 0)}/{v.get('jump', 0)}/{v.get('other', 0)}" for k, v in res.items()) + ".",
          **res)

# ---------------------------------------------------------------- C17 non-end cells always move; 4 and 5 are near-absorbing
m = MASK[1]
pred = obs = picked45 = moved45 = 0
pk = {4: 0.0, 5: 0.0}; mv = {4: 0, 5: 0}
for a in H:
    v = a[:-1][:, m]; ch = (a[1:] != a[:-1])[:, m]
    obs += ch.any(1).sum()
    pred += (~np.isin(v, [4, 5])).sum() / 100
    picked45 += np.isin(v, [4, 5]).sum() / 100
    moved45 += (ch & np.isin(v, [4, 5])).sum()
    for k in (4, 5):
        pk[k] += (v == k).sum() / 100; mv[k] += (ch & (v == k)).sum()
exit_rate = moved45 / picked45
exit4, exit5 = mv[4] / pk[4], mv[5] / pk[5]
end_share = np.mean([np.isin(r[-1][m], [4, 5]).mean() for r in D + H])
four_share = np.mean([(r[-1][m] == 4).mean() for r in D + H])
log.claim("C17", f"A picked Area 1 cell that is not 4 or 5 always moves; a picked 4 or 5 leaves only {exit_rate:.1%} of the time "
          f"(held-out; a 4 leaves {exit4:.1%}, a 5 leaves {exit5:.1%}). 'Every non-4/5 pick moves' plus that exit rate predicts {pred + moved45:,.0f} changing steps vs "
          f"{obs:,.0f} observed. By the last recorded step {end_share:.0%} of Area 1 cells are 4 or 5 ({four_share:.0%} are 4).",
          predicted=pred + moved45, observed=obs, exit_rate_4_5=exit_rate, exit_4=exit4, exit_5=exit5, end_share_4_5=end_share, end_share_4=four_share)


# ---------------------------------------------------------------- C18 / C19 when does a cell jump instead of advance?
def branch_table(runs):
    """Jump rate split by branch (low 0-3 / high 6-9) and first vs later move of the cell."""
    n = collections.Counter(); j = collections.Counter()
    for a in runs:
        seen = set()
        for t in range(len(a) - 1):
            for r, c in np.argwhere((a[t + 1] != a[t]) & m):
                o, nn = int(a[t, r, c]), int(a[t + 1, r, c])
                if o not in (4, 5):
                    key = ("low branch" if low(o) else "high branch") + (", first move" if (r, c) not in seen else ", later moves")
                    n[key] += 1; j[key] += low(o) != low(nn)
                seen.add((r, c))
    return {k: (j[k] / n[k], n[k]) for k in sorted(n)}


def history_table(runs):
    first = collections.Counter(); after = collections.Counter(); age = collections.Counter()
    nf = collections.Counter(); na = collections.Counter(); nage = collections.Counter()
    for a in runs:
        last_t = {}; last_kind = {}
        for t in range(len(a) - 1):
            for r, c in np.argwhere((a[t + 1] != a[t]) & m):
                o, n = int(a[t, r, c]), int(a[t + 1, r, c])
                j = low(o) != low(n)
                key = (r, c)
                if o not in (4, 5):
                    if key not in last_kind:
                        nf["first move"] += 1; first["first move"] += j
                    else:
                        lab = "after a jump" if last_kind[key] else "after an advance"
                        na[lab] += 1; after[lab] += j
                        b = min((t - last_t[key]) // 25 * 25, 300)
                        nage[b] += 1; age[b] += j
                last_t[key] = t; last_kind[key] = j
    hist = {k: (first[k] / nf[k], nf[k]) for k in nf}
    hist.update({k: (after[k] / na[k], na[k]) for k in na})
    ages = {int(b): (age[b] / nage[b], nage[b]) for b in sorted(nage)}
    return hist, ages


histD, agesD = history_table(D)
histH, agesH = history_table(H)
alljD = np.mean([low(e["old"]) != low(e["new"]) for e in single_events(D, 1) if e["old"] not in (4, 5)])
log.claim("C18", f"Whether a moving cell jumps is random: {alljD:.0%} of moves by non-4/5 cells are jumps (discovery). It depends on the "
          "cell's own history: jump rate on a cell's first move after Reset vs later moves = "
          + ", ".join(f"{k} {v[0]:.0%}" for k, v in histD.items()) + " (discovery); "
          + ", ".join(f"{k} {v[0]:.0%}" for k, v in histH.items()) + " (held-out).",
          overall_jump_share=alljD, discovery=histD, heldout=histH)
lab = lambda b: "300+" if b >= 300 else f"{b}-{b + 24}"
brD, brH = branch_table(D), branch_table(H)
log.claim("C18b", "Jumps are much more common on the high branch (5-9) than the low branch (0-4) (C30 shows this is why Area 1 ends mostly 4s): "
          + ", ".join(f"{k} {v[0]:.0%}" for k, v in brD.items()) + " (discovery); "
          + ", ".join(f"{k} {v[0]:.0%}" for k, v in brH.items()) + " (held-out).", discovery=brD, heldout=brH)
log.claim("C19", "The longer since a cell last moved, the more likely its next move is a jump (time since last move -> jump rate): "
          + ", ".join(f"{lab(b)} steps {v[0]:.0%}" for b, v in agesD.items() if b in (0, 50, 100, 200, 300)) + " (discovery); "
          + ", ".join(f"{lab(b)} steps {v[0]:.0%}" for b, v in agesH.items() if b in (0, 50, 100, 200, 300)) + " (held-out).",
          discovery=agesD, heldout=agesH)

# ---------------------------------------------------------------- C20 neighbours are associated, but never decide
by_other = collections.defaultdict(lambda: [0, 0])
dir_rates = []
for e in single_events(D, 1, VN, wrap=False):
    o, n = e["old"], e["new"]
    if o in (4, 5):
        continue
    j = low(o) != low(n)
    k = sum(low(x) != low(o) for x in e["nbrs"])
    by_other[k][0] += 1; by_other[k][1] += j
tab = {k: (v[1] / v[0], v[0]) for k, v in sorted(by_other.items())}
for d_i, d in enumerate(MOORE):
    c = collections.Counter(); s = collections.Counter()
    for e in single_events(D, 1, [d], wrap=False):
        if e["old"] in (4, 5) or not e["nbrs"]:
            continue
        x = e["nbrs"][0]; c[x] += 1; s[x] += low(e["old"]) != low(e["new"])
    dir_rates += [s[x] / c[x] for x in c if c[x] >= 200]
log.claim("C20", "Neighbours are associated with jumping but never decide it: jump rate by number of side-neighbours on the other branch = "
          + ", ".join(f"{k}: {v[0]:.0%}" for k, v in tab.items())
          + f"; across all 8 directions and all neighbour values the jump rate stays between {min(dir_rates):.0%} and {max(dir_rates):.0%} "
          "(a deciding neighbour would push it to 0% or 100%).",
          by_other_branch_side_neighbours=tab, per_direction_min=min(dir_rates), per_direction_max=max(dir_rates))

# ---------------------------------------------------------------- figure
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
ax = axs[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 4)
ax.set_title("Area 1: allowed moves")
for i, v in enumerate([0, 1, 2, 3, 4]):
    ax.text(1 + 1.8 * i, 3, str(v), ha="center", va="center", fontsize=14, fontweight="bold", color=INK,
            bbox=dict(boxstyle="circle", fc="#cde2fb" if v != 4 else AREA_COLOR[1], ec="none"))
for i, v in enumerate([9, 8, 7, 6, 5]):
    ax.text(1 + 1.8 * i, 1, str(v), ha="center", va="center", fontsize=14, fontweight="bold", color=INK,
            bbox=dict(boxstyle="circle", fc="#cde2fb" if v != 5 else AREA_COLOR[1], ec="none"))
for i in range(4):
    ax.annotate("", xy=(1 + 1.8 * (i + 1) - 0.35, 3), xytext=(1 + 1.8 * i + 0.35, 3), arrowprops=dict(arrowstyle="->", color=INK))
    ax.annotate("", xy=(1 + 1.8 * (i + 1) - 0.35, 1), xytext=(1 + 1.8 * i + 0.35, 1), arrowprops=dict(arrowstyle="->", color=INK))
ax.text(5, 2.0, "jump: any 0-4 -> 9      any 5-9 -> 0\n(4 and 5 are near-dead ends)", ha="center", va="center", fontsize=9, color=MUTED)
ax = axs[1]
bs = [b for b in agesD if b <= 300]
ax.plot(bs, [agesD[b][0] for b in bs], color=AREA_COLOR[1], marker="o", markersize=4, label="discovery")
bh = [b for b in agesH if b <= 300]
ax.plot(bh, [agesH[b][0] for b in bh], color="#86b6ef", marker="o", markersize=4, label="held-out")
ax.set_xlabel("steps since the cell last moved"); ax.set_ylabel("share of moves that are jumps")
ax.set_title("Jumps get likelier the longer a cell sits"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIGURES / "04_area1.png"); plt.close(fig)
log.save()
