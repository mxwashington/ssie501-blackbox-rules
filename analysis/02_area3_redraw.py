"""02 - Area 3: the picked cell gets a fresh random value from one fixed distribution.

Claims proven here: C06-C10 (see REPORT.md).
Run:  python analysis/02_area3_redraw.py
"""
import sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from scipy.stats import chi2_contingency
from blackbox.data import load_discovery, load_heldout
from blackbox.regions import MASK, MOORE
from blackbox.events import single_events
from blackbox.report import Log, FIGURES
from blackbox.plotting import plt, AREA_COLOR, INK, MUTED

log = Log("02_area3_redraw")
D, H = load_discovery(), load_heldout()
NV = 14  # values seen: 0..13


def transition_matrix(runs):
    """P(old -> new in one step) / (1/100): approximately P(new | cell picked) because one of 100 cells is picked per step."""
    m = MASK[3]
    occ = np.zeros(NV); T = np.zeros((NV, NV))
    for a in runs:
        v, n = a[:-1][:, m], a[1:][:, m]
        occ += np.bincount(v.ravel(), minlength=NV)[:NV]
        c = v != n
        np.add.at(T, (v[c], n[c]), 1)
    return 100 * T / np.maximum(occ, 1)[:, None], occ


def fit_draw_distribution(runs):
    """Maximum-likelihood p_k for 'picked cell redrawn from p': P(new=k | old=o, changed) = p_k / (1 - p_o)."""
    by_old = collections.defaultdict(collections.Counter)
    for e in single_events(runs, 3):
        by_old[e["old"]][e["new"]] += 1
    p = np.ones(NV) / NV
    for _ in range(500):
        num = np.zeros(NV); den = np.zeros(NV)
        for o, C in by_old.items():
            s = sum(C.values())
            for k, c in C.items():
                num[k] += c
            for k in range(NV):
                if k != o:
                    den[k] += s / (1 - p[o])
        p = np.where(den > 0, num / np.maximum(den, 1e-12), 0); p /= p.sum()
    return p, by_old


# ---------------------------------------------------------------- C06 new value ignores the old value
P, occ = transition_matrix(D)
log("   P(old -> new | picked) for old = 0..9 (rows), new = 0..6 (cols); diagonal = probability of redrawing the same value:")
for o in range(10):
    row = P[o].copy(); row[o] = np.nan
    log("   old", o, " ".join(f"{x:5.3f}" if not np.isnan(x) else "  -  " for x in row[:7]))
# For each target k, compare the rows (excluding the diagonal and rare rows): spread of P(-> k)
spreads = {}
for k in range(4):
    col = [P[o, k] for o in range(10) if o != k and occ[o] > 50000]
    spreads[k] = (min(col), max(col))
log.claim("C06", "In Area 3 the new value does not depend on the old value: the rows of the transition table are nearly identical. "
          + "; ".join(f"P(->{k}) ranges only {a:.2f}-{b:.2f} across old values" for k, (a, b) in spreads.items()) + ".",
          ranges=spreads)

# ---------------------------------------------------------------- C07 the fixed distribution
p, by_old = fit_draw_distribution(D)
ratios = [p[k + 1] / p[k] for k in range(1, 7)]
log.claim("C07", "The fixed redraw distribution (fitted on discovery data) is "
          + ", ".join(f"P({k})={p[k]:.3f}" for k in range(8))
          + f"; above 1 each value is about half as likely as the one before (ratios {', '.join(f'{r:.2f}' for r in ratios)}).",
          p=p, ratios=ratios)

# ---------------------------------------------------------------- C08 the new value ignores the neighbours
c = collections.Counter(); z = collections.Counter(); mean_new = collections.defaultdict(list)
for e in single_events(D, 3, MOORE, wrap=False):
    k = sum(x == 0 for x in e["nbrs"]); c[k] += 1; z[k] += e["new"] == 0
p0 = {k: z[k] / c[k] for k in sorted(c) if c[k] >= 500}
log.claim("C08", "The new value also ignores the neighbours: P(new = 0) stays "
          f"{min(p0.values()):.2f}-{max(p0.values()):.2f} whether 0 or 7 of the 8 neighbours are 0.",
          p_new0_by_zero_neighbours=p0)

# ---------------------------------------------------------------- C09 held-out test
def new_share(runs):
    cc = collections.Counter(e["new"] for e in single_events(runs, 3))
    s = sum(cc.values())
    return np.array([cc[k] for k in range(NV)]), s
cD, nD = new_share(D); cH, nH = new_share(H)
keep = np.where(cD + cH > 20)[0]
pval = chi2_contingency(np.vstack([cD[keep], cH[keep]]))[1]
# predicted number of steps with a change = sum over steps of (1/100) * sum over cells (1 - p[value])
pred = obs = 0.0
for a in H:
    v = a[:-1][:, MASK[3]]
    pred += (1 - p[np.clip(v, 0, NV - 1)]).sum() / 100
    obs += ((a[1:] != a[:-1])[:, MASK[3]]).any(1).sum()
log.claim("C09", f"Held-out test: the new-value distribution in the held-out runs matches discovery (chi-square p = {pval:.2f}), and "
          f"'one random cell redrawn from p each step' predicts {pred:,.0f} changing steps vs {obs:,.0f} observed "
          f"({(obs - pred) / pred:+.1%}).", chi2_p=pval, predicted=pred, observed=obs)

# ---------------------------------------------------------------- C10 values 10-13
seen = {}
for a_ in (1, 2, 3, 4):
    vals = np.concatenate([r[:, MASK[a_]].ravel() for r in D + H])
    seen[a_] = int((vals >= 10).sum())
share_big = sum(p[10:])
log.claim("C10", f"Values 10-13 appear only in Area 3 (cell-steps with a value >= 10 by area: {seen}). They are the rare "
          f"tail of the redraw distribution: P(new >= 10) = {share_big:.4f}.", by_area=seen, p_tail=share_big)

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ks = np.arange(10)
emp = cD[:10] / nD
ax.bar(ks - 0.2, p[:10], width=0.38, color=AREA_COLOR[3], label="fitted redraw distribution p (discovery)")
ax.bar(ks + 0.2, cH[:10] / nH, width=0.38, color="#9fd9c3", label="new values seen in held-out runs*")
ax.set_xticks(ks); ax.set_xlabel("new value given to the picked cell"); ax.set_ylabel("probability")
ax.set_title("Area 3: picked cell gets a fresh draw from one fixed distribution")
ax.legend(loc="upper right", fontsize=8)
ax.text(0.99, 0.55, "*counts only draws that changed the value,\n so it differs slightly from p by design",
        transform=ax.transAxes, ha="right", fontsize=7, color=MUTED)
fig.tight_layout(); fig.savefig(FIGURES / "02_area3_distribution.png"); plt.close(fig)
log.save()
