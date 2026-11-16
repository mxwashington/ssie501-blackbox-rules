"""08 - Build the printed tables of the written model and test the WRITTEN model on held-out runs.

Claims proven here: C31-C32 (see REPORT.md and the model specification document).
Writes results/model_tables.json = every number printed in the model specification.
Run:  python analysis/08_model_spec.py      (about 2 minutes)
"""
import sys, json, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_discovery, load_heldout
from blackbox.regions import MASK
from blackbox.spec_model import run, wait_bin
from blackbox.report import Log, FIGURES, RESULTS
from blackbox.plotting import plt, AREA_COLOR, INK

log = Log("08_model_spec")
D, H = load_discovery(), load_heldout()
BINS = ["first", "0-49", "50-99", "100-149", "150-199", "200+"]


def round_to_sum(probs, decimals):
    """Round a probability row so the printed numbers add up to exactly 1 (largest remainder)."""
    unit = 10 ** decimals
    raw = np.array(probs) * unit
    base = np.floor(raw).astype(int)
    for i in np.argsort(-(raw - base))[: unit - base.sum()]:
        base[i] += 1
    return (base / unit).tolist()


# ---------------------------------------------------------------- Region 1 tables (discovery runs)
m1 = MASK[1]
n = collections.Counter(); j = collections.Counter()
pk = {4: 0.0, 5: 0.0}; mv = {4: 0, 5: 0}
for a in D:
    last = {}
    v_all = a[:-1][:, m1]; ch_all = (a[1:] != a[:-1])[:, m1]
    for k in (4, 5):
        pk[k] += (v_all == k).sum() / 100; mv[k] += (ch_all & (v_all == k)).sum()
    for t in range(len(a) - 1):
        for r, c in np.argwhere((a[t + 1] != a[t]) & m1):
            o, nn = int(a[t, r, c]), int(a[t + 1, r, c])
            if o not in (4, 5):
                br = "low" if o <= 4 else "high"
                b = wait_bin(None if (r, c) not in last else t - last[(r, c)])
                n[(br, b)] += 1; j[(br, b)] += (o <= 4) != (nn <= 4)
            last[(r, c)] = t
region1 = {"jump": {br: {b: round(j[(br, b)] / n[(br, b)], 2) for b in BINS} for br in ("low", "high")},
           "jump_counts": {br: {b: n[(br, b)] for b in BINS} for br in ("low", "high")},
           "exit": {str(k): round(mv[k] / pk[k], 3) for k in (4, 5)}}

# ---------------------------------------------------------------- Region 2 tables (discovery runs)
m2 = MASK[2]
occ = np.zeros(10); T2 = np.zeros((10, 10))
for a in D:
    v, nn = a[:-1][:, m2], a[1:][:, m2]
    occ += np.bincount(v.ravel(), minlength=10)[:10]
    c = v != nn
    np.add.at(T2, (v[c], nn[c]), 1)
region2 = {"p_change": {str(o): round(min(1.0, 100 * T2[o].sum() / occ[o]), 2) for o in range(1, 10)},
           "target": {str(o): dict(zip([str(k) for k in range(10)], round_to_sum(T2[o] / T2[o].sum(), 2))) for o in range(1, 10)},
           "changes_seen": {str(o): int(T2[o].sum()) for o in range(1, 10)}}

# ---------------------------------------------------------------- Region 3 and 4 tables (from scripts 02 and 03)
met2 = json.loads((RESULTS / "metrics" / "02_area3_redraw.json").read_text())
met3 = json.loads((RESULTS / "metrics" / "03_area4_swaps.json").read_text())
p = np.array(met2["C07"]["p"])[:12]
p3 = round_to_sum(p / p.sum(), 4)
region3 = {"p": {str(k): p3[k] for k in range(12)}}
acc = met3["C14"]["discovery"]
region4 = {"accept": {str(k + 1): round(min(1.0, acc[k]), 2) for k in range(5)}}

TABLES = {"region1": region1, "region2": region2, "region3": region3, "region4": region4,
          "fitted_on": "discovery runs (Group 1, 61 runs x 2,000 steps)"}
(RESULTS / "model_tables.json").write_text(json.dumps(TABLES, indent=1))
log.claim("C31", "The printed model tables (all fitted on discovery runs only, saved to results/model_tables.json): "
          f"Region 1 jump rates low branch {region1['jump']['low']}, high branch {region1['jump']['high']}, exit rates {region1['exit']}; "
          f"Region 2 change rates {region2['p_change']}; Region 3 draw table {region3['p']}; Region 4 swap rates {region4['accept']}.",
          **TABLES)
for o in range(1, 10):
    log(f"   Region 2 target row for old value {o}: {region2['target'][str(o)]}  (from {region2['changes_seen'][str(o)]} changes)")


# ---------------------------------------------------------------- C32 the WRITTEN model vs held-out runs
def H_bits(v):
    cnt = np.bincount(np.asarray(v).ravel(), minlength=14); q = cnt[cnt > 0] / cnt.sum()
    return float(-(q * np.log2(q)).sum())


rng = np.random.default_rng(2026)
S = [run(h[0], len(h) - 1, TABLES, rng) for h in H]


def stats(runs):
    out = {}
    for a_ in (1, 2, 3, 4):
        m = MASK[a_]
        out[f"changes_area{a_}"] = float(np.mean([((r[1:] != r[:-1])[:, m]).sum() for r in runs]))
        out[f"end_entropy_area{a_}"] = float(np.mean([H_bits(r[-1][m]) for r in runs]))
        out[f"activity_area{a_}"] = np.mean([((r[1:] != r[:-1])[:, m]).sum(1) for r in runs], 0)
    out["area1_end_4"] = float(np.mean([(r[-1][MASK[1]] == 4).mean() for r in runs]))
    out["area1_end_45"] = float(np.mean([np.isin(r[-1][MASK[1]], [4, 5]).mean() for r in runs]))
    out["area2_end_0"] = float(np.mean([(r[-1][MASK[2]] == 0).mean() for r in runs]))
    out["area3_ge10_share"] = float(np.mean([(r[:, MASK[3]] >= 10).mean() for r in runs]))
    return out


real, spec = stats(H), stats(S)
keys = [k for k in real if not k.startswith("activity")]
rows = {k: (real[k], spec[k]) for k in keys}
for k, (rv, sv) in rows.items():
    log(f"   {k:20s} real {rv:9.4f}   written model {sv:9.4f}")
log.claim("C32", "The WRITTEN model (only the printed, rounded tables, implemented line by line from the specification in "
          "blackbox/spec_model.py), started from the 30 held-out starting grids: cell changes per run real vs model = "
          + "; ".join(f"Region {a_} {real[f'changes_area{a_}']:.0f} vs {spec[f'changes_area{a_}']:.0f}" for a_ in (1, 2, 3, 4))
          + "; final entropy (bits) = "
          + "; ".join(f"Region {a_} {real[f'end_entropy_area{a_}']:.2f} vs {spec[f'end_entropy_area{a_}']:.2f}" for a_ in (1, 2, 3, 4))
          + f"; Region 1 ends {real['area1_end_4']:.0%} 4s vs {spec['area1_end_4']:.0%}; Region 2 ends {real['area2_end_0']:.0%} zeros vs {spec['area2_end_0']:.0%}.",
          **{k: {"real": v[0], "model": v[1]} for k, v in rows.items()})

# figure: activity per region, real vs written model
fig, axs = plt.subplots(1, 4, figsize=(13, 3.0), sharex=True)
w = 50
for i, a_ in enumerate((1, 2, 3, 4)):
    ax = axs[i]
    ra, sa = real[f"activity_area{a_}"], spec[f"activity_area{a_}"]
    x = np.arange(0, len(ra) - w + 1, w) + w / 2
    ax.plot(x, [ra[k:k + w].mean() for k in range(0, len(ra) - w + 1, w)], color=AREA_COLOR[a_], label="real (held-out)")
    ax.plot(x, [sa[k:k + w].mean() for k in range(0, len(sa) - w + 1, w)], color=INK, linestyle="--", linewidth=1.4, label="written model")
    ax.set_title(f"Region {a_}: cells changed per step"); ax.set_ylim(bottom=0); ax.set_xlabel("step")
    if i == 0:
        ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIGURES / "08_written_model_vs_real.png"); plt.close(fig)
log.save()
