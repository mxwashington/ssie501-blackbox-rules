"""07 - Rebuild the box from the rules and compare it with the held-out runs it never saw.

Claims proven here: C28-C30 (see REPORT.md). Needs metrics from 02 and 03 (run them first).
Run:  python analysis/07_simulate.py      (about 1-2 minutes)
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_discovery, load_heldout
from blackbox.regions import MASK
from blackbox.simulator import fit_parameters, simulate
from blackbox.report import Log, FIGURES, RESULTS
from blackbox.plotting import plt, AREA_COLOR, INK, MUTED

log = Log("07_simulate")
D, H = load_discovery(), load_heldout()
met2 = json.loads((RESULTS / "metrics" / "02_area3_redraw.json").read_text())
met3 = json.loads((RESULTS / "metrics" / "03_area4_swaps.json").read_text())
P = fit_parameters(D, met2["C07"]["p"], met3["C14"]["discovery"])
log(f"   fitted on discovery: Area 1 jump on first move (low, high branch) {np.round(P['jump_first'], 2).tolist()}, exit of 4 / 5 {P['exit45'][4]:.3f} / {P['exit45'][5]:.3f}; "
    f"Area 2 P(change|picked) {np.round(P['p_change2'], 2).tolist()}")

rng = np.random.default_rng(501)
S = [simulate(h[0], len(h) - 1, P, rng) for h in H]   # same starting grids as the held-out runs


def H_bits(v):
    c = np.bincount(v.ravel(), minlength=14); p = c[c > 0] / c.sum()
    return float(-(p * np.log2(p)).sum())


def summary(runs):
    out = {}
    for a in (1, 2, 3, 4):
        m = MASK[a]
        act = np.mean([((r[1:] != r[:-1])[:, m]).sum(1) for r in runs], 0)
        ent = np.mean([[H_bits(r[t][m]) for t in range(0, len(r), 50)] for r in runs], 0)
        out[a] = dict(activity=act, entropy=ent,
                      changes_total=float(np.mean([((r[1:] != r[:-1])[:, m]).sum() for r in runs])),
                      end_entropy=float(ent[-1]))
    out["area1_end_share_4"] = float(np.mean([(r[-1][MASK[1]] == 4).mean() for r in runs]))
    out["area1_end_share_45"] = float(np.mean([np.isin(r[-1][MASK[1]], [4, 5]).mean() for r in runs]))
    out["area2_end_share_0"] = float(np.mean([(r[-1][MASK[2]] == 0).mean() for r in runs]))
    return out


real, sim = summary(H), summary(S)
rows = {}
for a in (1, 2, 3, 4):
    rows[a] = dict(real_changes=real[a]["changes_total"], sim_changes=sim[a]["changes_total"],
                   real_end_entropy=real[a]["end_entropy"], sim_end_entropy=sim[a]["end_entropy"],
                   entropy_curve_max_gap=float(np.abs(real[a]["entropy"] - sim[a]["entropy"]).max()))
for a, r in rows.items():
    log(f"   Area {a}: cell changes per run real {r['real_changes']:.0f} vs simulated {r['sim_changes']:.0f}; "
        f"final entropy real {r['real_end_entropy']:.2f} vs simulated {r['sim_end_entropy']:.2f}; "
        f"largest entropy-curve gap {r['entropy_curve_max_gap']:.2f} bits")
log.claim("C28", "A box rebuilt only from these rules, started from the 30 held-out starting grids, closely reproduces Areas 2-4 and "
          "roughly reproduces Area 1: "
          "cell changes per run (real vs simulated) "
          + "; ".join(f"Area {a} {r['real_changes']:.0f} vs {r['sim_changes']:.0f}" for a, r in rows.items())
          + "; final entropy " + "; ".join(f"Area {a} {r['real_end_entropy']:.2f} vs {r['sim_end_entropy']:.2f}" for a, r in rows.items())
          + f"; largest gap between the real and simulated entropy curves {max(r['entropy_curve_max_gap'] for r in rows.values()):.2f} bits.",
          areas=rows)
log.claim("C29", f"Final states: Area 1 ends {real['area1_end_share_4']:.0%} 4s (simulated {sim['area1_end_share_4']:.0%}) and "
          f"{real['area1_end_share_45']:.0%} 4-or-5 (simulated {sim['area1_end_share_45']:.0%}); Area 2 ends "
          f"{real['area2_end_share_0']:.0%} zeros (simulated {sim['area2_end_share_0']:.0%}). The rebuilt Area 1 leaves more cells "
          f"mid-branch, so the timing of Area 1 jumps is the least-understood part of the box.",
          real={k: v for k, v in real.items() if isinstance(k, str)}, sim={k: v for k, v in sim.items() if isinstance(k, str)})

# ---------------------------------------------------------------- C30 ablation: same jump rates on both branches
P_sym = fit_parameters(D, met2["C07"]["p"], met3["C14"]["discovery"], split_branches=False)
rng2 = np.random.default_rng(501)
S_sym = [simulate(h[0], len(h) - 1, P_sym, rng2) for h in H]
sym_share4 = float(np.mean([(r[-1][MASK[1]] == 4).mean() for r in S_sym]))
log.claim("C30", f"The branch asymmetry (C18b) is what makes Area 1 end mostly 4s: with branch-specific jump and exit rates the rebuilt box "
          f"ends {sim['area1_end_share_4']:.0%} 4s; with one shared rate it ends {sym_share4:.0%} 4s (real: {real['area1_end_share_4']:.0%}).",
          split=sim["area1_end_share_4"], shared=sym_share4, real=real["area1_end_share_4"])

fig, axs = plt.subplots(2, 4, figsize=(13, 5.4), sharex=True)
w = 50
for i, a in enumerate((1, 2, 3, 4)):
    ax = axs[0, i]
    ra = real[a]["activity"]; sa = sim[a]["activity"]
    x = np.arange(0, len(ra) - w + 1, w) + w / 2
    ax.plot(x, [ra[j:j + w].mean() for j in range(0, len(ra) - w + 1, w)], color=AREA_COLOR[a], label="real (held-out)")
    ax.plot(x, [sa[j:j + w].mean() for j in range(0, len(sa) - w + 1, w)], color=INK, linestyle="--", linewidth=1.4, label="rebuilt from rules")
    ax.set_title(f"Area {a}: cells changed per step"); ax.set_ylim(bottom=0)
    if i == 0:
        ax.legend(fontsize=8)
    ax = axs[1, i]
    xe = np.arange(0, len(H[0]), 50)
    ax.plot(xe, real[a]["entropy"], color=AREA_COLOR[a])
    ax.plot(xe, sim[a]["entropy"], color=INK, linestyle="--", linewidth=1.4)
    ax.set_title(f"Area {a}: entropy (bits)"); ax.set_ylim(0, 3.5); ax.set_xlabel("step")
fig.suptitle("Real held-out runs vs. a box rebuilt only from the discovered rules (30 runs each, same starting grids)",
             fontweight="bold", fontsize=11)
fig.tight_layout(); fig.savefig(FIGURES / "07_simulation_vs_real.png"); plt.close(fig)
log.save()
