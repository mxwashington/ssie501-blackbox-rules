"""06 - Why entropy falls in some areas and not others (answers Part 1 feedback "Entropy falls - not always, but why?").

Claims proven here: C26-C27 (see REPORT.md). Needs results/metrics/02_area3_redraw.json (run 02 first).
Run:  python analysis/06_entropy.py
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_discovery, load_heldout
from blackbox.regions import MASK
from blackbox.report import Log, FIGURES, RESULTS
from blackbox.plotting import plt, AREA_COLOR, INK, MUTED

log = Log("06_entropy")
D, H = load_discovery(), load_heldout()


def H_bits(vals):
    c = np.bincount(vals.ravel(), minlength=14)
    p = c[c > 0] / c.sum()
    return float(-(p * np.log2(p)).sum())


def curves(runs):
    T = min(len(r) for r in runs)
    out = {a: np.array([[H_bits(r[t][MASK[a]]) for t in range(T)] for r in runs]) for a in (1, 2, 3, 4)}
    out["grid"] = np.array([[H_bits(r[t]) for t in range(T)] for r in runs])
    return out


cD, cH = curves(D), curves(H)
summ = {}
for k in (1, 2, 3, 4, "grid"):
    summ[str(k)] = dict(start=float(cH[k][:, 0].mean()), end=float(cH[k][:, -1].mean()),
                        area4_max_change=float(np.abs(cH[k] - cH[k][:, :1]).max()) if k == 4 else None)
p = np.array(json.loads((RESULTS / "metrics" / "02_area3_redraw.json").read_text())["C07"]["p"])
pp = p[p > 0]
H_p = float(-(pp * np.log2(pp)).sum())
log.claim("C26", "Per-area entropy, held-out runs (bits, start -> last step): "
          + "; ".join(f"Area {k} {v['start']:.2f} -> {v['end']:.2f}" for k, v in summ.items() if k != "grid")
          + f"; whole grid {summ['grid']['start']:.2f} -> {summ['grid']['end']:.2f}. Area 4's entropy never moves "
          f"(largest change {summ['4']['area4_max_change']:.3f} bits).", **summ)
log.claim("C27", "Each area's entropy follows from its rule: Areas 1 and 2 fall because 4/5 and 0 trap cells (C17, C21); Area 4 is fixed "
          f"because swaps keep the value counts (C11); Area 3 settles at the entropy of its redraw distribution, H(p) = {H_p:.2f} bits "
          f"(C07), close to its observed level of {summ['3']['end']:.2f} bits (entropy measured on only 100 cells reads slightly low).", H_p=H_p, area3_end=summ["3"]["end"])

fig, ax = plt.subplots(figsize=(7, 3.8))
t = np.arange(cH[1].shape[1])
for a in (1, 2, 3, 4):
    y = cH[a].mean(0)
    ax.plot(t, y, color=AREA_COLOR[a], label=f"Area {a}")
    ax.text(t[-1] + 20, y[-1], f"Area {a}", color=INK, va="center", fontsize=9)
ax.axhline(H_p, color=MUTED, linestyle=":", linewidth=1.2)
ax.text(1050, H_p + 0.12, f"H(p) of Area 3's redraw distribution = {H_p:.2f}", color=MUTED, fontsize=8)
ax.set_xlim(0, t[-1] + 260); ax.set_ylim(0, 3.5)
ax.set_xlabel("step"); ax.set_ylabel("Shannon entropy (bits)")
ax.set_title("Entropy by area (mean of 30 held-out runs)")
ax.legend(loc="upper right", bbox_to_anchor=(0.88, 0.86), fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(FIGURES / "06_entropy_by_area.png"); plt.close(fig)
log.save()
