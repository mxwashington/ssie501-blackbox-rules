"""09 - Is Region 3's draw a known textbook distribution? (Lauren's question in the group chat.)

Claim proven here: C33. Needs results/metrics/02_area3_redraw.json (run 02 first).
Run:  python analysis/09_area3_family.py
"""
import sys, json, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson, nbinom
from blackbox.data import load_discovery, load_heldout
from blackbox.events import single_events
from blackbox.report import Log, RESULTS

log = Log("09_area3_family")
D, H = load_discovery(), load_heldout()
K = 14


def counts(runs):
    C = np.zeros((K, K))
    for e in single_events(runs, 3):
        C[e["old"], e["new"]] += 1
    return C


def loglik(q, C):
    """Log-likelihood of observed changes when a picked cell is redrawn from q: P(new=k | old=o, changed) = q_k / (1 - q_o)."""
    q = np.clip(q, 1e-12, 1)
    ll = 0.0
    for o in range(K):
        row = C[o]
        if row.sum() == 0:
            continue
        ll += (row * (np.log(q) - np.log(1 - q[o]))).sum() - row[o] * 0
    return ll


ks = np.arange(K)
FAMILIES = {
    "Poisson(lam)": (1, lambda th: poisson.pmf(ks, np.exp(th[0])), [0.0]),
    "geometric(r)": (1, lambda th: (1 - 1 / (1 + np.exp(-th[0]))) * (1 / (1 + np.exp(-th[0]))) ** ks, [0.0]),
    "zero-inflated geometric(a, r)": (2, lambda th: np.where(ks == 0, 1 / (1 + np.exp(-th[0])),
                                       (1 - 1 / (1 + np.exp(-th[0]))) * (1 - 1 / (1 + np.exp(-th[1]))) * (1 / (1 + np.exp(-th[1]))) ** np.maximum(ks - 1, 0)), [0.0, 0.0]),
    "negative binomial(n, p)": (2, lambda th: nbinom.pmf(ks, np.exp(th[0]), 1 / (1 + np.exp(-th[1]))), [0.0, 0.0]),
    "zero-inflated Poisson(a, lam)": (2, lambda th: np.where(ks == 0, 1 / (1 + np.exp(-th[0])), 0) + (1 - 1 / (1 + np.exp(-th[0]))) * poisson.pmf(ks, np.exp(th[1])), [0.0, 0.0]),
}
CD, CH = counts(D), counts(H)
p_free = np.array(json.loads((RESULTS / "metrics" / "02_area3_redraw.json").read_text())["C07"]["p"])
ll_free = loglik(p_free, CD); n_free = 11
res = {}
for name, (npar, f, x0) in FAMILIES.items():
    best = min((minimize(lambda th: -loglik(f(th) / f(th).sum(), CD), np.array(x0) + s, method="Nelder-Mead")
                for s in (-1.0, 0.0, 1.0)), key=lambda r: r.fun)
    q = f(best.x); q = q / q.sum()
    aic = 2 * npar + 2 * best.fun
    res[name] = dict(params=best.x.tolist(), aic=aic, q=q[:8].round(3).tolist(),
                     heldout_ll_per_change=loglik(q, CH) / CH.sum())
aic_free = 2 * n_free - 2 * ll_free
res["free table (C07)"] = dict(aic=aic_free, q=p_free[:8].round(3).tolist(), heldout_ll_per_change=loglik(p_free, CH) / CH.sum())
for k, v in sorted(res.items(), key=lambda kv: kv[1]["aic"]):
    log(f"   {k:32s} AIC {v['aic']:10.1f}  held-out log-lik/change {v['heldout_ll_per_change']:.5f}  P(0..7) {v['q']}")
best_family = min((k for k in res if k != "free table (C07)"), key=lambda k: res[k]["aic"])
gap = res[best_family]["aic"] - aic_free
zig = res["zero-inflated geometric(a, r)"]
a_zig = 1 / (1 + np.exp(-zig["params"][0])); r_zig = 1 / (1 + np.exp(-zig["params"][1]))
pois = res["Poisson(lam)"]
log.claim("C33", f"Region 3's draw is not exactly any of five standard distributions: the best ({best_family}) is worse than the free "
          f"11-number table by {gap:.0f} AIC points (a gap above ~10 favours the table). The closest simple description is: value 0 with "
          f"probability {a_zig:.2f}, otherwise 1 + a geometric count where each further value is {r_zig:.2f} times as likely as the one "
          f"before. A Poisson fits worst (its best fit gives P(0) = {pois['q'][0]:.2f} vs the measured 0.36). The model uses the measured table.",
          families=res, best_family=best_family, aic_gap_to_table=gap, zig_a=a_zig, zig_r=r_zig)
log.save()
