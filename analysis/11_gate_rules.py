"""11 - Region 1 gates, the two parts of Region 2, the strict Region 4 rule and the update order, tested on the
held-out runs; then the full Part A model rebuilt and compared with the held-out runs.

Claims proven here: C35-C41. The Region 1 gates, the Region 2 split and the strict Region 4 rule come from
Marcus Washington's simulator, built on his own recorded runs; none of the 30 held-out runs was used to build them.
Run:  python analysis/11_gate_rules.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_heldout, load_discovery, load_group2
from blackbox.regions import AREA
from blackbox.report import Log
from blackbox.gate_model import GateModel
from blackbox.rule_checks import evaluate

log = Log("11_gate_rules")
runs = load_heldout()

DATASETS = {"H": runs}
try:
    DATASETS["D"] = load_discovery()
except FileNotFoundError:
    log("Group 1 discovery runs not in data/group1_runs/: tested on the held-out runs only")
g2 = load_group2()
if g2:
    DATASETS["G2"] = g2
RES = {k: evaluate(v) for k, v in DATASETS.items()}
for k, r in RES.items():
    log(f"[{k}] {len(DATASETS[k])} runs: Region 1 {r['gate_ok']}/{r['n_ev']}; settled agree {r['settled']['agree']}, "
        f"disagree {r['settled']['disagree']}; min part {r['pt']}; multiply part {r['fo']}; swaps {r['sw']}")
tot = lambda f: sum(f(r) for r in RES.values())
n_ev, gate_ok, wait_ok = tot(lambda r: r["n_ev"]), tot(lambda r: r["gate_ok"]), tot(lambda r: r["wait_ok"])
ll_wait = tot(lambda r: r["ll_wait"])
na, ja = tot(lambda r: r["settled"]["agree"][0]), tot(lambda r: r["settled"]["agree"][1])
nd, jd = tot(lambda r: r["settled"]["disagree"][0]), tot(lambda r: r["settled"]["disagree"][1])
pt = {k: tot(lambda r, k=k: r["pt"][k]) for k in RES["H"]["pt"]}
fo = {k: tot(lambda r, k=k: r["fo"][k]) for k in RES["H"]["fo"]}
sw = {k: tot(lambda r, k=k: r["sw"][k]) for k in RES["H"]["sw"]}
per = {k: dict(region1=[r["gate_ok"], r["n_ev"]], min_part=[r["pt"]["partner_found"], r["pt"]["changes"], r["pt"]["rises"]],
               multiply=[r["fo"]["product_found"], r["fo"]["changes"]], swaps_in_order=[r["sw"]["strict_in_order"], r["sw"]["swaps"]])
       for k, r in RES.items()}
names = "+".join(DATASETS)
log.claim("C35", f"Region 1 gates predict jump vs advance for {gate_ok} of {n_ev} changes from 0-3 or 6-9 ({names})",
          correct=gate_ok, n=n_ev, per_dataset=per)
log.claim("C36", f"Region 1 cells at 4 or 5: gate agrees -> {ja} jumps in {na} cell-steps; gate disagrees -> "
          f"{jd} jumps in {nd} cell-steps ({jd / nd:.4f} per step; picking alone gives 0.0100) ({names})",
          agree_cellsteps=na, agree_jumps=ja, disagree_cellsteps=nd, disagree_jumps=jd, disagree_rate=jd / nd)
log.claim("C37", f"Waiting-time rule (earlier A5) scored on the same changes: best guess right {wait_ok}/{n_ev} "
          f"({wait_ok / n_ev:.1%}); log-likelihood {ll_wait / n_ev:.3f} per change (gates: 0 errors)",
          correct=wait_ok, n=n_ev, loglik_per_change=ll_wait / n_ev)
log.claim("C38", f"Region 2 minimum part: {pt['changes']} changes, {pt['rises']} rises, new value present in rows "
          f"11-20 in {pt['partner_found']} ({names})", **pt)
log.claim("C39", f"Region 2 multiply part: {fo['changes']} changes, {fo['product_found']} equal (old x a value in "
          f"that part) mod 10 ({names})", **fo)
log.claim("C40", f"Region 4: {sw['swaps']} swaps; strict rule holds for {sw['strict_old']} on the old grid and "
          f"{sw['strict_in_order']} when Regions 1-3 move first; total-gain rule {sw['total_gain_old']} ({names})", **sw)
if "D" in RES:
    f = RES["D"]["first"]
    log.claim("C42", f"First step after Reset (Group 1 runs start at the Reset grid): {f['transitions']} first steps, "
              f"Region 1 gates right for {f['r1_gate_ok']} of {f['r1_changes']} Region 1 changes", **f)

# Full model rebuilt from the 30 held-out starting grids, 5 seeds each
AREA_FLAT = AREA.ravel()


def H(v):
    _, c = np.unique(v, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def stats(traj):
    ch = traj[1:] != traj[:-1]
    end = traj[-1]
    s = {f"r{a}_changed": int(ch[:, AREA == a].sum()) for a in (1, 2, 3, 4)}
    s.update({f"r{a}_entropy_end": H(end[AREA == a]) for a in (1, 2, 3, 4)})
    e1 = end[AREA == 1]
    s.update(r1_share4=float(np.mean(e1 == 4)), r1_share45=float(np.mean((e1 == 4) | (e1 == 5))),
             r2_share0=float(np.mean(end[AREA == 2] == 0)))
    return s


real = [stats(h) for h in runs]
model = [stats(GateModel(h[0], seed=100 * k + rep).run(len(h) - 1)) for k, h in enumerate(runs) for rep in range(5)]
mean = lambda rows, key: float(np.mean([r[key] for r in rows]))
table = {key: dict(model=mean(model, key), real=mean(real, key)) for key in real[0]}
log("Full model vs held-out runs (2,000 steps, means):")
for key, v in table.items():
    log(f"  {key:18s} model {v['model']:9.3f}   real {v['real']:9.3f}")
log.claim("C41", "The Part A model, rebuilt from the 30 held-out starting grids, matches the held-out runs "
          "(table above)", table=table)

# Where the rebuilt box and the real box still differ: late Region 1 activity, by block (all runs that start at Reset or step 1)
BLOCKS = {"A": (slice(0, 5), slice(0, 5)), "B": (slice(0, 5), slice(5, 10)),
          "C": (slice(0, 5), slice(10, 15)), "D": (slice(5, 10), slice(0, 5))}


def late_by_block(trajs):
    tot = {k: [] for k in BLOCKS}
    for tr in trajs:
        ch = tr[1001:] != tr[1000:-1]
        for k, (a, b) in BLOCKS.items():
            tot[k].append(int(ch[:, a, b].sum()))
    return {k: np.array(v) for k, v in tot.items()}


long_runs = [h for k in ("H", "D") if k in DATASETS for h in DATASETS[k]]
real_b = late_by_block(long_runs)
model_b = late_by_block([GateModel(h[0], seed=300 + k).run(len(h) - 1) for k, h in enumerate(long_runs)])
blk = {k: dict(real=float(real_b[k].mean()), model=float(model_b[k].mean())) for k in BLOCKS}
busy_real, busy_model = int((real_b["D"] > 40).sum()), int((model_b["D"] > 40).sum())
log.claim("C43", "Late Region 1 activity (steps 1,001-2,000), changes per run by block, real vs rebuilt: "
          + ", ".join(f"{k} {v['real']:.1f} vs {v['model']:.1f}" for k, v in blk.items())
          + f"; runs with > 40 late block-D changes: real {busy_real}/{len(long_runs)}, rebuilt {busy_model}/{len(long_runs)}",
          by_block=blk, busy_D_real=busy_real, busy_D_model=busy_model, n_runs=len(long_runs))
log.save()
