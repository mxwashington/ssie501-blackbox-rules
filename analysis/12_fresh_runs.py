"""12 - Fresh runs recorded by Group 5 on 16 November 2026 (collect/recorder.js), plus checks across all datasets.

Claims proven here: C45-C48.
  C45  every Part A rule on the fresh single-step runs, and how the page labels the Reset grid
  C46  "Next n Step" with n = 10 behaves like ten single steps
  C47  the block-D difference of C43, retested on the fresh runs
  C48  the values 10-13 in Region 3, against the rate the Region 3 table predicts
Run:  python analysis/12_fresh_runs.py   (needs data/group5_runs/; see data/README.md)
"""
import sys
import json
import glob
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from blackbox.data import load_group5, load_heldout, load_discovery, DATA
from blackbox.regions import AREA
from blackbox.report import Log
from blackbox.gate_model import GateModel, R3_TABLE
from blackbox.rule_checks import evaluate

log = Log("12_fresh_runs")
single = load_group5(cycles=1)
tens = load_group5(cycles=10)
files = sorted(glob.glob(str(DATA / "group5_runs" / "blackbox_*.json")))
meta = [json.loads(Path(f).read_text()) for f in files]
meta = [m for m in meta if m["label"] != "pilot"]
if not single:
    raise SystemExit("No Group 5 runs in data/group5_runs/ (see data/README.md)")
log(f"Fresh runs: {len(single)} single-step runs of {len(single[0]) - 1} steps, {len(tens)} runs of 10-step presses; "
    f"recorded {min(m['started'] for m in meta)} to {max(m['ended'] for m in meta)} (UTC)")

# C45: every rule on the fresh single-step runs
r = evaluate(single)
reset_steps = sorted({m["reset_step"] for m in meta})
log.claim("C45", f"Fresh runs: Region 1 gates right for {r['gate_ok']}/{r['n_ev']} changes, first step after Reset "
          f"{r['first']['r1_gate_ok']}/{r['first']['r1_changes']}; settled cells jumped {r['settled']['agree'][1]} times in "
          f"{r['settled']['agree'][0]} cell-steps while the gate agreed; Region 2 minimum part {r['pt']['partner_found']}/"
          f"{r['pt']['changes']} ({r['pt']['rises']} rises), multiply part {r['fo']['product_found']}/{r['fo']['changes']}; "
          f"Region 4 {r['sw']['strict_in_order']}/{r['sw']['swaps']} swaps. The page labels the Reset grid "
          f"'Current step: {reset_steps}'",
          gate_ok=r["gate_ok"], n=r["n_ev"], first=r["first"], settled=r["settled"], pt=r["pt"], fo=r["fo"],
          sw=r["sw"], reset_step_labels=reset_steps, n_runs=len(single))


def H(v):
    _, c = np.unique(v, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())


def entropies(traj, every):
    """Entropy of each region at total steps 500, 1000, 2000 (frames are `every` steps apart)."""
    return {s: [H(traj[s // every][AREA == a]) for a in (1, 2, 3, 4)] for s in (500, 1000, 2000)}


# C46: 10-step presses vs single steps
if tens:
    e1 = [entropies(t, 1) for t in single]
    e10 = [entropies(t, 10) for t in tens]
    comp = {s: dict(single=np.mean([e[s] for e in e1], axis=0).round(3).tolist(),
                    ten=np.mean([e[s] for e in e10], axis=0).round(3).tolist()) for s in (500, 1000, 2000)}
    # within one press of 10 steps, Regions 1-3 can change at most 10 cells each
    maxcells = {a: max(int((t[k + 1] != t[k])[AREA == a].sum()) for t in tens for k in range(len(t) - 1)) for a in (1, 2, 3)}
    held = all(np.array_equal(np.sort(t[k][AREA == 4]), np.sort(t[0][AREA == 4])) for t in tens for k in range(len(t)))
    gap = max(abs(a - b) for s in comp for a, b in zip(comp[s]["single"], comp[s]["ten"]))
    log.claim("C46", f"10-step presses ({len(tens)} runs) vs single steps ({len(single)} runs): region entropies at "
              f"500/1000/2000 steps differ by at most {gap:.2f} bits; at most {max(maxcells.values())} cells of one "
              f"region changed in one press (10 steps allow 10); Region 4's value counts never changed: {held}",
              by_step=comp, max_cells_per_press=maxcells, region4_counts_constant=held, max_gap_bits=gap)

# C47: block D, retested on the fresh runs
B = {"A": (slice(0, 5), slice(0, 5)), "B": (slice(0, 5), slice(5, 10)),
     "C": (slice(0, 5), slice(10, 15)), "D": (slice(5, 10), slice(0, 5))}


def late(trajs):
    return {k: np.array([int((t[1001:] != t[1000:-1])[:, a, b].sum()) for t in trajs]) for k, (a, b) in B.items()}


rf = late(single)
mf = late([GateModel(t[0], seed=4000 + 10 * k + rep).run(len(t) - 1) for k, t in enumerate(single) for rep in range(5)])
blk = {k: dict(real=float(rf[k].mean()), model=float(mf[k].mean())) for k in B}
log.claim("C47", "Fresh runs, late Region 1 changes per run by block (real vs rebuilt from the same Reset grids, 5 runs each): "
          + ", ".join(f"{k} {v['real']:.1f} vs {v['model']:.1f}" for k, v in blk.items())
          + f"; runs with > 40 late block-D changes: real {(rf['D'] > 40).sum()}/{len(single)}, "
            f"rebuilt {(mf['D'] > 40).sum()}/{len(mf['D'])}",
          by_block=blk, busy_D_real=int((rf["D"] > 40).sum()), busy_D_model=int((mf["D"] > 40).sum()),
          n_real=len(single), n_model=len(mf["D"]))

# C48: values 10 and above in Region 3, all datasets, after step 500 (Region 3 has been redrawn many times by then)
pools = {"held-out": load_heldout(), "fresh": single}
try:
    pools["Group 1"] = load_discovery()
except FileNotFoundError:
    pass
p_hi = float(R3_TABLE[10:].sum())
counts = {}
for name, runs_ in pools.items():
    cells = sum(int(t[500:][:, AREA == 3].size) for t in runs_)
    hi = sum(int((t[500:][:, AREA == 3] >= 10).sum()) for t in runs_)
    other = sum(int((t[:, AREA != 3] >= 10).sum()) for t in runs_)
    counts[name] = dict(cell_steps=cells, values_10_plus=hi, share=hi / cells, outside_region3=other)
tot_c = sum(v["cell_steps"] for v in counts.values())
tot_h = sum(v["values_10_plus"] for v in counts.values())
log.claim("C48", f"Values 10-13: {tot_h} of {tot_c} Region 3 cell-steps after step 500 ({tot_h / tot_c:.5f}); the Region 3 "
          f"table predicts {p_hi:.5f}; outside Region 3: {sum(v['outside_region3'] for v in counts.values())}",
          by_dataset=counts, predicted=p_hi, observed=tot_h / tot_c)
log.save()
