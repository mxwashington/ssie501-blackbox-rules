"""Load the three recorded datasets as integer arrays of shape (frames, 20, 20).

Every dataset records the grid after every single step (n = 1), so frame t+1 is
exactly one step after frame t.

  DISCOVERY  - Group 1 scraper runs (61 runs x 2,001 frames). Used to FIND rules.
  HELDOUT    - Shared dataset compiled by W. Zhang (30 runs x 2,000 frames).
               Different starting grids; used only to TEST rules found on DISCOVERY.
  GROUP2     - Group 2 scraper frames (1 run x 41 frames). Small extra check.
"""
from pathlib import Path
import glob
import csv
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def load_discovery():
    files = sorted(glob.glob(str(DATA / "group1_runs" / "run_*.npy")))
    if not files:
        raise FileNotFoundError("Put Group 1 run_XX.npy files in data/group1_runs/ (see data/README.md)")
    return [np.load(f).astype(np.int64) for f in files]


def load_heldout():
    files = sorted(glob.glob(str(DATA / "shared_dataset" / "trajectory_*_FULLBOX.npz")))
    if not files:
        raise FileNotFoundError("Put trajectory_XXX_FULLBOX.npz files in data/shared_dataset/ (see data/README.md)")
    out = []
    for f in files:
        z = np.load(f)
        steps = z["recorded_steps"]
        assert set(np.diff(steps).tolist()) == {1}, f"{f}: frames are not consecutive single steps"
        out.append(z["full_history"].astype(np.int64))
    return out


def load_group2():
    path = DATA / "group2" / "frames.csv"
    if not path.exists():
        return []
    rows = list(csv.reader(open(path)))[1:]
    rows = [r for r in rows if r]
    steps = [int(r[0]) for r in rows]
    assert set(np.diff(steps).tolist()) == {1}
    return [np.array([[int(x) for x in r[1:]] for r in rows]).reshape(-1, 20, 20)]
