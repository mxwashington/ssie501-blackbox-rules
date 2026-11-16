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


def load_group5(cycles=1):
    """Group 5's own recordings (data/group5_runs/*.json, made with collect/recorder.js).

    Each file holds the Reset grid and, for every press of "Next n Step", the cells that changed.
    Returns the runs recorded with the given n as arrays of shape (presses + 1, 20, 20):
    frame 0 is the Reset grid, frame k the grid after k presses (k x n steps).
    """
    import json
    out = []
    for f in sorted(glob.glob(str(DATA / "group5_runs" / "blackbox_*.json"))):
        run = json.loads(Path(f).read_text())
        if run.get("cycles") != cycles or run.get("label") == "pilot":
            continue
        assert not run["errors"], f"{f}: recorder reported errors"
        assert run["steps"] == [run["reset_step"] + cycles * (k + 1) for k in range(len(run["steps"]))], f"{f}: step counter gaps"
        g = np.array(run["reset"], dtype=np.int64)
        frames = [g.reshape(20, 20).copy()]
        for ch in run["changes"]:
            for i, v in ch:
                g[i] = v
            frames.append(g.reshape(20, 20).copy())
        assert (frames[-1].ravel() == np.array(run["final"])).all(), f"{f}: replay does not reach the final grid"
        out.append(np.array(frames))
    return out
