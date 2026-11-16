"""Turn recorded runs into per-step change events."""
import numpy as np
from .regions import AREA, MASK, MOORE, neighbour_values


def changes_per_step(runs, area):
    """List with the number of cells of `area` that changed at each step, over all runs."""
    m = MASK[area]
    out = []
    for a in runs:
        out.extend(((a[1:] != a[:-1])[:, m]).sum(1).tolist())
    return out


def single_events(runs, area, nb=MOORE, wrap=False):
    """Yield one dict per step in which exactly one cell of `area` changed.

    keys: old, new, nbrs (neighbour values at the step's start), r, c, run, t
    """
    m = MASK[area]
    for ri, a in enumerate(runs):
        D = (a[1:] != a[:-1]) & m
        counts = D.reshape(len(D), -1).sum(1)
        for t in np.flatnonzero(counts == 1):
            r, c = np.argwhere(D[t])[0]
            g = a[t]
            yield dict(old=int(g[r, c]), new=int(a[t + 1][r, c]),
                       nbrs=neighbour_values(g, r, c, nb, wrap),
                       r=int(r), c=int(c), run=ri, t=int(t))


def swap_events(runs):
    """Yield (run, t, (r1,c1), (r2,c2)) for every step in which exactly two Area-4 cells changed."""
    m = MASK[4]
    for ri, a in enumerate(runs):
        D = (a[1:] != a[:-1]) & m
        counts = D.reshape(len(D), -1).sum(1)
        for t in np.flatnonzero(counts == 2):
            (r1, c1), (r2, c2) = np.argwhere(D[t])
            yield ri, int(t), (int(r1), int(c1)), (int(r2), int(c2))
