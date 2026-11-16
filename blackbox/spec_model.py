"""The Black Box model exactly as written in the model specification document (Part A).

This file is a deliberately literal, line-by-line translation of the written specification, kept
separate from blackbox/simulator.py. It reads ONLY the printed tables (results/model_tables.json),
so it tests whether the WRITTEN model - with its rounded numbers - reproduces the box.
"""
import numpy as np
from .regions import AREA, MASK

MOORE = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
CELLS = {a: [tuple(x) for x in np.argwhere(MASK[a])] for a in (1, 2, 3, 4)}


def wait_bin(wait):
    """Section 'Region 1', step 3: which column of the jump table to use."""
    if wait is None:
        return "first"
    for lo, hi, name in [(0, 49, "0-49"), (50, 99, "50-99"), (100, 149, "100-149"), (150, 199, "150-199")]:
        if lo <= wait <= hi:
            return name
    return "200+"


def draw(table, u):
    """Pick a value from a probability table {value: probability} with one uniform number u in [0, 1)."""
    acc = 0.0
    for value, prob in table:
        acc += prob
        if u < acc:
            return value
    return table[-1][0]


def matching_neighbours(g, r, c, value):
    """Count of the 8 surrounding cells (inside the 20x20 grid, any region) holding `value`."""
    n = 0
    for dr, dc in MOORE:
        rr, cc = r + dr, c + dc
        if 0 <= rr < 20 and 0 <= cc < 20 and g[rr][cc] == value:
            n += 1
    return n


def step(g, last_move, t, T, rng):
    """One step of the box. g = grid at the start of the step (list of lists); returns the new grid.
    last_move[(r, c)] = step number of the cell's last change (Region 1 only; missing = never moved)."""
    new = [row[:] for row in g]
    # ---- Region 1
    r, c = CELLS[1][rng.integers(100)]
    v = g[r][c]
    if v in (4, 5):
        if rng.random() < T["region1"]["exit"][str(v)]:
            new[r][c] = 9 if v == 4 else 0
            last_move[(r, c)] = t
    else:
        branch = "low" if v <= 3 else "high"
        wait = None if (r, c) not in last_move else t - last_move[(r, c)]
        p_jump = T["region1"]["jump"][branch][wait_bin(wait)]
        if rng.random() < p_jump:
            new[r][c] = 9 if branch == "low" else 0
        else:
            new[r][c] = v + 1 if branch == "low" else v - 1
        last_move[(r, c)] = t
    # ---- Region 2
    r, c = CELLS[2][rng.integers(100)]
    v = g[r][c]
    if v != 0 and rng.random() < T["region2"]["p_change"][str(v)]:
        row = [(int(k), p) for k, p in T["region2"]["target"][str(v)].items()]
        new[r][c] = draw(row, rng.random())
    # ---- Region 3
    r, c = CELLS[3][rng.integers(100)]
    new[r][c] = draw([(int(k), p) for k, p in T["region3"]["p"].items()], rng.random())
    # ---- Region 4
    i, j = rng.choice(100, 2, replace=False)
    (r1, c1), (r2, c2) = CELLS[4][i], CELLS[4][j]
    a, b = g[r1][c1], g[r2][c2]
    if a != b:
        before = matching_neighbours(g, r1, c1, a) + matching_neighbours(g, r2, c2, b)
        trial = [row[:] for row in g]
        trial[r1][c1], trial[r2][c2] = b, a
        after = matching_neighbours(trial, r1, c1, b) + matching_neighbours(trial, r2, c2, a)
        gain = after - before
        if gain >= 1 and rng.random() < T["region4"]["accept"][str(min(gain, 5))]:
            new[r1][c1], new[r2][c2] = b, a
    return new


def run(start, steps, T, rng):
    g = [list(map(int, row)) for row in start]
    last_move = {}
    frames = [np.array(g)]
    for t in range(steps):
        g = step(g, last_move, t, T, rng)
        frames.append(np.array(g))
    return np.array(frames)
