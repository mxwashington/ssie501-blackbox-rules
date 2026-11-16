"""The Part A model as code: Region 1 gates, Region 2 minimum/multiply parts, Region 3 redraw, Region 4 swaps.

This is a clean re-implementation of the rules in MODEL_SPECIFICATION.md, Part A. The Region 1 gates,
the Region 2 split and the strict Region 4 rule come from Marcus Washington's simulator
(SSIE501_Blackbox_Simulator_Explained_v2.ipynb); unused leftovers from his earlier versions are removed.
Rows and columns count from 0 here (row 0 is the top).

    from blackbox.gate_model import GateModel, r1_targets
    traj = GateModel(start_grid, seed=1).run(2000)     # (2001, 20, 20)
"""
import numpy as np

from blackbox.regions import AREA

N = 20
CELLS = {a: np.flatnonzero(AREA.ravel() == a) for a in (1, 2, 3, 4)}
PT = np.zeros((N, N), bool)
PT[0:10, 15:20] = True                      # Region 2, rows 1-10 x columns 16-20 (minimum part)
FO = (AREA == 2) & ~PT                      # Region 2, rows 6-10 x columns 6-15 (multiply part)
FO_CELLS = np.flatnonzero(FO.ravel())
PT_SET = set(np.flatnonzero(PT.ravel()).tolist())
BOTTOM = np.arange(10 * N, N * N)           # rows 11-20, any column: partners of the minimum part

# Region 3 redraw table (MODEL_SPECIFICATION A3), values 0..11
R3_TABLE = np.array([0.3634, 0.3292, 0.1690, 0.0756, 0.0338, 0.0159, 0.0076,
                     0.0032, 0.0015, 0.0005, 0.0002, 0.0001])
R3_TABLE = R3_TABLE / R3_TABLE.sum()

OFF8 = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)]
NB8 = [[(r + dr) * N + (c + dc) for dr, dc in OFF8 if 0 <= r + dr < N and 0 <= c + dc < N]
       for r in range(N) for c in range(N)]

# ---- Region 1 wiring -------------------------------------------------------------------------
# Blocks are the sixteen 5x5 squares, numbered 0-15 left to right, top to bottom.
# Region 1 = blocks 0, 1, 2 (rows 1-5) and 4 (rows 6-10, columns 1-5).
SOURCE = {0: 0, 1: 1, 2: 0, 4: 4}             # block whose cells feed the gates; block 2 reads block 0
CP_CQ = {0: (5, 6), 1: (8, 9), 2: (13, 14), 4: (3, 12)}   # blocks whose centre cells are inputs p and q
C0 = {0: 7, 1: 10, 2: 15, 4: 11}             # block whose centre cell is input c


def _band(v):
    return int(v >= 5)                        # LOW = 0 (values 0-4), HIGH = 1 (5 and up, including 10+)


def _centre(g, blk):
    return _band(g[(blk // 4) * 5 + 2, (blk % 4) * 5 + 2])


def r1_targets(g, b):
    """5x5 array of target halves (1 = HIGH, 0 = LOW) for Region 1 block b on grid g (20x20)."""
    sr, sc = (SOURCE[b] // 4) * 5, (SOURCE[b] % 4) * 5
    s = (g[sr:sr + 5, sc:sc + 5] >= 5).astype(int)

    def S(i, j):                              # source cell at block position (i, j), counted from 1
        return s[i - 1, j - 1]
    p, q, c = _centre(g, CP_CQ[b][0]), _centre(g, CP_CQ[b][1]), _centre(g, C0[b])
    n = lambda x: 1 - x
    T = [[n(S(3, 1)), S(1, 5) & n(S(1, 3)), S(3, 2), S(5, 2) & n(S(4, 5)), (S(1, 4) | S(2, 2) | S(2, 3)) & n(S(5, 3))],
         [S(2, 5), p & q, S(4, 4) & p, S(4, 1), S(1, 2) | S(3, 1)],
         [S(4, 2) & n(S(1, 3)), S(2, 4) & n(S(1, 4)), n(S(1, 1)), n(S(1, 5)) & (S(2, 2) | S(4, 4)), S(3, 4)],
         [S(3, 4) & S(5, 4), n(S(3, 2)) & (S(3, 1) | S(3, 5)), S(2, 5), c, S(3, 2) & n(S(4, 3))],
         [S(1, 4), S(1, 1) & S(2, 1) & n(S(4, 5)), S(3, 2), n(S(4, 3)) & (S(3, 2) | S(5, 4)), S(3, 2)]]
    return np.array(T, int)


def r1_target_grid(g):
    """20x20 array: target half of every Region 1 cell (-1 elsewhere)."""
    out = np.full((N, N), -1)
    for b in (0, 1, 2, 4):
        r0, c0 = (b // 4) * 5, (b % 4) * 5
        out[r0:r0 + 5, c0:c0 + 5] = r1_targets(g, b)
    return out


def r1_next(v, target_high):
    """New value of a picked Region 1 cell holding v (0-9) whose gate says target_high (0/1)."""
    high = v >= 5
    if target_high != high:
        return 0 if high else 9               # jump to the start of the other branch
    if v <= 3:
        return v + 1                          # advance up the low branch
    if v >= 6:
        return v - 1                          # advance down the high branch
    return v                                  # 4 or 5 in its target half: stays


def matches(g, i, v):
    return sum(1 for j in NB8[i] if g[j] == v)


class GateModel:
    def __init__(self, grid0, seed):
        self.rng = np.random.default_rng(seed)
        self.g = np.array(grid0, dtype=np.int64).ravel().copy()

    def step(self):
        g, rng = self.g, self.rng
        # Region 1
        i = CELLS[1][rng.integers(100)]
        if g[i] <= 9:
            r, c = divmod(i, N)
            b = (r // 5) * 4 + c // 5
            t = r1_targets(g.reshape(N, N), b)[r % 5, c % 5]
            g[i] = r1_next(int(g[i]), int(t))
        # Region 2
        i = CELLS[2][rng.integers(100)]
        if i in PT_SET:
            g[i] = min(g[i], g[BOTTOM[rng.integers(len(BOTTOM))]])
        else:
            g[i] = (g[i] * g[FO_CELLS[rng.integers(len(FO_CELLS))]]) % 10
        # Region 3
        g[CELLS[3][rng.integers(100)]] = rng.choice(len(R3_TABLE), p=R3_TABLE)
        # Region 4
        cells = CELLS[4]
        vals = g[cells]
        if np.all(vals == vals[0]):
            return
        while True:
            a, b = rng.integers(100, size=2)
            if a != b and vals[a] != vals[b]:
                break
        i, j = cells[a], cells[b]
        vi, vj = g[i], g[j]
        bi, bj = matches(g, i, vi), matches(g, j, vj)
        g[i], g[j] = vj, vi
        ai, aj = matches(g, i, vj), matches(g, j, vi)
        if not (ai >= bi and aj >= bj and ai + aj > bi + bj):
            g[i], g[j] = vi, vj

    def run(self, nsteps):
        out = np.empty((nsteps + 1, N, N), dtype=np.int64)
        out[0] = self.g.reshape(N, N)
        for t in range(1, nsteps + 1):
            self.step()
            out[t] = self.g.reshape(N, N)
        return out
