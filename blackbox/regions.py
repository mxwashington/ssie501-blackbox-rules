"""The four areas of the 20x20 grid (row, col are 0-based; row 0 is the top).

Each area is an L-shape of exactly 100 cells: the boundary moves 5 columns at the
half-way row. Numbering follows our Part 1 presentation:

  Area 1  top-left    rows 0-4 cols 0-14  + rows 5-9  cols 0-4    (freezes to 4s)
  Area 2  top-right   rows 0-4 cols 15-19 + rows 5-9  cols 5-19   (goes to 0s)
  Area 3  bottom-left rows 10-14 cols 0-14 + rows 15-19 cols 0-4  (never freezes)
  Area 4  bottom-right rows 10-14 cols 15-19 + rows 15-19 cols 5-19 (swaps)

(Group 1's notebook uses the same shapes but calls our Area 3 "R4" and our Area 4 "R3".)
"""
import numpy as np

AREA = np.zeros((20, 20), dtype=int)
AREA[0:5, 0:15] = 1
AREA[5:10, 0:5] = 1
AREA[0:5, 15:20] = 2
AREA[5:10, 5:20] = 2
AREA[10:15, 0:15] = 3
AREA[15:20, 0:5] = 3
AREA[10:15, 15:20] = 4
AREA[15:20, 5:20] = 4

MASK = {a: AREA == a for a in (1, 2, 3, 4)}

VN = [(-1, 0), (1, 0), (0, -1), (0, 1)]                      # 4 side neighbours
MOORE = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]  # 8 neighbours


def neighbour_values(g, r, c, nb=MOORE, wrap=False):
    """Values of the neighbours of (r, c). wrap=False drops cells off the grid edge."""
    out = []
    for dr, dc in nb:
        rr, cc = r + dr, c + dc
        if wrap:
            out.append(int(g[rr % 20, cc % 20]))
        elif 0 <= rr < 20 and 0 <= cc < 20:
            out.append(int(g[rr, cc]))
    return out
