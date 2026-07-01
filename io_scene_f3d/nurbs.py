"""B-spline (NURBS) parsing and evaluation for ASM ``nubs`` blocks.

ASM stores non-rational B-spline curves (``nubs``) and surfaces as a degree, a
knot vector given as (value, multiplicity) pairs, and a control net.  Two
conventions differ from the textbook form and are handled here:

* The first and last knot multiplicities are stored as ``degree`` rather than
  ``degree + 1``; we add one at each end to rebuild a clamped knot vector
  (validated: de Boor then reproduces edge endpoints exactly).

The public helpers evaluate a curve/surface with de Boor's algorithm and sample
them into polylines / point grids.
"""

from __future__ import annotations

import math

try:
    from . import sab
except ImportError:
    import sab


def _is_float(v):
    return isinstance(v, float)


def _build_knots(pairs):
    """(value, mult) pairs -> full clamped knot vector (+1 mult at each end)."""
    U = []
    last = len(pairs) - 1
    for idx, (val, mult) in enumerate(pairs):
        m = int(mult) + (1 if idx == 0 or idx == last else 0)
        U.extend([float(val)] * m)
    return U


def parse_nubs_curve(vals, start):
    """Parse a ``nubs`` curve starting at TypeName index ``start``.

    Returns ``(degree, knot_vector, control_points)`` or ``None``.
    """
    i = start + 1
    try:
        deg = int(vals[i]); i += 1
        i += 1                       # rational/other flag
        nk = int(vals[i]); i += 1
        pairs = []
        for _ in range(nk):
            pairs.append((vals[i], vals[i + 1])); i += 2
        cps = []
        while i + 2 < len(vals) and all(_is_float(vals[i + k]) for k in range(3)):
            cps.append((vals[i], vals[i + 1], vals[i + 2])); i += 3
    except (IndexError, TypeError, ValueError):
        return None
    if deg < 1 or len(cps) < 2:
        return None
    U = _build_knots(pairs)
    if len(U) != len(cps) + deg + 1:
        return None                  # knot/control mismatch -> not evaluable
    return deg, U, cps


def _find_span(deg, U, P, t):
    n = len(P) - 1
    if t <= U[deg]:
        return deg
    if t >= U[n + 1]:
        return n
    k = deg
    while k < n and U[k + 1] <= t:
        k += 1
    return k


def deboor(deg, U, P, t):
    """Evaluate a B-spline curve at parameter ``t`` (de Boor)."""
    k = _find_span(deg, U, P, t)
    d = [list(P[k - deg + j]) for j in range(deg + 1)]
    for r in range(1, deg + 1):
        for j in range(deg, r - 1, -1):
            i = k - deg + j
            den = U[i + deg - r + 1] - U[i]
            a = 0.0 if den == 0 else (t - U[i]) / den
            dj, dj1 = d[j], d[j - 1]
            for c in range(3):
                dj[c] = (1 - a) * dj1[c] + a * dj[c]
    return tuple(d[deg])


def sample_curve(deg, U, P, t0, t1, nseg):
    """Sample a B-spline curve from ``t0`` to ``t1`` into ``nseg+1`` points."""
    nseg = max(1, int(nseg))
    return [deboor(deg, U, P, t0 + (t1 - t0) * k / nseg) for k in range(nseg + 1)]


def curve_from_record(curve_rec):
    """Find and parse the ``nubs`` block inside a curve record; ``None`` if none."""
    for j, v in enumerate(curve_rec.values):
        if isinstance(v, sab.TypeName) and v.name == "nubs":
            return parse_nubs_curve(curve_rec.values, j)
    return None
