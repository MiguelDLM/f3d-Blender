"""B-spline parsing and evaluation for ASM ``nubs``/``nurbs`` blocks.

ASM stores B-spline curves as a degree, a knot vector given as
(value, multiplicity) pairs, and control points: plain (x, y, z) triples for
non-rational ``nubs``, homogeneous (x, y, z, w) quadruples for rational
``nurbs``.  Two conventions differ from the textbook form:

* The first and last knot multiplicities are stored as ``degree`` rather than
  ``degree + 1``; we add one at each end to rebuild a clamped knot vector
  (validated: de Boor then reproduces edge endpoints exactly).
* The number of control points is not stored; it follows from
  ``len(knots) - degree - 1``, which also guards against over-reading into
  trailing fields (e.g. the fit tolerance).

Everything is evaluated homogeneously (weight 1 for ``nubs``) with de Boor.
"""

from __future__ import annotations

try:
    from . import sab
except ImportError:
    import sab


def _build_knots(pairs):
    """(value, mult) pairs -> full clamped knot vector (+1 mult at each end)."""
    U = []
    last = len(pairs) - 1
    for idx, (val, mult) in enumerate(pairs):
        m = int(mult) + (1 if idx == 0 or idx == last else 0)
        U.extend([float(val)] * m)
    return U


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def parse_bs_curve(vals, start, rational):
    """Parse a ``nubs``/``nurbs`` curve at TypeName index ``start``.

    Returns ``(degree, knot_vector, control_points4, period)`` with homogeneous
    (x*w, y*w, z*w, w) control points, or ``None`` on any inconsistency.
    ``period`` is the knot-domain length for closed/periodic curves (flag=2)
    and ``None`` otherwise: edges on a closed curve may carry parameters
    outside (e.g. the negative mirror of) the stored domain, which must be
    wrapped modulo the period before evaluating.
    """
    i = start + 1
    try:
        deg = int(vals[i]); i += 1
        flag = vals[i]; i += 1       # 0 = open, 2 = closed/periodic
        nk = int(vals[i]); i += 1
        pairs = []
        for _ in range(nk):
            pairs.append((vals[i], vals[i + 1])); i += 2
        U = _build_knots(pairs)
        ncp = len(U) - deg - 1
        dim = 4 if rational else 3
        if deg < 1 or ncp < 2 or i + ncp * dim > len(vals):
            return None
        P = []
        for _ in range(ncp):
            comps = vals[i:i + dim]
            if not all(_is_num(c) for c in comps):
                return None
            i += dim
            if rational:
                x, y, z, w = (float(c) for c in comps)
                P.append((x * w, y * w, z * w, w))
            else:
                x, y, z = (float(c) for c in comps)
                P.append((x, y, z, 1.0))
    except (IndexError, TypeError, ValueError):
        return None
    period = (U[-1] - U[0]) if flag == 2 else None
    return deg, U, P, period


def _find_span(deg, U, n_ctrl, t):
    n = n_ctrl - 1
    if t <= U[deg]:
        return deg
    if t >= U[n + 1]:
        return n
    k = deg
    while k < n and U[k + 1] <= t:
        k += 1
    return k


def deboor(deg, U, P, t):
    """Evaluate a B-spline curve at ``t`` (homogeneous de Boor) -> (x, y, z)."""
    k = _find_span(deg, U, len(P), t)
    d = [list(P[k - deg + j]) for j in range(deg + 1)]
    for r in range(1, deg + 1):
        for j in range(deg, r - 1, -1):
            i = k - deg + j
            den = U[i + deg - r + 1] - U[i]
            a = 0.0 if den == 0 else (t - U[i]) / den
            dj, dj1 = d[j], d[j - 1]
            for c in range(4):
                dj[c] = (1 - a) * dj1[c] + a * dj[c]
    x, y, z, w = d[deg]
    if w != 0 and w != 1.0:
        return (x / w, y / w, z / w)
    return (x, y, z)


def sample_curve(deg, U, P, t0, t1, nseg, period=None):
    """Sample a B-spline curve from ``t0`` to ``t1`` into ``nseg+1`` points.

    For closed/periodic curves pass ``period``: each parameter is wrapped
    into the stored knot domain (edges may run e.g. over [-T, 0]).
    """
    nseg = max(1, int(nseg))
    ts = [t0 + (t1 - t0) * k / nseg for k in range(nseg + 1)]
    if period:
        lo = U[deg]
        ts = [lo + ((t - lo) % period) for t in ts]
    return [deboor(deg, U, P, t) for t in ts]


def curve_from_record(curve_rec):
    """Find and parse the first ``nubs``/``nurbs`` block in a curve record."""
    for j, v in enumerate(curve_rec.values):
        if isinstance(v, sab.TypeName) and v.name in ("nubs", "nurbs"):
            parsed = parse_bs_curve(curve_rec.values, j, v.name == "nurbs")
            if parsed is not None:
                return parsed
    return None
