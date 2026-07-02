"""Evaluable surface models extracted from ASM surface records.

Each class exposes the same tiny interface used by the tessellator:

* ``eval(u, v) -> (x, y, z)``  — point on the surface
* ``project(p) -> (u, v)``    — parameters of (a point near) ``p``
* ``periodic_u``              — whether u is an angle (needs unwrapping)

Supported: ``cone`` (cylinders, tapered cones, countersinks — analytic) and
``cyl_spl_sur`` (a B-spline directrix extruded along an axis).  Shared
surfaces (``ref N``) are resolved through :class:`RefPool`: ASM numbers every
subtype block (``{``) in stream order and ``ref N`` refers to the N-th one.
"""

from __future__ import annotations

import math

try:
    from . import sab, nurbs
except ImportError:
    import sab
    import nurbs


# --- vector helpers ---------------------------------------------------------
def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add3(a, b, c):
    return (a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2])


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _length(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    m = _length(a)
    return (a[0] / m, a[1] / m, a[2] / m) if m > 1e-30 else (0.0, 0.0, 0.0)


class Cone:
    """ACIS cone/cylinder: circle of radius ``r0`` swept along ``axis``.

    ``u`` is the angle around the axis, ``v`` the (signed) distance along it;
    the radius varies linearly: ``r(v) = r0 + v * taper``.  The taper sign is
    resolved empirically in :meth:`calibrate` because the stored sine/cosine
    sign convention is not obvious — the boundary points decide.
    """

    periodic_u = True

    def __init__(self, origin, axis, major, sine, cosine):
        self.origin = origin
        self.axis = _unit(axis)
        self.r0 = _length(major)
        self.M = _unit(major)
        self.N = _unit(_cross(self.axis, self.M))
        c = cosine if abs(cosine) > 1e-12 else 1.0
        self.taper = abs(sine / c)          # magnitude; sign set by calibrate()
        self._calibrated = self.taper == 0.0

    def eval(self, u, v):
        r = self.r0 + v * self.taper
        cu, su = math.cos(u), math.sin(u)
        return _add3(self.origin,
                     _scale(self.axis, v),
                     (r * (cu * self.M[0] + su * self.N[0]),
                      r * (cu * self.M[1] + su * self.N[1]),
                      r * (cu * self.M[2] + su * self.N[2])))

    def project(self, p):
        d = _sub(p, self.origin)
        v = _dot(d, self.axis)
        rad = _sub(d, _scale(self.axis, v))
        u = math.atan2(_dot(rad, self.N), _dot(rad, self.M))
        return (u, v)

    def calibrate(self, pts):
        """Fix the taper sign so ``r(v)`` fits the boundary points."""
        if self._calibrated or not pts:
            return
        err_pos = err_neg = 0.0
        for p in pts:
            d = _sub(p, self.origin)
            v = _dot(d, self.axis)
            r = _length(_sub(d, _scale(self.axis, v)))
            err_pos += abs(self.r0 + v * self.taper - r)
            err_neg += abs(self.r0 - v * self.taper - r)
        if err_neg < err_pos:
            self.taper = -self.taper
        self._calibrated = True

    def residual(self, p):
        u, v = self.project(p)
        return _length(_sub(self.eval(u, v), p))


class Torus:
    """ACIS torus: tube of radius ``r`` around a circle of radius ``R``.

    ``u`` is the angle around the main axis, ``v`` the angle around the tube
    (0 at the outer equator, +pi/2 towards the axis tip).  Handles lemon/apple
    tori (R < r) too -- the same equations apply.
    """

    periodic_u = True

    def __init__(self, center, axis, R, r, refdir):
        self.center = center
        self.axis = _unit(axis)
        self.R, self.r = R, r
        self.M = _unit(refdir)
        self.N = _unit(_cross(self.axis, self.M))

    def eval(self, u, v):
        rad = self.R + self.r * math.cos(v)
        cu, su = math.cos(u), math.sin(u)
        return _add3(self.center,
                     _scale(self.axis, self.r * math.sin(v)),
                     (rad * (cu * self.M[0] + su * self.N[0]),
                      rad * (cu * self.M[1] + su * self.N[1]),
                      rad * (cu * self.M[2] + su * self.N[2])))

    def project(self, p):
        d = _sub(p, self.center)
        h = _dot(d, self.axis)
        rad = _sub(d, _scale(self.axis, h))
        u = math.atan2(_dot(rad, self.N), _dot(rad, self.M))
        v = math.atan2(h, _length(rad) - self.R)
        return (u, v)

    def residual(self, p):
        u, v = self.project(p)
        return _length(_sub(self.eval(u, v), p))


class ExtrudedCurve:
    """``cyl_spl_sur``: B-spline directrix ``C(u)`` extruded along ``axis``.

    ``eval(u, v) = C(u) + v * axis``.  Projection finds the nearest directrix
    sample (dense precomputed table + local parabolic refinement).
    """

    periodic_u = False

    def __init__(self, deg, U, P, axis, samples=400):
        self.deg, self.U, self.P = deg, U, P
        self.axis = _unit(axis)
        t0, t1 = U[deg], U[len(P)]
        self._ts = [t0 + (t1 - t0) * k / samples for k in range(samples + 1)]
        self._pts = [nurbs.deboor(deg, U, P, t) for t in self._ts]

    def eval(self, u, v):
        c = nurbs.deboor(self.deg, self.U, self.P, u)
        return (c[0] + v * self.axis[0],
                c[1] + v * self.axis[1],
                c[2] + v * self.axis[2])

    def project(self, p):
        ax = self.axis
        best_k, best_d = 0, float("inf")
        for k, q in enumerate(self._pts):
            d = _sub(p, q)
            d = _sub(d, _scale(ax, _dot(d, ax)))   # distance ⟂ to axis
            dd = _dot(d, d)
            if dd < best_d:
                best_d, best_k = dd, k
        # local refinement between neighbours of the best sample
        lo = max(0, best_k - 1)
        hi = min(len(self._ts) - 1, best_k + 1)
        tlo, thi = self._ts[lo], self._ts[hi]
        for _ in range(24):
            tm1 = tlo + (thi - tlo) / 3
            tm2 = thi - (thi - tlo) / 3
            if self._perp_dist2(p, tm1) < self._perp_dist2(p, tm2):
                thi = tm2
            else:
                tlo = tm1
        u = (tlo + thi) / 2
        c = nurbs.deboor(self.deg, self.U, self.P, u)
        v = _dot(_sub(p, c), ax)
        return (u, v)

    def _perp_dist2(self, p, t):
        q = nurbs.deboor(self.deg, self.U, self.P, t)
        d = _sub(p, q)
        d = _sub(d, _scale(self.axis, _dot(d, self.axis)))
        return _dot(d, d)

    def residual(self, p):
        u, v = self.project(p)
        return _length(_sub(self.eval(u, v), p))


def _countable(name):
    """Whether a type name registers in the shared-geometry pool.

    ASM shares subtype data objects — interpolated curves (``*_int_cur``),
    parameter-space curves (``*_par_cur``) and spline surfaces (``*_spl_sur``).
    Each such definition, in stream order (attributes excluded), gets the next
    index; later ``ref N`` values refer back to definition N (0-based).
    Validated on the sample: the seven duplicate-body face surfaces resolve to
    the ``cyl_spl_sur`` definitions embedded in the first body's pcurves.
    """
    return name.endswith("_cur") or name.endswith("_sur")


class RefPool:
    """Resolves ``ref N`` shared-geometry references (see :func:`_countable`)."""

    def __init__(self, sabfile):
        self.defs = []                    # index -> (record, value_pos, name)
        for rec in sabfile.records:
            if rec.name == "ATTRIB_CUSTOM":
                continue
            for pos, v in enumerate(rec.values):
                if isinstance(v, sab.TypeName) and _countable(v.name):
                    self.defs.append((rec, pos, v.name))

    def get(self, n):
        if isinstance(n, int) and 0 <= n < len(self.defs):
            return self.defs[n]
        return None

    # -- resolution helpers -------------------------------------------------
    def curve_at(self, rec, pos, depth=0):
        """Parse the B-spline data of the curve definition at (rec, pos).

        Follows a ``ref N`` if the definition body holds a reference instead
        of inline data.  Returns ``(deg, U, P)`` or ``None``.
        """
        if depth > 4:
            return None
        end = _next_countable(rec, pos)
        for j in range(pos + 1, end):
            v = rec.values[j]
            if isinstance(v, sab.TypeName) and v.name in ("nubs", "nurbs"):
                parsed = nurbs.parse_bs_curve(rec.values, j, v.name == "nurbs")
                if parsed is not None:
                    return parsed
            if isinstance(v, sab.TypeName) and v.name == "ref":
                tgt = self.get(rec.values[j + 1] if j + 1 < len(rec.values) else None)
                if tgt is not None:
                    return self.curve_at(tgt[0], tgt[1], depth + 1)
        return None


def _next_countable(rec, pos):
    """Index one past this definition: the next countable name (or end)."""
    for j in range(pos + 1, len(rec.values)):
        v = rec.values[j]
        if isinstance(v, sab.TypeName) and _countable(v.name):
            return j
    return len(rec.values)


def _typenames(rec):
    return [v.name for v in rec.values if isinstance(v, sab.TypeName)]


def _extruded_from(rec, from_pos, pool):
    """Build an :class:`ExtrudedCurve` from a ``cyl_spl_sur`` definition.

    The directrix is either an inline nubs/nurbs or a ``ref`` into the pool;
    the extrusion axis is the first unit-length XYZ after the definition name.
    """
    parsed = None
    if pool is not None:
        parsed = pool.curve_at(rec, from_pos)
    if parsed is None:
        # fall back to any bs curve after the name (inline case, no pool)
        for j in range(from_pos + 1, len(rec.values)):
            v = rec.values[j]
            if isinstance(v, sab.TypeName) and v.name in ("nubs", "nurbs"):
                parsed = nurbs.parse_bs_curve(rec.values, j, v.name == "nurbs")
                if parsed is not None:
                    break
    if parsed is None:
        return None
    axis = None
    for v in rec.values[from_pos:]:
        if isinstance(v, tuple) and abs(_length(v) - 1.0) < 1e-9:
            axis = v
            break
    if axis is None:
        return None
    deg, U, P = parsed[:3]
    return ExtrudedCurve(deg, U, P, axis)


def _cone_from(rec):
    """Build a :class:`Cone` from a cone record's positional layout."""
    pos = rec.positions()
    if len(pos) < 3:
        return None
    origin, axis, major = pos[0], pos[1], pos[2]
    # sine/cosine: the two floats right after the ratio value that follows the
    # third XYZ (ratio is typically 1.0 for circular cones)
    floats = []
    seen = 0
    for v in rec.values:
        if isinstance(v, tuple):
            seen += 1
        elif seen >= 3 and isinstance(v, (int, float)) and not isinstance(v, bool):
            floats.append(float(v))
            if len(floats) >= 3:
                break
    if len(floats) < 3:
        return None
    _ratio, sine, cosine = floats[0], floats[1], floats[2]
    return Cone(origin, axis, major, sine, cosine)


def _torus_from(rec):
    """Build a :class:`Torus`: positions = center, axis, refdir; then R, r."""
    pos = rec.positions()
    if len(pos) < 3:
        return None
    center, axis, refdir = pos[0], pos[1], pos[2]
    # R and r are the two numbers between the axis and refdir tuples
    nums = []
    seen = 0
    for v in rec.values:
        if isinstance(v, tuple):
            seen += 1
            if seen >= 3:
                break
        elif seen == 2 and isinstance(v, (int, float)) and not isinstance(v, bool):
            nums.append(float(v))
    if len(nums) < 2:
        return None
    return Torus(center, axis, nums[0], nums[1], refdir)


def from_record(surf_rec, pool):
    """Build an evaluable surface from a face's surface record.

    Handles inline definitions and shared ``ref N`` ones.  Returns ``None``
    for unsupported kinds (planes are handled separately; blend surfaces have
    no closed form here and fall back to the tessellator's loft heuristic).
    """
    if surf_rec is None:
        return None
    if surf_rec.name == "cone":
        return _cone_from(surf_rec)
    if surf_rec.name == "torus":
        return _torus_from(surf_rec)
    names = _typenames(surf_rec)
    # dispatch on the FIRST *_spl_sur name: that is the face's own surface;
    # later ones are nested support surfaces (e.g. a blend's base cylinder)
    first_sur = next((n for n in names if n.endswith("_spl_sur")), None)
    if first_sur == "cyl_spl_sur":
        for j, v in enumerate(surf_rec.values):
            if isinstance(v, sab.TypeName) and v.name == "cyl_spl_sur":
                return _extruded_from(surf_rec, j, pool)
    elif first_sur is not None:
        return None                      # blends etc. -> loft heuristic
    if "ref" in names and pool is not None:
        for j, v in enumerate(surf_rec.values):
            if isinstance(v, sab.TypeName) and v.name == "ref":
                tgt = pool.get(surf_rec.values[j + 1]
                               if j + 1 < len(surf_rec.values) else None)
                if tgt is not None:
                    rec, pos, name = tgt
                    if name == "cyl_spl_sur":
                        return _extruded_from(rec, pos, pool)
                break
    return None
