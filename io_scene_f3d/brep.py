"""Reconstruct B-rep topology and pull geometry from parsed ASM records.

Consumes a :class:`sab.SabFile` (a flat list of records referencing each other
by index) and rebuilds the solid hierarchy::

    body -> lump -> shell -> face -> loop -> coedge -> edge -> vertex -> point

For every face it yields its boundary loops as ordered rings of 3D points
(outer loop first, then any inner/hole loops), plus the underlying surface kind
and — for analytic surfaces — its defining frame.  That is exactly what the
tessellator needs.

Field offsets below were verified against the sample by checking, for each
value position of every record of a kind, the *type* of record it points at
(see docs/FORMAT.md).  ``values[0]`` is always the entity :class:`sab.TypeName`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

try:  # package import (inside Blender add-on)
    from . import sab
except ImportError:  # standalone (tests)
    import sab


# --- small vector helpers -------------------------------------------------
def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _length(a):
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _unit(a):
    m = _length(a)
    return (a[0] / m, a[1] / m, a[2] / m) if m > 1e-30 else (0.0, 0.0, 0.0)

# --- verified value-position offsets --------------------------------------
BODY_LUMP = 4
LUMP_SHELL, LUMP_BODY = 5, 6
SHELL_FACE, SHELL_LUMP = 6, 8
FACE_NEXT, FACE_LOOP, FACE_SHELL, FACE_SURFACE = 4, 5, 6, 8
LOOP_NEXT, LOOP_COEDGE, LOOP_FACE = 4, 5, 6
COEDGE_NEXT, COEDGE_PREV, COEDGE_EDGE, COEDGE_SENSE, COEDGE_LOOP = 4, 5, 7, 8, 9
EDGE_VSTART, EDGE_VEND, EDGE_CURVE = 4, 6, 9
VERTEX_POINT = 6

# guards against malformed / cyclic data
_MAX_RING = 100000


def _ratio_value(rec) -> float:
    """Minor/major radius ratio of an ellipse: the number after its 3 XYZs."""
    seen_pos = 0
    for v in rec.values:
        if isinstance(v, tuple):
            seen_pos += 1
        elif seen_pos >= 3 and isinstance(v, (int, float)) and not isinstance(v, bool):
            return float(v)
    return 1.0


def _nubs_control_points(curve):
    """Extract the (x, y, z) control points of the curve's ``nubs`` block.

    In the flattened record the control points are the run of consecutive
    doubles that immediately precedes the trailing ``null_surface`` markers
    (the knot values are separated from them by integer multiplicities).
    """
    vals = curve.values
    start = None
    for i, v in enumerate(vals):
        if isinstance(v, sab.TypeName) and v.name == "nubs":
            start = i + 1
            break
    if start is None:
        return []
    # find the end: first TypeName after the nubs header
    end = len(vals)
    for i in range(start, len(vals)):
        if isinstance(vals[i], sab.TypeName):
            end = i
            break
    # the control-point block is the last contiguous run of floats in [start,end)
    run_end = end
    while run_end > start and not _is_num(vals[run_end - 1]):
        run_end -= 1
    run_start = run_end
    while run_start > start and _is_num(vals[run_start - 1]):
        run_start -= 1
    nums = [float(vals[i]) for i in range(run_start, run_end)]
    cps = [(nums[k], nums[k + 1], nums[k + 2]) for k in range(0, len(nums) - 2, 3)]
    return cps


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _sample_nubs(curve, arc_seg_fn):
    """Sample a spline (intcurve/nubs) edge.

    Disabled for now: a naive control-polygon approximation produces large
    overshoots (a B-spline's control points can lie far from the curve), so we
    fall back to the edge endpoints (a straight chord) which is safe.  Proper
    de Boor evaluation with a correctly-parsed knot vector is future work; the
    control points can be recovered via :func:`_nubs_control_points`.
    """
    return None


@dataclass
class Face:
    surface_kind: str                       # plane / cone / spline / ...
    loops: list = field(default_factory=list)   # list[list[(x,y,z)]]; [0]=outer
    normal: tuple | None = None             # analytic surface normal, if known


@dataclass
class Body:
    faces: list = field(default_factory=list)


class Brep:
    def __init__(self, sabfile: "sab.SabFile", deviation: float = 0.05,
                 min_arc_segments: int = 24):
        self.f = sabfile
        self.deviation = deviation        # max chord error in model units
        self.min_arc_segments = min_arc_segments  # per full 2*pi turn

    # -- low level helpers --------------------------------------------------
    def _rec(self, ref):
        if not isinstance(ref, sab.Ref) or ref.is_null:
            return None
        return self.f.resolve(ref)

    def _ref_at(self, rec, pos):
        if rec is None or pos >= len(rec.values):
            return None
        return self._rec(rec.values[pos])

    def _vertex_point(self, vertex):
        pt = self._ref_at(vertex, VERTEX_POINT)
        if pt is None:
            return None
        pos = pt.positions()
        return pos[0] if pos else None

    def _walk(self, first, next_pos):
        """Follow a ``next`` pointer chain, returning the ordered records."""
        out = []
        seen = set()
        cur = first
        while cur is not None and cur.index not in seen and len(out) < _MAX_RING:
            out.append(cur)
            seen.add(cur.index)
            cur = self._ref_at(cur, next_pos)
        return out

    # -- geometry -----------------------------------------------------------
    def _surface_info(self, surf):
        """Return (kind, normal|None) for a face's surface record."""
        if surf is None:
            return "unknown", None
        kind = surf.name          # first type name: plane / cone / spline
        normal = None
        if kind == "plane":
            pos = surf.positions()
            if len(pos) >= 2:
                normal = pos[1]           # origin, normal, u_ref
        return kind, normal

    def _arc_segments(self, radius, dtheta):
        """Number of chords for an arc of |dtheta| radians at given radius."""
        n_min = max(1, int(math.ceil(abs(dtheta) / (2 * math.pi)
                                     * self.min_arc_segments)))
        if radius > 1e-9 and self.deviation > 0:
            arg = 1.0 - self.deviation / radius
            if -1.0 < arg < 1.0:
                step = 2.0 * math.acos(arg)
                if step > 1e-6:
                    n_min = max(n_min, int(math.ceil(abs(dtheta) / step)))
        return min(max(n_min, 1), 512)

    def _sample_ellipse(self, ell, t0, t1):
        """Sample an ellipse/circle from param ``t0`` to ``t1`` (radians)."""
        pos = ell.positions()
        if len(pos) < 3:
            return []
        center, normal, major = pos[0], pos[1], pos[2]
        ratio = _ratio_value(ell)
        a = _length(major)
        if a < 1e-12:
            return []
        minor_dir = _unit(_cross(normal, _unit(major)))
        b = a * ratio
        n = self._arc_segments(max(a, b), t1 - t0)
        pts = []
        for k in range(n + 1):
            t = t0 + (t1 - t0) * k / n
            ct, st = math.cos(t), math.sin(t)
            pts.append(_add(center,
                            _add(_scale(major, ct),
                                 _scale(minor_dir, b * st))))
        return pts

    def _sample_edge(self, edge, reverse):
        """Return the polyline of ``edge`` traversed in the coedge direction.

        The first point is the coedge's start vertex, the last its end vertex;
        interior points come from evaluating the edge's curve.
        """
        v0 = self._vertex_point(self._ref_at(edge, EDGE_VSTART))
        v1 = self._vertex_point(self._ref_at(edge, EDGE_VEND))
        t0 = edge.values[5] if isinstance(edge.values[5], float) else 0.0
        t1 = edge.values[7] if isinstance(edge.values[7], float) else 1.0
        curve = self._ref_at(edge, EDGE_CURVE)
        kind = curve.name if curve is not None else "straight"

        pts = None
        if kind == "ellipse":
            pts = self._sample_ellipse(curve, t0, t1)
        elif kind in ("intcurve", "curve"):
            pts = _sample_nubs(curve, self._arc_segments)

        if not pts:  # straight or unevaluated curve -> just the endpoints
            pts = [p for p in (v0, v1) if p is not None]
        else:
            # snap sampled ends onto the exact vertices to avoid tiny gaps
            if v0 is not None:
                pts[0] = v0
            if v1 is not None:
                pts[-1] = v1

        if reverse:
            pts = list(reversed(pts))
        return pts

    def _loop_ring(self, loop):
        """Ordered ring of 3D points around one loop (coedge chain, sampled)."""
        coedges = self._walk(self._ref_at(loop, LOOP_COEDGE), COEDGE_NEXT)
        ring = []
        for ce in coedges:
            edge = self._ref_at(ce, COEDGE_EDGE)
            if edge is None:
                continue
            sense = ce.values[COEDGE_SENSE] if COEDGE_SENSE < len(ce.values) else False
            poly = self._sample_edge(edge, reverse=(sense is True))
            # append everything but the last point (the next coedge repeats it)
            ring.extend(poly[:-1] if len(poly) > 1 else poly)
        return ring

    def _face(self, face_rec):
        surf = self._ref_at(face_rec, FACE_SURFACE)
        kind, normal = self._surface_info(surf)
        loops = self._walk(self._ref_at(face_rec, FACE_LOOP), LOOP_NEXT)
        rings = []
        for lp in loops:
            ring = self._loop_ring(lp)
            if len(ring) >= 3:
                rings.append(ring)
        return Face(surface_kind=kind, loops=rings, normal=normal)

    # -- public -------------------------------------------------------------
    def bodies(self):
        """Return one :class:`Body` per solid, grouping faces by owner.

        Rather than chase lump/shell *next-sibling* pointers (whose slots vary),
        every face is walked up its verified owner chain
        ``face -> shell -> lump -> body`` and grouped by the owning body.  This
        needs only pointers that were confirmed against the sample.
        """
        groups: dict = {}
        order: list = []
        for face_rec in self.f.by_type("face"):
            shell = self._ref_at(face_rec, FACE_SHELL)
            lump = self._ref_at(shell, SHELL_LUMP)
            body = self._ref_at(lump, LUMP_BODY)
            key = body.index if body is not None else -1
            if key not in groups:
                groups[key] = []
                order.append(key)
            face = self._face(face_rec)
            if face.loops:
                groups[key].append(face)
        return [Body(faces=groups[k]) for k in order if groups[k]]
