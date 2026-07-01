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

from dataclasses import dataclass, field

try:  # package import (inside Blender add-on)
    from . import sab
except ImportError:  # standalone (tests)
    import sab

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


@dataclass
class Face:
    surface_kind: str                       # plane / cone / spline / ...
    loops: list = field(default_factory=list)   # list[list[(x,y,z)]]; [0]=outer
    normal: tuple | None = None             # analytic surface normal, if known


@dataclass
class Body:
    faces: list = field(default_factory=list)


class Brep:
    def __init__(self, sabfile: "sab.SabFile"):
        self.f = sabfile

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

    def _loop_ring(self, loop):
        """Ordered list of 3D points around one loop (via coedge chain)."""
        first_ce = self._ref_at(loop, LOOP_COEDGE)
        coedges = self._walk(first_ce, COEDGE_NEXT)
        ring = []
        for ce in coedges:
            edge = self._ref_at(ce, COEDGE_EDGE)
            if edge is None:
                continue
            v0 = self._ref_at(edge, EDGE_VSTART)
            v1 = self._ref_at(edge, EDGE_VEND)
            sense = ce.values[COEDGE_SENSE] if COEDGE_SENSE < len(ce.values) else False
            start_v = v1 if sense is True else v0
            p = self._vertex_point(start_v)
            if p is not None:
                ring.append(p)
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
