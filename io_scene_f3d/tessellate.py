"""Triangulate B-rep faces into meshes (pure Python, no Blender dependency).

For each face we have one or more boundary loops as rings of 3D points (outer
loop first, holes after).  We project them onto the face plane, run an
ear-clipping triangulation (with hole bridging), and lift the resulting
triangles back to 3D.

This is exact for planar faces.  Curved faces are approximated by the polygon of
their (optionally curve-sampled) boundary edges — good enough to visualise the
solid; finer surface tessellation is future work.
"""

from __future__ import annotations

import math

Vec = tuple


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    m = math.sqrt(_dot(a, a))
    if m < 1e-30:
        return (0.0, 0.0, 0.0)
    return (a[0] / m, a[1] / m, a[2] / m)


def newell_normal(pts) -> Vec:
    """Robust polygon normal via Newell's method."""
    nx = ny = nz = 0.0
    n = len(pts)
    for i in range(n):
        a = pts[i]
        b = pts[(i + 1) % n]
        nx += (a[1] - b[1]) * (a[2] + b[2])
        ny += (a[2] - b[2]) * (a[0] + b[0])
        nz += (a[0] - b[0]) * (a[1] + b[1])
    return _norm((nx, ny, nz))


def _basis(normal):
    n = _norm(normal)
    if n == (0.0, 0.0, 0.0):
        n = (0.0, 0.0, 1.0)
    a = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _norm(_cross(a, n))
    v = _cross(n, u)
    return u, v, n


def _area2(p2d):
    """Signed area*2 of a 2D polygon (positive == CCW)."""
    s = 0.0
    n = len(p2d)
    for i in range(n):
        x0, y0 = p2d[i]
        x1, y1 = p2d[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return s


def _point_in_tri(p, a, b, c):
    d = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
    if abs(d) < 1e-30:
        return False
    a1 = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / d
    a2 = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / d
    a3 = 1.0 - a1 - a2
    return a1 >= -1e-12 and a2 >= -1e-12 and a3 >= -1e-12


def _ear_clip(poly2d):
    """Ear-clipping triangulation of a simple polygon.

    ``poly2d`` is a list of (x, y).  Returns a list of (i, j, k) index triples
    into ``poly2d``.  The polygon is made CCW internally.
    """
    n = len(poly2d)
    if n < 3:
        return []
    idx = list(range(n))
    if _area2(poly2d) < 0:
        idx.reverse()

    def is_convex(i0, i1, i2):
        a, b, c = poly2d[i0], poly2d[i1], poly2d[i2]
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0

    tris = []
    guard = 0
    while len(idx) > 3 and guard < 10000:
        guard += 1
        clipped = False
        m = len(idx)
        for k in range(m):
            i0, i1, i2 = idx[(k - 1) % m], idx[k], idx[(k + 1) % m]
            if not is_convex(i0, i1, i2):
                continue
            a, b, c = poly2d[i0], poly2d[i1], poly2d[i2]
            if any(
                j not in (i0, i1, i2) and _point_in_tri(poly2d[j], a, b, c)
                for j in idx
            ):
                continue
            tris.append((i0, i1, i2))
            idx.pop(k)
            clipped = True
            break
        if not clipped:
            break  # degenerate; stop with what we have
    if len(idx) == 3:
        tris.append((idx[0], idx[1], idx[2]))
    return tris


def _bridge_holes(outer2d, holes2d):
    """Merge holes into the outer polygon by bridge edges.

    Returns (merged2d, mapping) where mapping[i] gives the original
    (loop_index, vertex_index) for each merged vertex.  Uses Eberly's
    right-most-vertex bridging; falls back to ignoring a hole it cannot bridge.
    """
    # start with outer, oriented CCW
    outer = list(range(len(outer2d)))
    pts = list(outer2d)
    mapping = [(0, i) for i in range(len(outer2d))]
    ring = outer  # indices into pts, current merged polygon order
    merged = list(range(len(pts)))

    def order_ccw(indices):
        poly = [pts[i] for i in indices]
        if _area2(poly) < 0:
            indices = list(reversed(indices))
        return indices

    merged = order_ccw(merged)

    for h, hole in enumerate(holes2d, start=1):
        base = len(pts)
        hpts = list(hole)
        # holes wind opposite to outer
        hidx = list(range(len(hpts)))
        if _area2(hpts) > 0:
            hidx.reverse()
        for j in hidx:
            pts.append(hpts[j])
            mapping.append((h, j))
        hole_indices = [base + off for off in range(len(hidx))]
        # rightmost hole vertex
        m_pos = max(range(len(hidx)), key=lambda t: pts[base + t][0])
        m_idx = base + m_pos
        # connect to the nearest merged-polygon vertex (simple, robust-ish)
        best = min(
            range(len(merged)),
            key=lambda t: (pts[merged[t]][0] - pts[m_idx][0]) ** 2
            + (pts[merged[t]][1] - pts[m_idx][1]) ** 2,
        )
        # splice: outer ... best, hole(rotated to start at m), m, best ...
        rot = hole_indices[m_pos:] + hole_indices[:m_pos]
        merged = (
            merged[: best + 1]
            + rot
            + [rot[0], merged[best]]
            + merged[best + 1:]
        )
    return pts, merged, mapping


def tessellate_face(face):
    """Triangulate one :class:`brep.Face`.

    Returns ``(verts3d, tris)`` where ``verts3d`` is a list of (x, y, z) and
    ``tris`` a list of (a, b, c) indices into it.
    """
    loops = face.loops
    if not loops:
        return [], []
    outer3d = loops[0]
    normal = face.normal or newell_normal(outer3d)
    u, v, _ = _basis(normal)
    origin = outer3d[0]

    def to2d(p):
        d = _sub(p, origin)
        return (_dot(d, u), _dot(d, v))

    if len(loops) == 1:
        # Orient the ring CCW in 2D so ear-clip indices line up with verts3d.
        idx = list(range(len(outer3d)))
        if _area2([to2d(p) for p in outer3d]) < 0:
            idx.reverse()
        verts3d = [outer3d[i] for i in idx]
        tris = _ear_clip([to2d(p) for p in verts3d])
        return verts3d, tris

    outer2d = [to2d(p) for p in outer3d]
    holes2d = [[to2d(p) for p in lp] for lp in loops[1:]]
    pts2d, merged, mapping = _bridge_holes(outer2d, holes2d)
    merged_pts2d = [pts2d[i] for i in merged]
    local_tris = _ear_clip(merged_pts2d)
    # build 3d verts in merged order
    verts3d = []
    for i in merged:
        loop_i, vtx_i = mapping[i]
        verts3d.append(loops[loop_i][vtx_i])
    return verts3d, local_tris
