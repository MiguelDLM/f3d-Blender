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


def _resample(poly, n):
    """Resample a polyline to exactly ``n`` points by arc length."""
    if len(poly) == 1:
        return [poly[0]] * n
    seglen = [_dist(poly[i], poly[i + 1]) for i in range(len(poly) - 1)]
    total = sum(seglen)
    if total < 1e-12:
        return [poly[0]] * n
    out = [poly[0]]
    step = total / (n - 1)
    d = 0.0
    i = 0
    acc = 0.0
    for k in range(1, n - 1):
        target = k * step
        while i < len(seglen) and acc + seglen[i] < target:
            acc += seglen[i]
            i += 1
        if i >= len(seglen):
            out.append(poly[-1])
            continue
        t = (target - acc) / seglen[i] if seglen[i] > 1e-12 else 0.0
        a, b = poly[i], poly[i + 1]
        out.append(tuple(a[c] + (b[c] - a[c]) * t for c in range(3)))
    out.append(poly[-1])
    return out


def _dist(a, b):
    return math.sqrt(sum((a[c] - b[c]) ** 2 for c in range(3)))


def _loft(rail_a, rail_b):
    """Triangulate the strip between two equal-length point lists."""
    n = len(rail_a)
    verts = list(rail_a) + list(rail_b)
    tris = []
    for i in range(n - 1):
        a0, a1 = i, i + 1
        b0, b1 = n + i, n + i + 1
        tris.append((a0, a1, b1))
        tris.append((a0, b1, b0))
    return verts, tris


def _tessellate_swept(face):
    """Tessellate a swept/extruded/ruled face by lofting between two profiles.

    Handles a face with two boundary loops (e.g. a cylindrical hole wall
    between two circles) and a single-loop face split by two straight,
    mutually-parallel rail edges into two profile chains (e.g. the rod cradle
    or a fillet).  The sweep direction is taken from ``face.axis`` when known,
    otherwise inferred from the rails, so it also works for shared (``ref``)
    surfaces that carry no inline geometry.  Returns ``None`` if the face does
    not match this pattern.
    """
    # Case A: two loops -> loft between them directly (no axis needed).
    if len(face.loops) == 2:
        a, b = face.loops[0], face.loops[1]
        if len(a) < 3 or len(b) < 3:
            return None
        n = max(len(a), len(b), 8)
        ra = _resample(a + [a[0]], n)
        rb = _align_ring(ra, _resample(b + [b[0]], n))
        return _loft(ra, rb)

    # Case B: single loop split by two parallel straight rails.
    if len(face.loops) == 1 and face.loop_edges:
        edges = face.loop_edges[0]
        if len(edges) < 3:
            return None
        straights = [k for k, e in enumerate(edges) if _is_straight(e)]
        rails = _pick_rails(edges, straights, face.axis)
        if not rails:
            return None
        i, j = rails
        prof1 = [p for e in edges[i + 1:j] for p in e]
        prof2 = [p for e in (edges[j + 1:] + edges[:i]) for p in e]
        if len(prof1) < 2 or len(prof2) < 2:
            return None
        n = max(len(prof1), len(prof2), 6)
        ra = _resample(prof1, n)
        rb = _resample(list(reversed(prof2)), n)   # opposite traversal
        return _loft(ra, rb)

    return None


def _edge_dir(e):
    d = _sub(e[-1], e[0])
    L = math.sqrt(_dot(d, d))
    return (d[0] / L, d[1] / L, d[2] / L) if L > 1e-12 else None


def _is_straight(e):
    L = _dist(e[0], e[-1])
    if L < 1e-9:
        return False
    return all(_point_line_dist(p, e[0], e[-1]) < 0.02 * L for p in e)


def _pick_rails(edges, straights, axis):
    """Choose the two rail edges (parallel, opposite sides of the loop)."""
    if len(straights) < 2:
        return None
    dirs = {k: _edge_dir(edges[k]) for k in straights}
    axis_u = _norm(axis) if axis else None
    best = None
    for a in range(len(straights)):
        for b in range(a + 1, len(straights)):
            ka, kb = straights[a], straights[b]
            da, db = dirs[ka], dirs[kb]
            if da is None or db is None:
                continue
            par = abs(_dot(da, db))            # parallel rails
            if par < 0.95:
                continue
            score = par
            if axis_u is not None:
                score += abs(_dot(da, axis_u))
            # prefer rails that are on opposite sides (indices spread out)
            if best is None or score > best[0]:
                best = (score, ka, kb)
    if best is None:
        return None
    _, i, j = best
    return (min(i, j), max(i, j))


def _point_line_dist(p, a, b):
    ab = _sub(b, a)
    L2 = _dot(ab, ab)
    if L2 < 1e-18:
        return _dist(p, a)
    t = _dot(_sub(p, a), ab) / L2
    proj = tuple(a[c] + ab[c] * t for c in range(3))
    return _dist(p, proj)


def _align_ring(ref, ring):
    """Rotate/flip closed ``ring`` to best match ``ref`` point order."""
    n = len(ring)

    def cost(seq):
        return sum(_dist(ref[i], seq[i]) for i in range(n))

    best = ring
    best_c = cost(ring)
    rev = list(reversed(ring))
    for cand in (ring, rev):
        for s in range(n):
            rot = cand[s:] + cand[:s]
            c = cost(rot)
            if c < best_c:
                best_c = c
                best = rot
    return best


def tessellate_face(face):
    """Triangulate one :class:`brep.Face`.

    Returns ``(verts3d, tris)`` where ``verts3d`` is a list of (x, y, z) and
    ``tris`` a list of (a, b, c) indices into it.
    """
    loops = face.loops
    if not loops:
        return [], []

    # Curved faces (anything but a plane) are lofted between their profiles;
    # fall back to planar triangulation when that doesn't apply.
    if face.surface_kind != "plane" or len(loops) >= 2:
        swept = _tessellate_swept(face)
        if swept and swept[1]:
            return swept

    return triangulate_planar(loops, face.normal)


def _mathutils_tessellate(loops2d):
    """Use Blender's robust polygon tessellator if available; else ``None``."""
    try:
        from mathutils.geometry import tessellate_polygon
        from mathutils import Vector
    except ImportError:
        return None
    chains = [[Vector((x, y, 0.0)) for (x, y) in lp] for lp in loops2d]
    return tessellate_polygon(chains)


def triangulate_planar(loops, normal=None):
    """Triangulate planar loops (outer + holes); returns (verts3d, tris).

    Uses Blender's ``tessellate_polygon`` (robust, hole-aware) when running
    inside Blender, and falls back to an ear-clip + hole-bridge otherwise so
    the standalone tests still work.
    """
    outer3d = loops[0]
    n = normal or newell_normal(outer3d)
    u, v, _ = _basis(n)
    origin = outer3d[0]

    def to2d(p):
        d = _sub(p, origin)
        return (_dot(d, u), _dot(d, v))

    # concatenated vertex list across all loops, with 2D projections per loop
    verts3d = [p for lp in loops for p in lp]
    loops2d = [[to2d(p) for p in lp] for lp in loops]

    tris = _mathutils_tessellate(loops2d)
    if tris is not None:
        # tessellate_polygon indexes into the concatenated point list
        offsets = []
        off = 0
        for lp in loops:
            offsets.append(off)
            off += len(lp)
        flat = [(offsets[li] + vi) for li, lp in enumerate(loops2d) for vi in range(len(lp))]
        return verts3d, [tuple(flat[i] for i in t) for t in tris]

    # --- pure-Python fallback (tests) ---
    if len(loops) == 1:
        idx = list(range(len(outer3d)))
        if _area2([to2d(p) for p in outer3d]) < 0:
            idx.reverse()
        vv = [outer3d[i] for i in idx]
        return vv, _ear_clip([to2d(p) for p in vv])

    outer2d = [to2d(p) for p in outer3d]
    holes2d = [[to2d(p) for p in lp] for lp in loops[1:]]
    pts2d, merged, mapping = _bridge_holes(outer2d, holes2d)
    tris = _ear_clip([pts2d[i] for i in merged])
    vv = [loops[mapping[i][0]][mapping[i][1]] for i in merged]
    return vv, tris
