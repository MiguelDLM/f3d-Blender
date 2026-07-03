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


def _arclen_fractions(pts, closed):
    """Cumulative arc-length fraction of each point along a chain/ring."""
    n = len(pts)
    seg = [_dist(pts[i], pts[i + 1]) for i in range(n - 1)]
    if closed:
        seg.append(_dist(pts[-1], pts[0]))
    total = sum(seg) or 1.0
    fr = [0.0]
    acc = 0.0
    for s in seg[:len(pts) - 1]:
        acc += s
        fr.append(acc / total)
    return fr


def _bridge(chain_a, chain_b, closed):
    """Triangulate between two chains WITHOUT moving any point.

    Walks both point lists by arc-length fraction (like Blender's bridge
    loops), so the original vertices — shared exactly with neighbouring
    faces — are preserved and the final mesh welds watertight.
    """
    fa = _arclen_fractions(chain_a, closed)
    fb = _arclen_fractions(chain_b, closed)
    na, nb = len(chain_a), len(chain_b)
    verts = list(chain_a) + list(chain_b)
    tris = []
    ia = ib = 0
    steps = na + nb + 2
    for _ in range(steps):
        a_done = ia >= na - (0 if closed else 1)
        b_done = ib >= nb - (0 if closed else 1)
        if a_done and b_done:
            break
        next_fa = fa[ia + 1] if ia + 1 < na else 1.0
        next_fb = fb[ib + 1] if ib + 1 < nb else 1.0
        if not a_done and (b_done or next_fa <= next_fb):
            a2 = (ia + 1) % na
            tris.append((ia, a2, na + ib))
            ia += 1
        else:
            b2 = (ib + 1) % nb
            tris.append((ia % na, na + b2, na + ib))
            ib += 1
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
    # Case A: two loops -> bridge between them directly (no axis needed).
    if len(face.loops) == 2:
        a, b = face.loops[0], face.loops[1]
        if len(a) < 3 or len(b) < 3:
            return None
        return _bridge(a, _phase_align(a, b), closed=True)

    # Case B: single loop split by two parallel straight rails.
    if len(face.loops) == 1 and face.loop_edges:
        edges = face.loop_edges[0]
        if len(edges) < 3:
            return None
        # Seam-cut annulus (e.g. a fillet ring): the same edge appears twice
        # (once per direction); the remaining edges form the two closed
        # profiles -> loft them ring-to-ring.  The profiles are rebuilt by
        # endpoint chaining rather than by their position between the seam
        # copies, because greedy loop chaining can put the two copies next to
        # each other (their endpoints coincide exactly).
        seam = _find_seam_pair(edges)
        if seam:
            i, j = seam
            rest = [e for k, e in enumerate(edges) if k not in (i, j)]
            comps = _closed_components(rest)
            if len(comps) == 2 and all(len(c) >= 3 for c in comps):
                prof1, prof2 = comps
                strip = _ring_strip(prof1, prof2, closed=True)
                if strip:
                    return strip
                return _bridge(prof1, _phase_align(prof1, prof2), closed=True)
        straights = [k for k, e in enumerate(edges) if _is_straight(e)]
        rails = _pick_rails(edges, straights, face.axis)
        if not rails:
            return None
        i, j = rails
        prof1 = [p for e in edges[i + 1:j] for p in e]
        prof2 = [p for e in (edges[j + 1:] + edges[:i]) for p in e]
        if len(prof1) < 2 or len(prof2) < 2:
            return None
        prof2 = list(reversed(prof2))     # opposite traversal
        return _bridge(prof1, prof2, closed=False)

    return None


def _pt_seg_dist(p, a, b):
    """Distance from ``p`` to segment ``ab`` and the closest point."""
    ab = _sub(b, a)
    L2 = _dot(ab, ab)
    t = 0.0 if L2 < 1e-18 else max(0.0, min(1.0, _dot(_sub(p, a), ab) / L2))
    q = tuple(a[c] + ab[c] * t for c in range(3))
    return _dist(p, q), q


def _pt_polyline_nearest(p, poly, closed):
    """Nearest point on a polyline (segment-interpolated)."""
    best = (float("inf"), poly[0])
    n = len(poly)
    last = n if closed else n - 1
    for i in range(last):
        d, q = _pt_seg_dist(p, poly[i], poly[(i + 1) % n])
        if d < best[0]:
            best = (d, q)
    return best


def _ring_strip(prof_a, prof_b, closed):
    """Quad-strip between two equal-count vertex-paired profiles.

    Used for seam-cut blend rings whose two contact profiles were sampled
    from the same periodic curve family, so their vertices correspond 1:1
    once phase-aligned.  Pairing vertices directly (instead of walking by
    arc length as :func:`_bridge` does) keeps the strip untwisted, which
    matters on blends: Fusion evaluates their cross-sections as shallow
    spline sections (NOT circular rolling-ball arcs -- verified against the
    reference OBJ, whose sections deviate from a straight chord by ~5% of
    its length), so the honest approximation is the ruled strip between
    corresponding points.  Returns ``None`` when the counts differ.
    """
    if len(prof_a) != len(prof_b):
        return None
    if closed:
        prof_b = _phase_align(prof_a, prof_b)
    elif _dist(prof_a[0], prof_b[0]) > _dist(prof_a[0], prof_b[-1]):
        prof_b = list(reversed(prof_b))
    n = len(prof_a)
    verts = list(prof_a) + list(prof_b)
    tris = []
    last = n if closed else n - 1
    for i in range(last):
        j = (i + 1) % n
        tris.append((i, j, n + j))
        tris.append((i, n + j, n + i))
    return verts, tris


def _closed_components(polys):
    """Chain edge polylines into closed rings by shared-endpoint matching.

    Loop edges share exact vertex coordinates, so matching at 1e-9 is safe.
    Each returned ring has its duplicate closing point dropped.
    """
    remaining = [list(p) for p in polys]
    comps = []
    while remaining:
        chain = remaining.pop(0)
        grown = True
        while grown and _dist(chain[0], chain[-1]) > 1e-9:
            grown = False
            for k, poly in enumerate(remaining):
                if _dist(chain[-1], poly[0]) < 1e-9:
                    chain += poly[1:]
                elif _dist(chain[-1], poly[-1]) < 1e-9:
                    chain += list(reversed(poly))[1:]
                else:
                    continue
                remaining.pop(k)
                grown = True
                break
        if len(chain) > 1 and _dist(chain[0], chain[-1]) < 1e-9:
            chain.pop()
        comps.append(chain)
    return comps


def _find_seam_pair(edges):
    """Indices (i, j) of two polylines with identical (or reversed) geometry."""
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            a, b = edges[i], edges[j]
            if len(a) != len(b):
                continue
            if all(_dist(p, q) < 1e-9 for p, q in zip(a, b)) or \
               all(_dist(p, q) < 1e-9 for p, q in zip(a, reversed(b))):
                return (i, j)
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


def _phase_align(ref, ring):
    """Rotate/flip closed ``ring`` (without moving points) to match ``ref``.

    Chooses the rotation starting nearest ``ref[0]`` and the direction whose
    quarter-way point matches ``ref``'s quarter-way point best, so a
    subsequent arc-length bridge walks both rings in the same sense.
    """
    n = len(ring)
    start = min(range(n), key=lambda i: _dist(ring[i], ref[0]))
    fwd = ring[start:] + ring[:start]
    rev = list(reversed(fwd))
    rev = [rev[-1]] + rev[:-1]          # keep the matched start point first
    qr = ref[len(ref) // 4]
    if _dist(fwd[n // 4], qr) <= _dist(rev[n // 4], qr):
        return fwd
    return rev


def _wrap_pi(a):
    """Wrap an angle difference into (-pi, pi]."""
    return (a + math.pi) % (2 * math.pi) - math.pi


def _point_in_poly(x, y, poly):
    """Even-odd point-in-polygon test in 2D."""
    inside = False
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xt = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if xt > x:
                inside = not inside
    return inside


def _cut_full_wrap(surf, uv_loops, rings3, w_u, full):
    """Join two full-period rings into one simple uv polygon via a seam.

    Ring A is traversed with increasing u, ring B (rotated to start where A
    ends, one period up) with decreasing u; the two are connected by a right
    seam chain and its exact copy one period down on the left, so both seam
    sides evaluate to identical 3D points and the mesh welds shut.
    Returns ``(outer_uv, outer_3d)`` with the lists aligned point-for-point.
    """
    two_pi = 2 * math.pi
    ia, ib = full
    A = list(uv_loops[ia])
    B = list(uv_loops[ib])
    A3 = list(rings3[ia])
    B3 = list(rings3[ib])
    if w_u[ia] < 0:
        A.reverse()
        A3.reverse()
    if w_u[ib] > 0:
        B.reverse()
        B3.reverse()
    # cut where no hole is: rotate A to start at the boundary point whose u
    # (mod 2*pi) is furthest from every hole sample, so the seam line does
    # not cross a hole ring (which would make the outer polygon invalid)
    hole_us = [u for i, uv in enumerate(uv_loops) if i not in full
               for (u, _v) in uv]
    if hole_us:
        def clearance(u):
            return min(abs(_wrap_pi(u - hu)) for hu in hole_us)
        j = max(range(len(A)), key=lambda i: clearance(A[i][0]))
        if j:
            A = A[j:] + [(u + two_pi, v) for (u, v) in A[:j]]
            A3 = A3[j:] + A3[:j]
    # bring B onto A's period branch
    mu_a = sum(p[0] for p in A) / len(A)
    mu_b = sum(p[0] for p in B) / len(B)
    k = round((mu_a - mu_b) / two_pi)
    if k:
        B = [(u + two_pi * k, v) for (u, v) in B]
    # rotate B (winds downward) so it starts near A's end, one period up
    target = A[0][0] + two_pi
    j = min(range(len(B)), key=lambda i: abs(B[i][0] - target))
    B = B[j:] + [(u - two_pi, v) for (u, v) in B[:j]]
    B3 = B3[j:] + B3[:j]

    a_end = (A[0][0] + two_pi, A[0][1])     # same 3D point as A[0]
    b_end = (B[0][0] - two_pi, B[0][1])     # same 3D point as B[0]

    # seam samples at ring-A boundary density; evaluated once and shared by
    # both seam sides so they weld exactly
    seg = sum(_dist(A3[i], A3[i + 1]) for i in range(len(A3) - 1))
    seg /= max(1, len(A3) - 1)
    span = _dist(surf.eval(*a_end), surf.eval(*B[0]))
    m = int(min(200, span / max(seg, 1e-9)))
    seam_uv = []
    seam_3d = []
    for t in range(1, m + 1):
        f = t / (m + 1)
        uv = (a_end[0] + (B[0][0] - a_end[0]) * f,
              a_end[1] + (B[0][1] - a_end[1]) * f)
        seam_uv.append(uv)
        seam_3d.append(surf.eval(*uv))

    outer_uv = (A + [a_end] + seam_uv + B + [b_end]
                + [(u - two_pi, v) for (u, v) in reversed(seam_uv)])
    outer_3d = (A3 + [A3[0]] + seam_3d + B3 + [B3[0]]
                + list(reversed(seam_3d)))
    return outer_uv, outer_3d


def _tessellate_on_surface(face):
    """Exact tessellation of a face with an evaluable surface.

    The boundary loops are mapped into the surface's (u, v) parameter space,
    triangulated there with Blender's constrained Delaunay (plus a grid of
    interior Steiner points at boundary-sampling density), and the interior
    vertices are lifted back with ``surface.eval``.  Boundary vertices reuse
    the exact 3D ring points so adjacent faces weld watertight.

    Returns ``None`` when not applicable (no surface, full period wrap, or
    running without Blender's ``mathutils``).
    """
    surf = getattr(face, "surface", None)
    if surf is None:
        return None
    try:
        from mathutils.geometry import delaunay_2d_cdt
        from mathutils import Vector
    except ImportError:
        return None

    # --- project rings to (u, v), unwrapping periodic u/v along the loop ---
    periodic_v = getattr(surf, "periodic_v", False)
    uv_loops = []
    w_u = []
    for ring in face.loops:
        uv = []
        prev = None
        for p in ring:
            u, v = surf.project(p)
            if prev is not None:
                if surf.periodic_u:
                    u = prev[0] + _wrap_pi(u - prev[0])
                if periodic_v:
                    v = prev[1] + _wrap_pi(v - prev[1])
            uv.append((u, v))
            prev = (u, v)
        wu = wv = 0.0
        if len(uv) > 1:
            if surf.periodic_u:
                wu = (uv[-1][0] + _wrap_pi(uv[0][0] - uv[-1][0])) - uv[0][0]
            if periodic_v:
                wv = (uv[-1][1] + _wrap_pi(uv[0][1] - uv[-1][1])) - uv[0][1]
        if abs(wv) > math.pi:
            return None      # loop wraps the full v period -> loft instead
        uv_loops.append(uv)
        w_u.append(wu)

    # Rings unwrap independently, so a hole may land one period away from
    # the outer loop; shift each ring by whole periods onto the first ring.
    if len(uv_loops) > 1:
        ref_u = sum(p[0] for p in uv_loops[0]) / len(uv_loops[0])
        ref_v = sum(p[1] for p in uv_loops[0]) / len(uv_loops[0])
        for li in range(1, len(uv_loops)):
            uv = uv_loops[li]
            du = dv = 0.0
            if surf.periodic_u:
                mu = sum(p[0] for p in uv) / len(uv)
                du = 2 * math.pi * round((ref_u - mu) / (2 * math.pi))
            if periodic_v:
                mv = sum(p[1] for p in uv) / len(uv)
                dv = 2 * math.pi * round((ref_v - mv) / (2 * math.pi))
            if du or dv:
                uv_loops[li] = [(u + du, v + dv) for (u, v) in uv]

    # --- full-period faces (a tube wall): cut the seam open ---
    # A face bounded by two rings that each wrap the whole u period (e.g. a
    # cylinder wall with its two end circles, possibly with holes where other
    # parts join) has no simple uv polygon.  Cut it at one u value: join the
    # two rings with two seam chains that share the same 3D points, so the
    # mesh welds shut across the cut.
    full = [i for i, w in enumerate(w_u) if abs(w) > math.pi]
    if full:
        if not surf.periodic_u or len(full) != 2:
            return None      # single seam-cut rings are handled by the loft
        cut = _cut_full_wrap(surf, uv_loops, face.loops, w_u, full)
        if cut is None:
            return None
        outer_uv, outer_3d = cut
        loops_uv = [outer_uv]
        loops_3d = [outer_3d]
        dom_u = sum(p[0] for p in outer_uv) / len(outer_uv)
        for i in range(len(uv_loops)):
            if i in full:
                continue
            uv = uv_loops[i]
            mu = sum(p[0] for p in uv) / len(uv)
            k = round((dom_u - mu) / (2 * math.pi))
            loops_uv.append([(u + 2 * math.pi * k, v) for (u, v) in uv])
            loops_3d.append(list(face.loops[i]))
    else:
        loops_uv = uv_loops
        loops_3d = [list(r) for r in face.loops]

    # --- scale u and v to be roughly isometric with 3D distance ---
    def _axis_scale(k):
        ratios = []
        for uv, ring in zip(loops_uv, loops_3d):
            for i in range(len(uv) - 1):
                dp = abs(uv[i + 1][k] - uv[i][k])
                d3 = _dist(ring[i], ring[i + 1])
                if dp > 1e-9 and d3 > 1e-12:
                    ratios.append(d3 / dp)
        s = sorted(ratios)[len(ratios) // 2] if ratios else 1.0
        return max(s, 1e-6)

    su = _axis_scale(0)
    sv = _axis_scale(1) if periodic_v else 1.0

    pts2 = []
    exact3d = []
    edges = []
    polys2 = []
    for uv, ring in zip(loops_uv, loops_3d):
        base = len(pts2)
        n = len(uv)
        pts2.extend(Vector((u * su, v * sv)) for (u, v) in uv)
        exact3d.extend(ring)
        edges.extend((base + i, base + (i + 1) % n) for i in range(n))
        polys2.append([(u * su, v * sv) for (u, v) in uv])
    n_boundary = len(pts2)

    # --- interior Steiner grid at boundary sampling density ---
    seg = sorted((pts2[a] - pts2[b]).length for a, b in edges)
    h = max(seg[len(seg) // 2] if seg else 0.1, 1e-4)
    xs = [p.x for p in pts2]
    ys = [p.y for p in pts2]
    # keep the interior grid bounded (~1500 points): visual smoothness is set
    # by the boundary sampling; the interior only needs comparable density
    area_box = (max(xs) - min(xs)) * (max(ys) - min(ys))
    while area_box / (h * h) > 1500:
        h *= 1.5
    grid = []
    if True:
        x = min(xs) + h * 0.5
        while x < max(xs):
            y = min(ys) + h * 0.5
            while y < max(ys):
                if _point_in_poly(x, y, polys2[0]) and not any(
                    _point_in_poly(x, y, hp) for hp in polys2[1:]
                ):
                    grid.append(Vector((x, y)))
                y += h
            x += h

    try:
        vco, _e, ofaces, overts, _oe, _of = delaunay_2d_cdt(
            pts2 + grid, edges, [], 1, 1e-6
        )
    except Exception:
        return None

    verts3 = []
    for i, co in enumerate(vco):
        orig = [k for k in overts[i] if k < n_boundary]
        if orig:
            verts3.append(exact3d[orig[0]])
        else:
            verts3.append(surf.eval(co.x / su, co.y / sv))

    tris = []
    for fverts in ofaces:
        for k in range(1, len(fverts) - 1):
            a, b, c = fverts[0], fverts[k], fverts[k + 1]
            # drop triangles whose centroid lies outside the trimmed region
            cx = (vco[a].x + vco[b].x + vco[c].x) / 3
            cy = (vco[a].y + vco[b].y + vco[c].y) / 3
            if _point_in_poly(cx, cy, polys2[0]) and not any(
                _point_in_poly(cx, cy, hp) for hp in polys2[1:]
            ):
                tris.append((a, b, c))
    return (verts3, tris) if tris else None


def tessellate_face(face):
    """Triangulate one :class:`brep.Face`.

    Returns ``(verts3d, tris)`` where ``verts3d`` is a list of (x, y, z) and
    ``tris`` a list of (a, b, c) indices into it.
    """
    loops = face.loops
    if not loops:
        return [], []

    if face.surface_kind != "plane":
        exact = _tessellate_on_surface(face)
        if exact:
            return exact
        # fallback: blends (no evaluable surface), faces wrapping the full
        # period (hole walls -> exact ring-to-ring loft), or no Blender
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
