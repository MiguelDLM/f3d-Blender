"""Blender-side glue: turn a ``.f3d`` file into mesh objects.

Pipeline: :mod:`container` (unzip) -> :mod:`sab` (tokenize ASM) -> :mod:`brep`
(topology + geometry) -> :mod:`tessellate` (triangles) -> Blender ``bmesh``.

Only this module imports ``bpy``/``bmesh``; everything upstream is plain Python.
Fusion stores model coordinates in centimetres, so the default ``global_scale``
of 0.01 maps them to Blender metres (1 model cm -> 0.01 m).
"""

from __future__ import annotations

import os

import bpy
import bmesh
import mathutils

try:
    from . import container, sab, brep, tessellate
except ImportError:  # pragma: no cover - allow flat import
    import container, sab, brep, tessellate


def _weld_t_junctions(bm, tol):
    """Split boundary edges at boundary vertices lying on them, then weld.

    ASM stores the same model ring as different edge records on each side of
    some face pairs, sampled at different densities; after welding, the odd
    vertices of one side sit ON (or within chord error of) the other side's
    boundary edges, leaving hairline cracks.  Splitting the edge at the
    stray vertex and snapping the new vertex onto it closes the crack
    exactly instead of papering over it with fill faces.
    """
    border = [e for e in bm.edges if e.is_boundary]
    if not border:
        return
    bverts = list({v for e in border for v in e.verts})
    kd = mathutils.kdtree.KDTree(len(bverts))
    for i, v in enumerate(bverts):
        kd.insert(v.co, i)
    kd.balance()
    for e in border:
        if not e.is_valid:
            continue
        hits = []
        a = e.verts[0].co.copy()
        b = e.verts[1].co.copy()
        ab = b - a
        length2 = ab.length_squared
        if length2 < 1e-16:
            continue
        for (co, idx, _d) in kd.find_range((a + b) / 2, ab.length / 2 + tol):
            v = bverts[idx]
            if not v.is_valid or v in e.verts:
                continue
            t = (v.co - a).dot(ab) / length2
            if t < 1e-3 or t > 1.0 - 1e-3:
                continue
            if ((a + ab * t) - v.co).length <= tol:
                hits.append((t, v))
        if not hits:
            continue
        # sort by t only -- ties must not fall through to BMVert comparison
        hits.sort(key=lambda h: h[0], reverse=True)  # far end first
        cur = e
        v0 = e.verts[0]
        prev_t = 1.0
        for t, v in hits:
            # fac is a fraction (from v0) of the remaining [0, prev_t] span
            # of the original edge; after the split, continue on whichever
            # half still contains v0
            new_edge, new_vert = bmesh.utils.edge_split(cur, v0, t / prev_t)
            new_vert.co = v.co      # exact weld target
            cur = new_edge if v0 in new_edge.verts else cur
            if v0 not in cur.verts:
                break
            prev_t = t
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-6)


def _build_mesh(name, body, scale, chord_tol=0.0):
    """Build a welded, normal-consistent mesh from one :class:`brep.Body`."""
    verts = []
    faces = []
    index = {}

    def vid(p):
        key = (round(p[0], 6), round(p[1], 6), round(p[2], 6))
        i = index.get(key)
        if i is None:
            i = len(verts)
            index[key] = i
            verts.append((p[0] * scale, p[1] * scale, p[2] * scale))
        return i

    for face in body.faces:
        vlist, tris = tessellate.tessellate_face(face)
        for a, b, c in tris:
            ia, ib, ic = vid(vlist[a]), vid(vlist[b]), vid(vlist[c])
            if ia != ib and ib != ic and ia != ic:
                faces.append((ia, ib, ic))

    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.validate(verbose=False)

    # Weld coincident geometry and make normals consistent (outward).
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-6)
    # A solid B-rep tessellates watertight, so any remaining boundary edge is
    # a defect -- typically a sliver where ASM stores the same ring as
    # different edge records on each side (sampled at different densities).
    # Close such small cracks; genuine model holes have walls, not boundaries.
    _weld_t_junctions(bm, max(chord_tol, 1e-5))
    border = [e for e in bm.edges if e.is_boundary]
    if border:
        bmesh.ops.holes_fill(bm, edges=border, sides=8)
        bmesh.ops.triangulate(bm, faces=bm.faces)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    return mesh


def _body_signature(body):
    """A coarse fingerprint of a body's placement (for de-duplication).

    Fusion bakes several coincident bodies into one file (configuration-table
    states / identity-transform instances) that occupy the same space.  They
    may differ slightly in tessellation, so we fingerprint by rounded bounding
    box + face count rather than exact vertices; coincident bodies collapse to
    one and only genuinely separate parts survive.
    """
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    for face in body.faces:
        for ring in face.loops:
            for p in ring:
                for c in range(3):
                    lo[c] = min(lo[c], p[c])
                    hi[c] = max(hi[c], p[c])
    return (
        tuple(round(x, 2) for x in lo),
        tuple(round(x, 2) for x in hi),
        len(body.faces),
    )


def load(context, filepath, deviation=0.1, join_bodies=False,
         global_scale=0.01, merge_duplicates=True):
    """Import ``filepath`` and return the number of objects created."""
    base = os.path.splitext(os.path.basename(filepath))[0]

    with container.open_f3d(filepath) as c:
        blobs = c.brep_blobs()

    if not blobs:
        raise ValueError("no B-rep geometry found in .f3d (unsupported variant?)")

    all_bodies = []
    for blob in blobs:
        parsed = brep.Brep(sab.parse(blob), deviation=deviation)
        all_bodies.extend(parsed.bodies())

    if merge_duplicates:
        unique = []
        seen = set()
        for body in all_bodies:
            sig = _body_signature(body)
            if sig in seen:
                continue
            seen.add(sig)
            unique.append(body)
        all_bodies = unique

    objects = []
    for bi, body in enumerate(all_bodies):
        name = base if len(all_bodies) == 1 else f"{base}.{bi:03d}"
        mesh = _build_mesh(name, body, global_scale,
                           chord_tol=deviation * global_scale)
        if len(mesh.polygons) == 0:
            bpy.data.meshes.remove(mesh)
            continue
        obj = bpy.data.objects.new(name, mesh)
        context.collection.objects.link(obj)
        objects.append(obj)

    if not objects:
        raise ValueError("geometry parsed but nothing could be tessellated")

    # Selection / active object for a nice post-import state.
    for o in context.selected_objects:
        o.select_set(False)
    for o in objects:
        o.select_set(True)
    context.view_layer.objects.active = objects[0]

    if join_bodies and len(objects) > 1:
        ctx = {"active_object": objects[0], "selected_editable_objects": objects}
        with context.temp_override(**ctx):
            bpy.ops.object.join()
        objects = [objects[0]]

    return len(objects)
