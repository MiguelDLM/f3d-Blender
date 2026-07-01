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

try:
    from . import container, sab, brep, tessellate
except ImportError:  # pragma: no cover - allow flat import
    import container, sab, brep, tessellate


def _build_mesh(name, body, scale):
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
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    return mesh


def _body_signature(body):
    """A hashable fingerprint of a body's geometry (for de-duplication).

    Fusion may bake several geometrically identical, coincident bodies into one
    file (configuration states / identity-transform instances).  Bodies with the
    same rounded vertex set collapse to one signature.
    """
    coords = set()
    for face in body.faces:
        for ring in face.loops:
            for p in ring:
                coords.add((round(p[0], 5), round(p[1], 5), round(p[2], 5)))
    return frozenset(coords)


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
        all_bodies.extend(brep.Brep(sab.parse(blob)).bodies())

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
        mesh = _build_mesh(name, body, global_scale)
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
