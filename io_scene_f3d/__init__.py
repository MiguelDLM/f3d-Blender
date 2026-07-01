"""Blender add-on: import Autodesk Fusion 360 ``.f3d`` archives.

Registers a File > Import entry that reads the ASM B-rep out of the ``.f3d``
ZIP container, reconstructs its topology and tessellates it into a Blender
mesh.  All parsing lives in the pure-Python sibling modules (``container``,
``sab``, ``brep``, ``tessellate``) so they can be unit-tested without Blender.
"""

bl_info = {
    "name": "Import Autodesk Fusion (.f3d)",
    "author": "f3d-blender contributors",
    "version": (0, 1, 0),
    "blender": (4, 2, 0),
    "location": "File > Import > Autodesk Fusion (.f3d)",
    "description": "Import Fusion 360 .f3d archives (ASM B-rep) as meshes",
    "category": "Import-Export",
}

# Blender imports are only available when running inside Blender; guarded so the
# module can still be imported by the standalone test-suite.
try:
    import bpy
    from bpy.props import StringProperty, FloatProperty, BoolProperty
    from bpy.types import Operator
    from bpy_extras.io_utils import ImportHelper
    _HAVE_BPY = True
except ImportError:  # pragma: no cover - exercised only outside Blender
    _HAVE_BPY = False


if _HAVE_BPY:

    class IMPORT_OT_f3d(Operator, ImportHelper):
        """Import an Autodesk Fusion 360 .f3d file."""

        bl_idname = "import_scene.f3d"
        bl_label = "Import Autodesk Fusion (.f3d)"
        bl_options = {"REGISTER", "UNDO"}

        filename_ext = ".f3d"
        filter_glob: StringProperty(default="*.f3d", options={"HIDDEN"})

        deviation: FloatProperty(
            name="Chord deviation",
            description="Max distance between a curved surface and its mesh (in model units). Smaller = finer mesh",
            default=0.1, min=0.001, max=10.0,
        )
        global_scale: FloatProperty(
            name="Scale",
            description="Scale factor. Fusion works in cm; 0.01 maps model cm to Blender m",
            default=0.01, min=1e-6, max=1000.0,
        )
        merge_duplicates: BoolProperty(
            name="Merge duplicates",
            description="Collapse geometrically identical coincident bodies "
                        "(Fusion configuration states / identity instances)",
            default=True,
        )
        join_bodies: BoolProperty(
            name="Join bodies",
            description="Merge all solid bodies into a single mesh object",
            default=False,
        )

        def execute(self, context):
            from . import importer
            try:
                n = importer.load(
                    context, self.filepath,
                    deviation=self.deviation,
                    join_bodies=self.join_bodies,
                    global_scale=self.global_scale,
                    merge_duplicates=self.merge_duplicates,
                )
            except Exception as exc:  # surface a clean error in the UI
                self.report({"ERROR"}, f"F3D import failed: {exc}")
                return {"CANCELLED"}
            self.report({"INFO"}, f"Imported {n} object(s) from F3D")
            return {"FINISHED"}

    def _menu_func_import(self, context):
        self.layout.operator(
            IMPORT_OT_f3d.bl_idname, text="Autodesk Fusion (.f3d)"
        )

    _classes = (IMPORT_OT_f3d,)

    def register():
        for cls in _classes:
            bpy.utils.register_class(cls)
        bpy.types.TOPBAR_MT_file_import.append(_menu_func_import)

    def unregister():
        bpy.types.TOPBAR_MT_file_import.remove(_menu_func_import)
        for cls in reversed(_classes):
            bpy.utils.unregister_class(cls)

else:

    def register():
        raise RuntimeError("io_scene_f3d must be registered from within Blender")

    def unregister():
        pass
