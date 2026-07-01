"""End-to-end geometry test (no Blender): container -> sab -> brep -> tessellate.

Run with:  python3 tests/test_pipeline.py
"""

import os
import sys

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "io_scene_f3d"))

import container  # noqa: E402
import sab         # noqa: E402
import brep        # noqa: E402
import tessellate  # noqa: E402

SAMPLE = os.path.join(ROOT, "samples", "soporte de corttinas.f3d")


def main() -> int:
    with container.open_f3d(SAMPLE) as c:
        blobs = c.brep_blobs()

    total_bodies = 0
    total_faces = 0
    total_tris = 0
    failed = 0
    for blob in blobs:
        bodies = brep.Brep(sab.parse(blob)).bodies()
        total_bodies += len(bodies)
        for body in bodies:
            for face in body.faces:
                total_faces += 1
                verts, tris = tessellate.tessellate_face(face)
                total_tris += len(tris)
                if not tris:
                    failed += 1

    print(f"bodies={total_bodies} faces={total_faces} "
          f"triangles={total_tris} failed={failed}")

    assert total_bodies == 8, f"expected 8 bodies, got {total_bodies}"
    assert total_faces > 100, "too few faces"
    assert total_tris > 300, "too few triangles"
    assert failed <= 2, f"too many faces failed to tessellate: {failed}"
    print("OK ✔  full pipeline produced a mesh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
