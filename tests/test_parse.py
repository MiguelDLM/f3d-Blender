"""Standalone sanity tests for the container + SAB parser (no Blender needed).

Run with:  python3 tests/test_parse.py
"""

import os
import sys
import collections

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "io_scene_f3d"))

import container  # noqa: E402
import sab         # noqa: E402

SAMPLE = os.path.join(ROOT, "samples", "soporte de corttinas.f3d")


def main() -> int:
    assert os.path.exists(SAMPLE), f"sample missing: {SAMPLE}"

    with container.open_f3d(SAMPLE) as c:
        print("asset root :", c.asset_root)
        print("brep parts :", c.brep_names)
        print("preview    :", c.preview_name)
        blobs = c.brep_blobs()

    assert blobs, "no BREP parts found"

    total = collections.Counter()
    for blob in blobs:
        f = sab.parse(blob)
        for r in f.records:
            total[r.name] += 1
        print(f"\nheader: {f.header.get('strings', [])[:2]} "
              f"resabs={f.header.get('resabs')}")
        print(f"records: {len(f.records)}")

    print("\nentity counts:")
    for k, v in total.most_common():
        print(f"  {v:5d}  {k}")

    # Basic topological expectations for the known sample.
    for kind in ("body", "shell", "face", "loop", "coedge", "edge", "vertex"):
        assert total[kind] > 0, f"missing {kind} records"
    print("\nOK ✔  topology present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
