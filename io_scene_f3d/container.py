"""Read the ``.f3d`` ZIP container and locate its parts.

A Fusion ``.f3d`` is a ZIP archive.  This module opens it and exposes the
interesting members: the ASM B-rep blob(s), the preview thumbnail and the raw
design/metadata streams.  It performs no geometry parsing (see :mod:`sab`).
"""

from __future__ import annotations

import zipfile
from dataclasses import dataclass, field


@dataclass
class F3DContainer:
    path: str
    asset_root: str = ""                    # e.g. "FusionAssetName[Active]"
    brep_names: list = field(default_factory=list)   # names of BREP.*.smb
    preview_name: str | None = None
    _zip: zipfile.ZipFile | None = None

    def read(self, name: str) -> bytes:
        assert self._zip is not None
        return self._zip.read(name)

    def brep_blobs(self) -> list:
        """Return the raw bytes of every ASM B-rep part in the archive."""
        return [self.read(n) for n in self.brep_names]

    def close(self) -> None:
        if self._zip is not None:
            self._zip.close()
            self._zip = None

    def __enter__(self) -> "F3DContainer":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def open_f3d(path: str) -> F3DContainer:
    """Open a ``.f3d`` archive and index its notable members."""
    zf = zipfile.ZipFile(path, "r")
    names = zf.namelist()

    # Newer files ship two ASM parts: BREP.*.smbh holds the CURRENT model
    # (bodies at their final world positions), while BREP.*.smb holds the
    # per-feature design history (dozens of intermediate body snapshots).
    # Prefer the model; fall back to the history stream for older files
    # that only carry an .smb (whose last snapshots are the final bodies).
    model = sorted(
        n for n in names
        if n.rsplit("/", 1)[-1].upper().startswith("BREP.")
        and n.lower().endswith(".smbh")
    )
    breps = model or sorted(
        n for n in names
        if n.rsplit("/", 1)[-1].upper().startswith("BREP.")
        and n.lower().endswith(".smb")
    )
    preview = next(
        (n for n in names if n.lower().endswith("previews/small.png")), None
    )
    asset_root = ""
    if breps:
        # ".../Breps.BlobParts/BREP.xxx.smb" -> asset root is two levels up
        parts = breps[0].split("/")
        if len(parts) >= 3:
            asset_root = parts[0]

    return F3DContainer(
        path=path,
        asset_root=asset_root,
        brep_names=breps,
        preview_name=preview,
        _zip=zf,
    )


def is_f3d(path: str) -> bool:
    """Cheap check: a valid ZIP whose first entry name starts with ``Fusion``."""
    try:
        with zipfile.ZipFile(path, "r") as zf:
            names = zf.namelist()
    except (zipfile.BadZipFile, OSError):
        return False
    return any(
        n.startswith("Fusion") or n.startswith("Manifest.dat") for n in names
    )
