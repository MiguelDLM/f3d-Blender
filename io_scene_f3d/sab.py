"""Tokenizer / record parser for the ASM ``BinaryFile`` (ACIS-SAB style) stream.

The geometry inside a Fusion ``.f3d`` lives in a ``BREP.*.smb`` entry whose
payload is Autodesk Shape Manager (ASM) binary.  It uses the same tokenized
record layout as Spatial's "Standard ACIS Binary" (``.sab``): a stream of
tagged values grouped into records, each record terminated by an
end-of-record tag.  Records reference each other by their ordinal index
(0-based, in stream order); ``-1`` means a null reference.

This module turns the raw bytes into a list of :class:`Record` objects.  It is
deliberately schema-agnostic: it preserves the raw ordered list of values per
record so higher layers (:mod:`brep`) can interpret them per entity type.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from typing import Any

MAGIC = b"ASM BinaryFile"

# --- token tags -----------------------------------------------------------
TAG_INT = 0x04       # signed int32
TAG_FLOAT = 0x05     # 32-bit float (rare)
TAG_DOUBLE = 0x06    # 64-bit IEEE double
TAG_STRING = 0x07    # u8-length-prefixed byte string
TAG_STRING2 = 0x08   # u8-length-prefixed byte string (alt use)
TAG_ENUM_A = 0x0A    # 0-byte enum / logical
TAG_ENUM_B = 0x0B    # 0-byte enum / logical
TAG_POINTER = 0x0C   # signed int32 record index (-1 == null)
TAG_ENTITY = 0x0D    # u8-length entity type name; starts a record body
TAG_SUBTYPE_OPEN = 0x0E   # u8-length subtype name; opens a nested subtype
TAG_SUBTYPE_CLOSE = 0x0F  # closes a nested subtype
TAG_MARKER = 0x10    # 0-byte marker (optional-extension flag)
TAG_EOR = 0x11       # end of record
TAG_POSITION = 0x13  # a 3D position: three packed doubles (x, y, z)
TAG_VECTOR = 0x14    # a 3D vector/direction: three packed doubles
TAG_U32 = 0x15       # unsigned int32 (flags / counts)


class SabError(Exception):
    """Raised when the SAB stream cannot be tokenized."""


@dataclass
class Ref:
    """A pointer to another record, by ordinal index (-1 == null)."""

    index: int

    @property
    def is_null(self) -> bool:
        return self.index < 0


@dataclass
class TypeName:
    """A polymorphic type-name token appearing inside a record body.

    ACIS/ASM stores a derived entity by first writing its base class name and
    then the derived name(s); a single record can therefore contain several
    of these (e.g. ``curve`` -> ``exact_int_cur`` -> ``nubs``).  Only the
    *first* one is used as the record's :attr:`Record.name`.
    """

    name: str


@dataclass
class Record:
    """One ASM entity record: an ordered list of tokens between two EORs.

    ``values`` entries are plain Python objects: ``int``, ``float``, ``str``,
    ``tuple`` (a 3D position/vector), :class:`Ref`, :class:`TypeName`, or
    ``None`` (for 0-byte enum/marker tags).  ``name`` is the first
    :class:`TypeName` seen (the entity's most-derived leaf type comes last).
    """

    index: int
    name: str
    values: list = field(default_factory=list)

    def refs(self) -> list:
        return [v for v in self.values if isinstance(v, Ref)]

    def doubles(self) -> list:
        return [v for v in self.values if isinstance(v, float)]

    def type_names(self) -> list:
        return [v.name for v in self.values if isinstance(v, TypeName)]

    @property
    def leaf_type(self) -> str:
        """The most-derived (last) type name, i.e. the concrete geometry."""
        names = self.type_names()
        return names[-1] if names else self.name

    def positions(self) -> list:
        return [v for v in self.values if isinstance(v, tuple)]


@dataclass
class SabFile:
    header: dict
    records: list
    _ordinals: list = field(default=None, repr=False)

    def by_type(self, name: str) -> list:
        return [r for r in self.records if r.name == name]

    def resolve(self, ref) -> "Record | None":
        """Resolve a pointer to its record.

        Pointer ordinals do NOT count every physical record: history journal
        files (``.smbh``) interleave ``delta_state`` and ``Begin`` markers
        that the writer skips when numbering entities (validated: with them
        excluded, 100% of face->surface pointers resolve to surfaces).  Files
        without such markers are unaffected (the mapping is the identity).
        """
        idx = ref.index if isinstance(ref, Ref) else ref
        if idx is None or idx < 0:
            return None
        if self._ordinals is None:
            self._ordinals = [
                r for r in self.records
                if r.name not in ("delta_state", "Begin")
            ]
        if idx >= len(self._ordinals):
            return None
        return self._ordinals[idx]


class _Reader:
    __slots__ = ("b", "i", "n")

    def __init__(self, b: bytes, i: int = 0):
        self.b = b
        self.i = i
        self.n = len(b)

    def u8(self) -> int:
        v = self.b[self.i]
        self.i += 1
        return v

    def i32(self) -> int:
        v = struct.unpack_from("<i", self.b, self.i)[0]
        self.i += 4
        return v

    def u32(self) -> int:
        v = struct.unpack_from("<I", self.b, self.i)[0]
        self.i += 4
        return v

    def f32(self) -> float:
        v = struct.unpack_from("<f", self.b, self.i)[0]
        self.i += 4
        return v

    def f64(self) -> float:
        v = struct.unpack_from("<d", self.b, self.i)[0]
        self.i += 8
        return v

    def string(self) -> str:
        ln = self.u8()
        s = self.b[self.i:self.i + ln]
        self.i += ln
        return s.decode("latin1")

    def vec3(self) -> tuple:
        v = struct.unpack_from("<3d", self.b, self.i)
        self.i += 24
        return v


SUBTYPE_OPEN = object()   # sentinel appended for a 0x0E subtype-open marker
SUBTYPE_CLOSE = object()  # sentinel appended for a 0x0F subtype-close marker


def _read_value(r: _Reader, tag: int):
    """Read one non-structural value (returns object) or a control sentinel."""
    if tag == TAG_INT:
        return r.i32()
    if tag == TAG_POINTER:
        return Ref(r.i32())
    if tag == TAG_DOUBLE:
        return r.f64()
    if tag == TAG_FLOAT:
        return r.f32()
    if tag == TAG_U32:
        return r.u32()
    if tag in (TAG_STRING, TAG_STRING2):
        return r.string()
    if tag in (TAG_POSITION, TAG_VECTOR):
        return r.vec3()
    if tag == TAG_ENUM_A:
        return False   # logical / enum, first value (e.g. forward)
    if tag == TAG_ENUM_B:
        return True    # logical / enum, second value (e.g. reversed)
    if tag == TAG_MARKER:
        return None
    raise SabError(f"unhandled value tag 0x{tag:02x} at offset {r.i - 1}")


def parse(data: bytes) -> SabFile:
    """Tokenize an ASM ``BinaryFile`` blob into a :class:`SabFile`.

    Records are delimited *only* by the end-of-record tag (0x11).  Entity type
    names (0x0D) may appear several times within one record (polymorphic
    base/derived chain), so they are stored as :class:`TypeName` tokens rather
    than used as record delimiters.  Subtype open/close (0x0E/0x0F) are kept as
    flat structural markers, which keeps the byte stream perfectly in sync
    without having to model ASM's nested save/restore semantics.
    """
    if not data.startswith(MAGIC):
        raise SabError("not an ASM BinaryFile (bad magic)")

    r = _Reader(data)
    r.i = 15   # past "ASM BinaryFile<n>"
    r.i += 16  # fixed preamble: two 8-byte / int32 fields

    records: list = []
    cur_vals: list = []
    cur_name: str | None = None

    def flush():
        nonlocal cur_vals, cur_name
        if cur_vals or cur_name is not None:
            records.append(Record(len(records), cur_name or "", cur_vals))
        cur_vals = []
        cur_name = None

    while r.i < r.n:
        tag = r.u8()
        if tag == TAG_EOR:
            flush()
        elif tag == TAG_ENTITY:
            name = r.string()
            if cur_name is None:
                cur_name = name
            cur_vals.append(TypeName(name))
        elif tag == TAG_SUBTYPE_OPEN:
            name = r.string()
            if cur_name is None:
                cur_name = name
            cur_vals.append(TypeName(name))
            cur_vals.append(SUBTYPE_OPEN)
        elif tag == TAG_SUBTYPE_CLOSE:
            cur_vals.append(SUBTYPE_CLOSE)
        else:
            cur_vals.append(_read_value(r, tag))
    flush()

    header: dict = {}
    if records and records[0].name == "asmheader":
        h = records[0]
        strs = [v for v in h.values if isinstance(v, str)]
        dbls = [v for v in h.values if isinstance(v, float)]
        header = {"strings": strs, "doubles": dbls}
        if len(dbls) >= 3:
            header["resabs"] = dbls[1]  # absolute tolerance (observed 1e-6)

    return SabFile(header=header, records=records)


def parse_file(path: str) -> SabFile:
    with open(path, "rb") as fh:
        return parse(fh.read())
