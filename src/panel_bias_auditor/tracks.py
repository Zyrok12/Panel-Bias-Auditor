from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import InputFormatError
from .models import Interval
from .parsers import parse_bed


VALID_TRACK_ROLES = {"critical", "difficult"}


@dataclass(frozen=True)
class TrackDefinition:
    name: str
    role: str
    path: Path
    description: str = ""
    source: str = ""
    version: str = ""
    license: str = ""

    @classmethod
    def from_dict(cls, raw: dict[str, Any], base_dir: Path) -> "TrackDefinition":
        name = str(raw.get("name", "")).strip()
        role = str(raw.get("role", "")).strip().lower()
        path_raw = str(raw.get("path", "")).strip()
        if not name:
            raise InputFormatError("Track manifest entry is missing required field: name")
        if role not in VALID_TRACK_ROLES:
            raise InputFormatError(f"Track '{name}' has invalid role '{role}'. Expected one of {sorted(VALID_TRACK_ROLES)}")
        if not path_raw:
            raise InputFormatError(f"Track '{name}' is missing required field: path")
        path = Path(path_raw)
        if not path.is_absolute():
            path = base_dir / path
        if not path.exists():
            raise InputFormatError(f"Track '{name}' points to a missing BED file: {path}")
        return cls(
            name=name,
            role=role,
            path=path,
            description=str(raw.get("description", "")),
            source=str(raw.get("source", "")),
            version=str(raw.get("version", "")),
            license=str(raw.get("license", "")),
        )

    def to_metadata(self) -> dict[str, str]:
        return {
            "name": self.name,
            "role": self.role,
            "path": str(self.path),
            "description": self.description,
            "source": self.source,
            "version": self.version,
            "license": self.license,
        }


@dataclass(frozen=True)
class TrackBundle:
    name: str
    version: str
    genome_build: str
    description: str
    tracks: tuple[TrackDefinition, ...]

    @property
    def metadata(self) -> dict[str, object]:
        return {
            "name": self.name,
            "version": self.version,
            "genome_build": self.genome_build,
            "description": self.description,
            "tracks": [track.to_metadata() for track in self.tracks],
        }


def load_track_manifest(path: str | Path) -> TrackBundle:
    manifest_path = Path(path)
    if not manifest_path.exists():
        raise InputFormatError(f"Track manifest does not exist: {manifest_path}")
    if not manifest_path.is_file():
        raise InputFormatError(f"Track manifest path is not a file: {manifest_path}")
    try:
        raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise InputFormatError(f"{manifest_path}: invalid JSON track manifest: {exc}") from exc
    base_dir = manifest_path.parent
    tracks_raw = raw.get("tracks", [])
    if not isinstance(tracks_raw, list) or not tracks_raw:
        raise InputFormatError("Track manifest must contain a non-empty 'tracks' list")

    tracks = tuple(TrackDefinition.from_dict(item, base_dir) for item in tracks_raw)
    roles = {track.role for track in tracks}
    if "critical" not in roles:
        raise InputFormatError("Track manifest must include at least one critical track")

    return TrackBundle(
        name=str(raw.get("name", manifest_path.stem)),
        version=str(raw.get("version", "")),
        genome_build=str(raw.get("genome_build", "unknown")),
        description=str(raw.get("description", "")),
        tracks=tracks,
    )


def load_tracks_from_bundle(bundle: TrackBundle) -> tuple[list[Interval], list[Interval]]:
    critical: list[Interval] = []
    difficult: list[Interval] = []
    for track in bundle.tracks:
        intervals = parse_bed(track.path, source=f"{track.role}:{track.name}")
        if track.role == "critical":
            critical.extend(intervals)
        elif track.role == "difficult":
            difficult.extend(intervals)
    return critical, difficult
