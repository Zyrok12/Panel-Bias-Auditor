from __future__ import annotations

from dataclasses import dataclass, field


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome labels for interval comparison."""
    value = chrom.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper()


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    name: str = "."
    source: str = "unknown"
    attrs: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"Interval start must be >= 0: {self}")
        if self.end <= self.start:
            raise ValueError(f"Interval end must be greater than start: {self}")

    @property
    def norm_chrom(self) -> str:
        return normalize_chrom(self.chrom)

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlaps(self, other: "Interval") -> bool:
        return self.norm_chrom == other.norm_chrom and self.start < other.end and other.start < self.end

    def overlap_span(self, other: "Interval") -> tuple[int, int] | None:
        if not self.overlaps(other):
            return None
        return max(self.start, other.start), min(self.end, other.end)


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    identifier: str = "."
    qual: str = "."
    filt: str = "."
    info: dict[str, str] = field(default_factory=dict)

    @property
    def norm_chrom(self) -> str:
        return normalize_chrom(self.chrom)

    @property
    def start(self) -> int:
        return self.pos - 1

    @property
    def end(self) -> int:
        if "END" in self.info:
            try:
                return int(self.info["END"])
            except ValueError:
                pass
        return self.start + max(len(self.ref), 1)

    @property
    def interval(self) -> Interval:
        label = self.identifier if self.identifier != "." else f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"
        return Interval(self.chrom, self.start, self.end, label, source="vcf")


@dataclass
class RegionCoverage:
    region: Interval
    covered_bases: int

    @property
    def missing_bases(self) -> int:
        return self.region.length - self.covered_bases

    @property
    def coverage_fraction(self) -> float:
        return self.covered_bases / self.region.length if self.region.length else 0.0

    @property
    def status(self) -> str:
        if self.covered_bases == 0:
            return "missing"
        if self.covered_bases == self.region.length:
            return "covered"
        return "partial"

