from __future__ import annotations

from pathlib import Path

from .errors import InputFormatError
from .models import Interval, Variant


def parse_attrs(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for token in raw.replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            attrs[key.strip()] = value.strip()
        else:
            attrs[token] = "true"
    return attrs


def parse_bed(path: str | Path, source: str) -> list[Interval]:
    intervals: list[Interval] = []
    path = Path(path)
    if not path.exists():
        raise InputFormatError(f"BED file does not exist: {path}")
    if not path.is_file():
        raise InputFormatError(f"BED path is not a file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise InputFormatError(f"{path}:{line_number}: BED line requires at least 3 columns: chrom start end")
            chrom = fields[0]
            if not chrom:
                raise InputFormatError(f"{path}:{line_number}: BED chromosome cannot be empty")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise InputFormatError(f"{path}:{line_number}: BED start/end must be integers") from exc
            if start < 0:
                raise InputFormatError(f"{path}:{line_number}: BED start must be >= 0")
            if end <= start:
                raise InputFormatError(f"{path}:{line_number}: BED end must be greater than start")
            name = fields[3] if len(fields) >= 4 else f"{chrom}:{start}-{end}"
            attrs = parse_attrs(fields[4]) if len(fields) >= 5 else {}
            try:
                intervals.append(Interval(chrom, start, end, name=name, source=source, attrs=attrs))
            except ValueError as exc:
                raise InputFormatError(f"{path}:{line_number}: invalid BED interval: {exc}") from exc
    if not intervals:
        raise InputFormatError(f"BED file contains no intervals after comments/blank lines: {path}")
    return intervals


def parse_vcf_info(info: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    if info in ("", "."):
        return parsed
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = "true"
    return parsed


def parse_vcf(path: str | Path) -> list[Variant]:
    variants: list[Variant] = []
    path = Path(path)
    if not path.exists():
        raise InputFormatError(f"VCF file does not exist: {path}")
    if not path.is_file():
        raise InputFormatError(f"VCF path is not a file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise InputFormatError(f"{path}:{line_number}: VCF line requires at least 8 tab-delimited columns")
            chrom, pos_raw, identifier, ref, alt, qual, filt, info_raw = fields[:8]
            if not chrom:
                raise InputFormatError(f"{path}:{line_number}: VCF CHROM cannot be empty")
            try:
                pos = int(pos_raw)
            except ValueError as exc:
                raise InputFormatError(f"{path}:{line_number}: VCF POS must be an integer") from exc
            if pos <= 0:
                raise InputFormatError(f"{path}:{line_number}: VCF POS must be >= 1")
            if not ref or ref == ".":
                raise InputFormatError(f"{path}:{line_number}: VCF REF must be present")
            if not alt or alt == ".":
                raise InputFormatError(f"{path}:{line_number}: VCF ALT must be present")
            variants.append(
                Variant(
                    chrom=chrom,
                    pos=pos,
                    identifier=identifier,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    filt=filt,
                    info=parse_vcf_info(info_raw),
                )
            )
    if not variants:
        raise InputFormatError(f"VCF file contains no variants after headers/comments: {path}")
    return variants
