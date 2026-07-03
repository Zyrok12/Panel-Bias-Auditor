from __future__ import annotations


class AuditorError(Exception):
    """Base class for user-facing auditor errors."""


class InputFormatError(AuditorError):
    """Raised when an input file is malformed or cannot be used."""

