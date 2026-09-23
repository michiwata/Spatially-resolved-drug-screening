from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class PreprocessingError(Exception):
    pass


@dataclass(frozen=True, slots=True)
class InputSchemaError(PreprocessingError):
    path: Path
    reason: str

    def __str__(self) -> str:
        return f"{self.path}: {self.reason}"


@dataclass(frozen=True, slots=True)
class BarcodeMismatchError(PreprocessingError):
    sample_id: str
    reason: str

    def __str__(self) -> str:
        return f"{self.sample_id}: {self.reason}"


@dataclass(frozen=True, slots=True)
class GeneMappingError(PreprocessingError):
    reason: str

    def __str__(self) -> str:
        return self.reason


@dataclass(frozen=True, slots=True)
class SampleNotFoundError(PreprocessingError):
    sample_id: str
    raw_root: Path

    def __str__(self) -> str:
        return f"sample {self.sample_id!r} was not found under {self.raw_root}"


@dataclass(frozen=True, slots=True)
class OptionalDependencyError(PreprocessingError):
    package: str
    purpose: str

    def __str__(self) -> str:
        return f"install {self.package} to {self.purpose}"
