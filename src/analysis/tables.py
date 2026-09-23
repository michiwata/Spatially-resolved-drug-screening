from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from pathlib import Path

from .types import MissingColumnError, OutputRow, Table


def read_tsv(path: Path) -> Table:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        rows = [{key: value for key, value in row.items()} for row in reader]
    return Table(path=path, header=header, rows=rows)


def require_columns(table: Table, columns: Iterable[str]) -> None:
    missing = tuple(column for column in columns if column not in table.header)
    if missing:
        raise MissingColumnError(path=table.path, columns=missing)


def hsa_columns(table: Table) -> list[str]:
    return [column for column in table.header if column.startswith("hsa:")]


def write_tsv(path: Path, rows: Sequence[OutputRow], header: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(header), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def read_required(path: Path, columns: Iterable[str]) -> Table:
    table = read_tsv(path)
    require_columns(table, columns)
    return table


def sort_numeric_text(values: Iterable[str]) -> list[str]:
    return sorted(values, key=lambda value: int(value) if value.isdigit() else value)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(text, encoding="utf-8")
