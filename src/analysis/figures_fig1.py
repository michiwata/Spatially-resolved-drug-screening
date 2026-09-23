from __future__ import annotations

from pathlib import Path

from .types import FigureConfig


def copy_fig1(config: FigureConfig, outdir: Path) -> None:
    static = config.figure1_static
    if static.exists() and static.is_file():
        dest = outdir / "fig1_overview.png"
        dest.write_bytes(static.read_bytes())
        return
    raise FileNotFoundError(f"Figure 1 static asset not found: {static}")
