from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

from .analysis.figures import render_figure_pack
from .analysis.types import AnalysisError, FigureConfig


def main(argv: Sequence[str] | None = None) -> int:
    args = tuple(sys.argv[1:] if argv is None else argv)
    if args in (("--help",), ("-h",)):
        print("usage: python -m st_virtual_screening.figures_cli --config PAPER_TOML")
        return 0
    if len(args) != 2 or args[0] != "--config":
        print(
            "usage: python -m st_virtual_screening.figures_cli --config PAPER_TOML",
            file=sys.stderr,
        )
        return 2
    try:
        config = FigureConfig.from_toml(Path(args[1]))
        render_figure_pack(config)
    except AnalysisError as error:
        print(error, file=sys.stderr)
        return 2
    print(config.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
