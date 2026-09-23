from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

from .analysis.pipeline import run_standardized_analysis
from .analysis.types import AnalysisError, StandardizedAnalysisConfig


def main(argv: Sequence[str] | None = None) -> int:
    args = tuple(sys.argv[1:] if argv is None else argv)
    if args in (("--help",), ("-h",)):
        print("usage: python -m st_virtual_screening.analyze_cli --config PAPER_TOML")
        return 0
    if len(args) != 2 or args[0] != "--config":
        print(
            "usage: python -m st_virtual_screening.analyze_cli --config PAPER_TOML",
            file=sys.stderr,
        )
        return 2
    try:
        config = StandardizedAnalysisConfig.from_toml(Path(args[1]))
        run_standardized_analysis(config)
    except AnalysisError as error:
        print(error, file=sys.stderr)
        return 2
    print(config.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
