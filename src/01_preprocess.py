from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from .preprocessing.errors import PreprocessingError, SampleNotFoundError
from .preprocessing.export import export_cmapwide
from .preprocessing.io import find_samples, load_counts_tsv, load_visium_h5
from .preprocessing.models import ExportConfig, SamplePaths
from .preprocessing.spatial import OrderedFixtureClusterer, ScanpyLeidenClusterer


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m st_virtual_screening.preprocess_cli",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subcommands = parser.add_subparsers(dest="command", required=True)
    audit = subcommands.add_parser(
        "audit-raw",
        help="list Visium RAW samples discoverable in a flat GEO directory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    audit.add_argument("--raw-root", required=True)
    audit.set_defaults(handler=_audit_raw)
    export = subcommands.add_parser(
        "export-wide",
        help="export one 092A Visium sample as a CMap-wide spot table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    export.add_argument("--raw-root", required=True)
    export.add_argument("--sample-id", required=True)
    export.add_argument("--cmap-signatures", required=True)
    export.add_argument("--mapping-dir", required=True)
    export.add_argument("--out-root", required=True)
    export.add_argument("--counts-tsv", default="", help="small spot-by-feature TSV fixture; omit to load the 10x H5")
    export.add_argument("--keep-out-of-tissue", action="store_true")
    export.add_argument("--leiden-resolution", type=float, default=0.6, help="paper Leiden resolution; default: %(default)s")
    export.add_argument("--random-state", type=int, default=0, help="paper random seed; default: %(default)s")
    export.set_defaults(handler=_export_wide)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except PreprocessingError as error:
        print(error, file=sys.stderr)
        return 2


def _audit_raw(args: argparse.Namespace) -> int:
    for sample in find_samples(Path(args.raw_root)):
        print("\t".join((sample.sample_id, str(sample.h5), str(sample.positions), str(sample.scalefactors))))
    return 0


def _export_wide(args: argparse.Namespace) -> int:
    raw_root = Path(args.raw_root)
    sample = _sample_by_id(find_samples(raw_root), args.sample_id, raw_root)
    fixture_counts = str(args.counts_tsv)
    matrix = load_counts_tsv(Path(fixture_counts)) if fixture_counts else load_visium_h5(sample.h5)
    clusterer = OrderedFixtureClusterer() if fixture_counts else ScanpyLeidenClusterer()
    result = export_cmapwide(ExportConfig(
        sample=sample,
        matrix=matrix,
        cmap_signature_path=Path(args.cmap_signatures),
        mapping_dir=Path(args.mapping_dir),
        out_root=Path(args.out_root),
        clusterer=clusterer,
        in_tissue_only=not bool(args.keep_out_of_tissue),
        leiden_resolution=float(args.leiden_resolution),
        random_state=int(args.random_state),
    ))
    print(result.wide_tsv)
    return 0


def _sample_by_id(samples: list[SamplePaths], sample_id: str, raw_root: Path) -> SamplePaths:
    for sample in samples:
        if sample.sample_id == sample_id:
            return sample
    raise SampleNotFoundError(sample_id=sample_id, raw_root=raw_root)


if __name__ == "__main__":
    raise SystemExit(main())
