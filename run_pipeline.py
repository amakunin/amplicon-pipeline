#!/usr/bin/env python3
"""Entrypoint scaffold for the fixed-order amplicon pipeline."""

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "ONT amplicon pipeline scaffold. Analysis execution is not "
            "implemented; no outputs will be created."
        )
    )
    parser.add_argument("run_fastq_gz", type=Path, help="Raw run FASTQ (.fastq.gz)")
    parser.add_argument("manifest_tsv", type=Path, help="Sample manifest (.tsv)")
    parser.add_argument("config_dir", type=Path, help="Run configuration directory")
    parser.add_argument("output_dir", type=Path, help="Analysis output directory")
    parser.add_argument(
        "--profile",
        choices=("local", "lsf"),
        default="local",
        help="Execution environment (default: local; not yet implemented)",
    )
    parser.parse_args()
    parser.exit(
        status=1,
        message=(
            "Pipeline execution is not implemented. This scaffold does not "
            "read inputs, launch jobs or create outputs.\n"
        ),
    )


if __name__ == "__main__":
    main()
