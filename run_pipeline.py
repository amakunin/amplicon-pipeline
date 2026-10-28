#!/usr/bin/env python3
"""Entrypoint scaffold for the fixed-order amplicon pipeline."""

import argparse
from pathlib import Path
from typing import Any

import yaml


def load_config(config_dir: Path) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    config_dir = config_dir.resolve()
    with (config_dir / "pipeline.yaml").open() as handle:
        config = yaml.safe_load(handle)

    demux = config["demux"]
    if demux["method"] not in ("hybrid", "illumina"):
        raise ValueError("demux.method must be 'hybrid' or 'illumina'")
    for key in (
        "reverse_barcode_fasta",
        "reverse_barcode_arrangement",
        "ont_barcode_mapping",
        "front_barcode_mapping",
        "reverse_barcode_mapping",
    ):
        if key in demux:
            demux[key] = (config_dir / demux[key]).resolve()

    targets = {}
    for name, filename in config["targets"].items():
        path = (config_dir / filename).resolve()
        config["targets"][name] = path
        with path.open() as handle:
            parameters = yaml.safe_load(handle)
        if "sintax_ref" in parameters:
            parameters["sintax_ref"] = str(
                (config_dir / parameters["sintax_ref"]).resolve()
            )
        targets[name] = parameters
    return config, targets


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
    args = parser.parse_args()
    load_config(args.config_dir)
    parser.exit(
        status=1,
        message=(
            "Pipeline execution is not implemented. Configuration was loaded; "
            "no jobs were launched or outputs created.\n"
        ),
    )


if __name__ == "__main__":
    main()
