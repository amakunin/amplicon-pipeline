#!/usr/bin/env python3
"""Entrypoint scaffold for the fixed-order amplicon pipeline."""

import argparse
import csv
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
        "dorado_exec",
        "seqtk_exec",
        "reverse_barcode_fasta",
        "reverse_barcode_arrangement",
        "ont_barcode_mapping",
        "front_barcode_mapping",
        "reverse_barcode_mapping",
    ):
        if key in demux and demux[key] is not None:
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


def load_barcode_mapping(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {
            row["tag_sequence"].casefold(): row["index"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_manifest(
    manifest_path: Path, config: dict[str, Any]
) -> list[dict[str, Any]]:
    method = config["demux"]["method"]
    barcode_mappings = (
        {
            "ont_barcode": load_barcode_mapping(
                config["demux"]["ont_barcode_mapping"]
            ),
            "rear_barcode": load_barcode_mapping(
                config["demux"]["reverse_barcode_mapping"]
            ),
        }
        if method == "hybrid"
        else {
            "front_barcode": load_barcode_mapping(
                config["demux"]["front_barcode_mapping"]
            ),
            "rear_barcode": load_barcode_mapping(
                config["demux"]["reverse_barcode_mapping"]
            ),
        }
    )
    required_columns = ("sample_id", *barcode_mappings)
    samples = []
    sample_ids = set()
    barcode_pairs = set()

    with manifest_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for column in required_columns:
                row[column]
            sample_id = row["sample_id"]
            if sample_id in sample_ids:
                raise ValueError(f"Duplicate sample_id: {sample_id}")
            sample_ids.add(sample_id)

            indices = {}
            for column, mapping in barcode_mappings.items():
                barcode = row[column]
                try:
                    indices[column] = mapping[barcode.casefold()]
                except KeyError:
                    raise ValueError(
                        f"Unknown {column} {barcode!r} for sample {sample_id}"
                    ) from None

            barcode_pair = tuple(indices.values())
            if barcode_pair in barcode_pairs:
                raise ValueError(
                    f"Duplicate barcode pair for sample {sample_id}: {barcode_pair}"
                )
            barcode_pairs.add(barcode_pair)
            row["barcode_indices"] = indices
            samples.append(row)

    return samples


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
    config, _ = load_config(args.config_dir)
    read_manifest(args.manifest_tsv, config)
    parser.exit(
        status=1,
        message=(
            "Pipeline execution is not implemented. Configuration was loaded; "
            "no jobs were launched or outputs created.\n"
        ),
    )


if __name__ == "__main__":
    main()
