#!/usr/bin/env python3
"""Entrypoint scaffold for the fixed-order amplicon pipeline."""

import argparse
import csv
import re
import shlex
import subprocess
import sys
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
        "barcode_dir",
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
        for key in ("sintax_ref", "sintax_ref_tax_custom"):
            if key in parameters:
                parameters[key] = str((config_dir / parameters[key]).resolve())
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


def demux_command(
    config: dict[str, Any], run_fastq: Path, output_dir: Path
) -> list[str]:
    demux = config["demux"]
    command = [
        "bash",
        str(Path(__file__).resolve().parent / f"scripts/demux_{demux['method']}.sh"),
        "--input", str(run_fastq.resolve()),
        "--output-dir", str(output_dir.resolve() / "demux"),
        "--barcode-dir", str(demux["barcode_dir"]),
        "--threads", str(config["threads"]),
    ]
    if demux.get("dorado_exec"):
        command += ["--dorado", str(demux["dorado_exec"])]
    return command


def check_nextflow() -> None:
    try:
        output = subprocess.run(
            ["nextflow", "-version"], capture_output=True, text=True, check=True
        ).stdout
    except FileNotFoundError:
        raise SystemExit("nextflow not found on PATH (version 26 or higher required)")
    version = re.search(r"version (\d+)", output)
    if version is None or int(version.group(1)) < 26:
        raise SystemExit(f"nextflow 26 or higher is required by ampliseq 3.0.0:{output}")


def write_resource_config(target_dir: Path, threads: int) -> Path:
    path = target_dir / "resources.config"
    path.write_text(
        f"executor.cpus = {threads}\n"
        f"process.resourceLimits = [cpus: {threads}]\n"
    )
    return path


def ampliseq_command(
    config_dir: Path, target_dir: Path, profile: str
) -> list[str]:
    return [
        "nextflow", "run", "nf-core/ampliseq", "-r", "3.0.0",
        "-profile", profile,
        "-c", str(config_dir.resolve() / "ampliseq.config"),
        "-c", str(target_dir / "resources.config"),
        "-params-file", str(target_dir / "params.yaml"),
    ]


def run_ampliseq(
    config: dict[str, Any],
    args: argparse.Namespace,
) -> None:
    ampliseq_dir = args.output_dir.resolve() / "ampliseq"
    for target in config["targets"]:
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent / "scripts/prepare_inputs.py"),
                "--manifest", str(args.manifest_tsv),
                "--config-dir", str(args.config_dir),
                "--demux-dir", str(args.output_dir.resolve() / "demux"),
                "--target", target,
                "--output-dir", str(ampliseq_dir),
            ],
            check=True,
        )
        target_dir = ampliseq_dir / target
        write_resource_config(target_dir, config["threads"])
        subprocess.run(
            ampliseq_command(args.config_dir, target_dir, args.profile),
            cwd=target_dir,
            check=True,
        )


def print_plan(config: dict[str, Any], samples: list[dict[str, Any]]) -> None:
    print(f"Method: {config['demux']['method']}")
    print(f"Targets: {', '.join(config['targets'])}")
    print("Samples:")
    for sample in samples:
        indices = ", ".join(
            f"{name}={index}"
            for name, index in sample["barcode_indices"].items()
        )
        print(f"  {sample['sample_id']}: {indices}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "ONT amplicon pipeline scaffold. Only demultiplexing and ampliseq "
            "run; later stages are not implemented."
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
        help="Nextflow profile from ampliseq.config (default: local)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the loaded plan without creating outputs or running analysis",
    )
    args = parser.parse_args()
    config, _ = load_config(args.config_dir)
    samples = read_manifest(args.manifest_tsv, config)
    command = demux_command(config, args.run_fastq_gz, args.output_dir)
    if args.dry_run:
        print_plan(config, samples)
        print(f"Demux command: {shlex.join(command)}")
        for target in config["targets"]:
            target_dir = args.output_dir.resolve() / "ampliseq" / target
            ampliseq = ampliseq_command(args.config_dir, target_dir, args.profile)
            print(f"Ampliseq command ({target}): {shlex.join(ampliseq)}")
        return
    check_nextflow()
    subprocess.run(command, check=True)
    run_ampliseq(config, args)
    parser.exit(
        status=1,
        message="Stages after ampliseq are not implemented.\n",
    )


if __name__ == "__main__":
    main()
