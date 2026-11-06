#!/usr/bin/env python3
"""Run Porechop ABI and Chopper once per expected demultiplexed sample."""

import argparse
import csv
import gzip
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_pipeline import load_config, read_manifest
from scripts.prepare_inputs import demux_filename


def preprocess_reads(
    manifest_path: Path,
    config_dir: Path,
    demux_dir: Path,
    output_dir: Path,
) -> None:
    config, _ = load_config(config_dir)
    samples = read_manifest(manifest_path, config)
    demux_dir = demux_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")

    preprocess = config.get("preprocess", {})
    porechop_options = shlex.split(preprocess.get("porechop_abi_options", ""))
    chopper_options = shlex.split(preprocess.get("chopper_options", ""))
    threads = str(config["threads"])
    included = []
    excluded = []
    method = config["demux"]["method"]
    for sample in samples:
        filename = demux_filename(method, sample["barcode_indices"])
        reads = demux_dir / filename
        if not reads.exists():
            excluded.append((sample["sample_id"], "no_demux_file"))
            continue
        with gzip.open(reads, "rb") as handle:
            has_reads = bool(handle.read(1))
        if not has_reads:
            excluded.append((sample["sample_id"], "empty_fastq"))
            continue
        included.append((filename, reads))

    output_dir.mkdir(parents=True)
    with (output_dir / "excluded_samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "reason"])
        writer.writerows(excluded)

    for filename, reads in included:
        output = output_dir / filename
        with tempfile.TemporaryDirectory() as temporary:
            porechop_output = Path(temporary) / "porechop.fastq"
            subprocess.run(
                [
                    "porechop_abi", "-abi", *porechop_options,
                    "-t", threads, "-i", str(reads), "-o", str(porechop_output),
                ],
                check=True,
            )
            with gzip.open(output, "wb") as handle:
                subprocess.run(
                    [
                        "chopper", *chopper_options, "-t", threads,
                        "-i", str(porechop_output),
                    ],
                    stdout=handle,
                    check=True,
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--demux-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    preprocess_reads(
        args.manifest, args.config_dir, args.demux_dir, args.output_dir,
    )


if __name__ == "__main__":
    main()
