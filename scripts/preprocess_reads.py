#!/usr/bin/env python3
"""Run Porechop ABI and Chopper once per expected demultiplexed sample."""

import argparse
import csv
import gzip
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_pipeline import load_config, read_manifest
from scripts.prepare_inputs import demux_filename


def run_chopper(
    input_path: Path, output_path: Path, options: list[str], threads: str,
) -> None:
    processes: list[subprocess.Popen[bytes]] = []
    try:
        with output_path.open("wb") as handle:
            zcat = subprocess.Popen(["zcat", str(input_path)], stdout=subprocess.PIPE)
            processes.append(zcat)
            assert zcat.stdout is not None
            chopper = subprocess.Popen(
                ["chopper", *options, "-t", threads],
                stdin=zcat.stdout, stdout=subprocess.PIPE,
            )
            processes.append(chopper)
            zcat.stdout.close()
            assert chopper.stdout is not None
            gzip_process = subprocess.Popen(
                ["gzip"], stdin=chopper.stdout, stdout=handle,
            )
            processes.append(gzip_process)
            chopper.stdout.close()
            for process in reversed(processes):
                process.wait()
            for process in reversed(processes):
                if process.returncode:
                    raise subprocess.CalledProcessError(process.returncode, process.args)
    finally:
        for process in processes:
            if process.stdout is not None:
                process.stdout.close()
            if process.poll() is None:
                process.terminate()
        for process in processes:
            process.wait()


def preprocess_reads(
    manifest_path: Path,
    config_dir: Path,
    demux_dir: Path,
    output_dir: Path,
    keep_intermediates: bool = False,
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

    work = output_dir / "work"
    for filename, reads in included:
        output = output_dir / filename
        sample_work = work / filename.removesuffix(".fastq.gz")
        sample_work.mkdir(parents=True)
        porechop_output = sample_work / "porechop.fastq.gz"
        subprocess.run(
            [
                "porechop_abi", *porechop_options,
                "-t", threads, "-i", str(reads), "-o", str(porechop_output),
                "--temp_dir", str(sample_work / "tmp"),
            ],
            check=True,
        )
        run_chopper(porechop_output, output, chopper_options, threads)
    if not keep_intermediates and work.exists():
        shutil.rmtree(work)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--demux-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Retain preprocessing work after success; failures always retain work",
    )
    args = parser.parse_args()
    preprocess_reads(
        args.manifest, args.config_dir, args.demux_dir, args.output_dir,
        args.keep_intermediates,
    )


if __name__ == "__main__":
    main()
