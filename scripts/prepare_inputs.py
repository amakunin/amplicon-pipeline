#!/usr/bin/env python3
"""Write the ampliseq samplesheet, excluded-sample list and parameter file for one target."""

import argparse
import csv
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_pipeline import load_config, read_manifest


def demux_filename(method: str, indices: dict[str, str]) -> str:
    first = "ont_barcode" if method == "hybrid" else "barcode"
    first_index, rear_index = indices.values()
    return f"{first}{first_index}_barcode{rear_index}.fastq.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--preprocessed-dir", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    config, targets = load_config(args.config_dir)
    samples = read_manifest(args.manifest, config)
    method = config["demux"]["method"]
    preprocessed_dir = args.preprocessed_dir.resolve()
    with (preprocessed_dir / "excluded_samples.tsv").open(newline="") as handle:
        upstream_exclusions = {
            row["sample_id"]: row["reason"]
            for row in csv.DictReader(handle, delimiter="\t")
        }

    included, excluded = [], []
    for sample in samples:
        sample_id = sample["sample_id"]
        if sample_id in upstream_exclusions:
            excluded.append((sample_id, upstream_exclusions[sample_id]))
            continue
        reads = preprocessed_dir / demux_filename(
            method, sample["barcode_indices"]
        )
        reads.stat()
        included.append((sample_id, reads))

    target_dir = args.output_dir.resolve() / args.target
    target_dir.mkdir(parents=True, exist_ok=True)

    with (target_dir / "samplesheet.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "fastq_1"])
        writer.writerows(included)

    with (target_dir / "excluded_samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "reason"])
        writer.writerows(excluded)

    parameters = {
        **targets[args.target],
        "input": str(target_dir / "samplesheet.tsv"),
        "outdir": str(target_dir / "results"),
    }
    with (target_dir / "params.yaml").open("w") as handle:
        yaml.safe_dump(parameters, handle, sort_keys=False)


if __name__ == "__main__":
    main()
