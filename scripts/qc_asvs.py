#!/usr/bin/env python3
"""Annotate ASV lengths and six-frame coding evidence without filtering."""

import argparse
import csv
from pathlib import Path

import yaml


def check_frames(sequence: str) -> tuple[list[str], list[str]]:
    sequence = sequence.upper()
    reverse = sequence.translate(
        str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")
    )[::-1]
    stop_free: list[str] = []
    uncertain: list[str] = []
    for strand, oriented in (("+", sequence), ("-", reverse)):
        for shift in range(3):
            codons = [
                oriented[start:start + 3]
                for start in range(shift, len(oriented) - 2, 3)
            ]
            # NCBI table 5: TGA encodes tryptophan, not a stop.
            if any(codon in ("TAA", "TAG") for codon in codons):
                continue
            frame = f"{strand}{shift}"
            if any(set(codon) - set("ACGT") for codon in codons):
                uncertain.append(frame)
            else:
                stop_free.append(frame)
    return stop_free, uncertain


def annotate_asvs(asvs_path: Path, config_dir: Path, output_dir: Path) -> None:
    with (config_dir / "pipeline.yaml").open() as handle:
        config = yaml.safe_load(handle)
    summary = config["summary"]
    settings = {}
    for target in config["ampliseq"]["targets"]:
        parameters = {**summary["common"], **summary["targets"].get(target, {})}
        if "genetic_code" in parameters and parameters["genetic_code"] != 5:
            raise ValueError(f"Only genetic code 5 is supported: {target}")
        settings[target] = parameters

    annotations = []
    seen = set()
    with asvs_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target, asv_id = row["target"], row["asv_id"]
            key = (target, asv_id)
            if key in seen:
                raise ValueError(f"Duplicate target-ASV: {target}, {asv_id}")
            seen.add(key)
            parameters = settings[target]
            sequence = row["sequence"]
            if not sequence:
                raise ValueError(f"Empty sequence: {target}, {asv_id}")
            length = len(sequence)
            lower_bounds = []
            upper_bounds = []
            if "allowed_lengths" in parameters:
                lower_bounds.append(min(parameters["allowed_lengths"]))
                upper_bounds.append(max(parameters["allowed_lengths"]))
            if "min_length" in parameters:
                lower_bounds.append(parameters["min_length"])
            if "max_length" in parameters:
                upper_bounds.append(parameters["max_length"])
            if lower_bounds and upper_bounds and max(lower_bounds) > min(upper_bounds):
                raise ValueError(f"Conflicting length bounds: {target}")
            if any(length < bound for bound in lower_bounds):
                length_status = "too_short"
            elif any(length > bound for bound in upper_bounds):
                length_status = "too_long"
            elif "allowed_lengths" in parameters and length not in parameters["allowed_lengths"]:
                length_status = "unexpected_length"
            else:
                length_status = "pass" if lower_bounds or upper_bounds else "not_checked"
            stop_free, uncertain = [], []
            coding_status = "not_checked"
            if "genetic_code" in parameters:
                stop_free, uncertain = check_frames(sequence)
                if uncertain:
                    coding_status = "uncertain"
                elif not stop_free:
                    coding_status = "no_stop_free_frame"
                elif len(stop_free) == 1:
                    coding_status = "one_stop_free_frame"
                else:
                    coding_status = "multiple_stop_free_frames"
            annotations.append({
                "target": target, "asv_id": asv_id, "length": length,
                "length_status": length_status,
                "stop_free_frames": ",".join(stop_free),
                "uncertain_frames": ",".join(uncertain),
                "coding_status": coding_status,
            })

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "asv_qc.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
            "target", "asv_id", "length", "length_status",
            "stop_free_frames", "uncertain_frames", "coding_status",
        ])
        writer.writeheader()
        writer.writerows(annotations)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asvs", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    annotate_asvs(args.asvs, args.config_dir, args.output_dir)


if __name__ == "__main__":
    main()
