#!/usr/bin/env python3
"""Annotate sample-ASV assignments with sequence QC and assignment evidence."""

import argparse
import csv
from pathlib import Path

import yaml


SEQUENCE_QC_COLUMNS = [
    "length", "length_status", "stop_free_frames", "uncertain_frames",
    "coding_status",
]


def annotate_assignments(
    assignments_path: Path,
    asv_qc_path: Path,
    config_dir: Path,
    output_dir: Path,
) -> None:
    with (config_dir / "pipeline.yaml").open() as handle:
        config = yaml.safe_load(handle)
    summary = config["summary"]
    settings = {}
    for target in config["ampliseq"]["targets"]:
        settings[target] = {
            **summary["common"],
            **summary["targets"].get(target, {}),
        }
        obsolete = {"min_contam_source_reads", "min_contam_read_ratio"} & settings[target].keys()
        if obsolete:
            raise ValueError(f"Obsolete contamination settings: {', '.join(sorted(obsolete))}")

    sequence_qc = {}
    with asv_qc_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequence_qc[(row["target"], row["asv_id"])] = row

    with assignments_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assignment_columns = list(reader.fieldnames or [])
        assignments = list(reader)

    grouped = {}
    for row in assignments:
        key = (row["target"], row["asv_id"])
        row.update({
            column: sequence_qc[key][column]
            for column in SEQUENCE_QC_COLUMNS
        })
        grouped.setdefault(key, []).append(row)

    report_rows = []
    for (target, asv_id), rows in grouped.items():
        parameters = settings[target]
        min_total_reads = parameters.get("min_contam_total_reads", 10000)
        max_contam_prop = parameters.get("max_contam_prop", 0.01)
        max_locov = parameters.get("max_locov", 100)
        total_reads = sum(int(row["read_count"]) for row in rows)
        cutoff = total_reads * max_contam_prop
        affected = 0
        for row in rows:
            count = int(row["read_count"])
            labels = []
            if total_reads > min_total_reads and count / total_reads < max_contam_prop:
                labels.append("contam_possible")
                affected += 1
            if count < max_locov:
                labels.append("locov")
            row["qc_labels"] = ",".join(labels)
            row["contam_total_reads"] = str(total_reads)
            row["contam_target_cutoff"] = str(cutoff)
        if affected:
            report_rows.append({
                "target": target,
                "asv_id": asv_id,
                "contam_total_reads": str(total_reads),
                "contam_target_cutoff": str(cutoff),
                "affected_samples": affected,
            })

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "assignment_qc.tsv"
    fields = [
        *assignment_columns, *SEQUENCE_QC_COLUMNS, "qc_labels",
        "contam_total_reads", "contam_target_cutoff",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(assignments)

    report_path = output_dir / "assignment_qc.log"
    report_fields = [
        "target", "asv_id", "contam_total_reads", "contam_target_cutoff",
        "affected_samples",
    ]
    report_rows.sort(key=lambda row: (row["target"], row["asv_id"]))
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=report_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)
    print(report_path.read_text(), end="")
    print(f"Log: {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--asv-qc", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    annotate_assignments(
        args.assignments, args.asv_qc, args.config_dir, args.output_dir,
    )


if __name__ == "__main__":
    main()
