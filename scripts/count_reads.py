#!/usr/bin/env python3
"""Report hierarchical read counts without filtering or running analysis."""

import argparse
import csv
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_pipeline import load_config, read_manifest
from scripts.preprocess_reads import count_fastq_reads
from scripts.prepare_inputs import demux_filename
from scripts.summarise_asvs import read_asv_table


def read_rows(
    path: Path, keys: tuple[str, ...], required: tuple[str, ...] = (),
) -> dict[tuple[str, ...], dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for column in (*keys, *required):
            if column not in (reader.fieldnames or []):
                raise KeyError(f"Missing column {column}: {path}")
        rows = {}
        for row in reader:
            key = tuple(row[column] for column in keys)
            if key in rows:
                raise ValueError(f"Duplicate observation {key}: {path}")
            rows[key] = row
    return rows


def count_reads(
    run_fastq: Path, manifest_path: Path, config_dir: Path, output_dir: Path,
) -> None:
    config, targets = load_config(config_dir)
    samples = read_manifest(manifest_path, config)
    counts = read_rows(
        output_dir / "preprocessed/read_counts.tsv", ("sample_id",),
        ("input_reads", "porechop_reads", "chopper_reads"),
    )
    asvs = read_rows(
        output_dir / "summary/asvs.tsv", ("target", "asv_id"),
        ("length_status", "coding_status", "total_reads"),
    )
    assignments = read_rows(
        output_dir / "summary/assignments.tsv", ("sample_id", "target", "asv_id"),
        ("read_count", "contam_possible", "locov"),
    )
    summaries = {}
    features = {}
    for target in targets:
        directory = output_dir / "ampliseq" / target / "results"
        summaries[target] = read_rows(
            directory / "overall_summary.tsv", ("sample",),
            ("cutadapt_total_processed", "cutadapt_passing_filters", "savont_output"),
        )
        columns, table = read_asv_table(directory / "savont/pooled_feature-table.tsv")
        features[target] = {
            sample: sum(int(row[sample]) for row in table.values())
            for sample in columns[1:]
        }
    totals: dict[tuple[str, str], list[int]] = {}
    asv_totals: dict[tuple[str, str], int] = {}
    length_pass = {"pass", "not_checked"}
    coding_pass = {"not_checked", "one_stop_free_frame", "multiple_stop_free_frames"}
    length_statuses = length_pass | {"too_short", "too_long", "unexpected_length"}
    coding_statuses = coding_pass | {"uncertain", "no_stop_free_frame"}
    for row in assignments.values():
        asv = asvs[(row["target"], row["asv_id"])]
        if asv["length_status"] not in length_statuses or asv["coding_status"] not in coding_statuses:
            raise ValueError(f"Unknown sequence QC status: {asv}")
        if any(row[key] not in {"true", "false"} for key in ("contam_possible", "locov")):
            raise ValueError(f"Invalid assignment QC flag: {row}")
        reads = int(row["read_count"])
        key = (row["target"], row["asv_id"])
        asv_totals[key] = asv_totals.get(key, 0) + reads
        sequence_pass = asv["length_status"] in length_pass and asv["coding_status"] in coding_pass
        both_pass = sequence_pass and row["contam_possible"] == row["locov"] == "false"
        total = totals.setdefault((row["sample_id"], row["target"]), [0, 0, 0])
        total[0] += reads
        total[1] += reads if sequence_pass else 0
        total[2] += reads if both_pass else 0

    fields = [
        "node_id", "parent_id", "stage", "first_barcode", "second_barcode",
        "sample_id", "target", "label", "read_count",
    ]
    report = []

    def node(parent: str, stage: str, reads: str | int, first: str = "",
             second: str = "", sample: str = "", target: str = "", label: str = "") -> str:
        identifier = str(len(report) + 1)
        report.append(dict(zip(fields, (
            identifier, parent, stage, first, second, sample, target, label, str(reads),
        ))))
        return identifier

    def number(value: str | int) -> int | None:
        return None if value in ("", "NA") else int(value)

    def reconcile(label: str, before: str | int, after: str | int, equality: bool = True) -> None:
        a, b = number(before), number(after)
        if a is not None and b is not None and (a != b if equality else b > a):
            print(f"Warning: {label}: {before} -> {after}", file=sys.stderr)

    for key, asv in asvs.items():
        reconcile(f"{key} ASV/assignment total", asv["total_reads"], asv_totals.get(key, 0))
    pairs = {}
    for sample in samples:
        name = demux_filename(config["demux"]["method"], sample["barcode_indices"])
        pairs[name.removesuffix(".fastq.gz")] = sample["sample_id"]
    prefix = "ont_barcode" if config["demux"]["method"] == "hybrid" else "barcode"
    pattern = re.compile(rf"({prefix}[0-9]+)_(barcode[0-9]+|unclassified)")
    branches: dict[str, dict[str, int | str]] = {}
    demux_total = 0
    with (output_dir / "demux/demux_stats.txt").open(newline="") as handle:
        seen = set()
        for name, value in csv.reader(handle, delimiter="\t"):
            if name in seen:
                raise ValueError(f"Duplicate demux observation: {name}")
            seen.add(name)
            reads = int(value)
            demux_total += reads
            match = pattern.fullmatch(name)
            if name == "unclassified":
                first, second = "unclassified", ""
            elif match:
                first, second = match.groups()
            else:
                raise ValueError(f"Unexpected demux stats name: {name}")
            branches.setdefault(first, {})[second] = reads
    for name in pairs:
        match = pattern.fullmatch(name)
        assert match is not None
        first, second = match.groups()
        branches.setdefault(first, {}).setdefault(second, "")
    raw = count_fastq_reads(run_fastq)
    run = node("", "run", raw)
    reconcile("raw/demux total", raw, demux_total)
    expected_first = {name.split("_barcode")[0] for name in pairs}
    for first, children in sorted(branches.items()):
        label = "unmatched" if first == "unclassified" else (
            "expected" if first in expected_first else "unexpected"
        )
        observed_children = [int(v) for v in children.values() if v != ""]
        first_node = node(
            run, "first_barcode", sum(observed_children) if observed_children else "",
            first=first, label=label,
        )
        for second, reads in sorted(children.items()):
            if first == "unclassified":
                continue
            name = f"{first}_{second}"
            sample = pairs.get(name, "")
            label = "unmatched" if second == "unclassified" else ("expected" if sample else "unexpected")
            pair_node = node(first_node, "second_barcode", reads, first, second, label=label)
            if not sample:
                continue
            measured = counts[(sample,)]
            parent = node(pair_node, "sample", measured["input_reads"], first, second, sample)
            reconcile(f"{sample} demux/input", reads, measured["input_reads"])
            previous = measured["input_reads"]
            for stage, key in (("porechop", "porechop_reads"), ("chopper", "chopper_reads")):
                reconcile(f"{sample} {stage}", previous, measured[key], equality=False)
                parent = node(parent, stage, measured[key], first, second, sample)
                previous = measured[key]
            primer_inputs = []
            for target in targets:
                target_node = node(parent, "target", "", first, second, sample, target)
                summary = summaries[target].get((sample,), {})
                primer_input = summary.get("cutadapt_total_processed", "")
                reconcile(f"{sample}/{target} Chopper/primer input", previous, primer_input)
                primer_inputs.append((target, primer_input))
                current = target_node
                last = primer_input
                for stage, key in (
                    ("primer_input", "cutadapt_total_processed"),
                    ("primer_passing", "cutadapt_passing_filters"),
                    ("savont_output", "savont_output"),
                ):
                    value = summary.get(key, "")
                    reconcile(f"{sample}/{target} {stage}", last, value, equality=False)
                    current = node(current, stage, value, first, second, sample, target)
                    last = value
                feature_count = features[target].get(sample, "")
                reconcile(f"{sample}/{target} Savont/features", last, feature_count)
                current = node(current, "assigned_reads", feature_count, first, second, sample, target)
                observed = totals.get(
                    (sample, target), [0, 0, 0] if feature_count != "" else ["", "", ""],
                )
                reconcile(f"{sample}/{target} features/assignments", feature_count, observed[0])
                for stage, value in zip(("sequence_qc_pass", "both_qc_pass"), observed[1:]):
                    current = node(current, stage, value, first, second, sample, target)
            for index, (target_a, a) in enumerate(primer_inputs):
                for target_b, b in primer_inputs[index + 1:]:
                    reconcile(f"{sample} primer inputs {target_a}/{target_b}", a, b)
    with (output_dir / "summary/read_flow.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(report)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-fastq", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    count_reads(args.run_fastq, args.manifest, args.config_dir, args.output_dir)


if __name__ == "__main__":
    main()
