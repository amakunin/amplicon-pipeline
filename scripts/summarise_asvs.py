#!/usr/bin/env python3
"""Aggregate pooled Savont/SINTAX evidence and annotate QC without filtering."""

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


RAW_COLUMNS = [
    "ASV_ID", "sintax_raw_taxonomy", "sintax_strand", "sintax_filtered_taxonomy",
]
SEQUENCE_QC_COLUMNS = [
    "length_status", "stop_free_frames", "uncertain_frames", "coding_status",
]


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


def annotate_asvs(
    asvs: list[dict[str, str | int]], settings: dict[str, dict[str, Any]],
) -> None:
    for target, parameters in settings.items():
        if "genetic_code" in parameters and parameters["genetic_code"] != 5:
            raise ValueError(f"Only genetic code 5 is supported: {target}")
    for row in asvs:
        target = str(row["target"])
        parameters = settings[target]
        sequence = str(row["sequence"])
        if not sequence:
            raise ValueError(f"Empty sequence: {target}, {row['asv_id']}")
        length = len(sequence)
        row["length"] = length
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
        row.update({
            "length_status": length_status,
            "stop_free_frames": ",".join(stop_free),
            "uncertain_frames": ",".join(uncertain),
            "coding_status": coding_status,
        })


def annotate_assignments(
    asvs: list[dict[str, str | int]],
    assignments: list[dict[str, str | int]],
    settings: dict[str, dict[str, Any]],
) -> list[str]:
    for parameters in settings.values():
        obsolete = {"min_contam_source_reads", "min_contam_read_ratio"} & parameters.keys()
        if obsolete:
            raise ValueError(f"Obsolete contamination settings: {', '.join(sorted(obsolete))}")
    grouped: dict[tuple[str, str], list[dict[str, str | int]]] = {}
    for row in assignments:
        key = (str(row["target"]), str(row["asv_id"]))
        grouped.setdefault(key, []).append(row)
    report = []
    for asv in sorted(asvs, key=lambda row: (str(row["target"]), str(row["asv_id"]))):
        target, asv_id = str(asv["target"]), str(asv["asv_id"])
        rows = grouped.get((target, asv_id), [])
        parameters = settings[target]
        total = sum(int(row["read_count"]) for row in rows)
        proportion = parameters.get("max_contam_prop", 0.01)
        cutoff = total * proportion
        asv["total_reads"] = total
        asv["contam_target_cutoff"] = str(cutoff)
        affected = 0
        for row in rows:
            count = int(row["read_count"])
            contam = total > parameters.get("min_contam_total_reads", 10000) and count / total < proportion
            row["contam_possible"] = "true" if contam else "false"
            row["locov"] = "true" if count < parameters.get("max_locov", 100) else "false"
            affected += int(contam)
        if affected:
            report.append(
                f"Contamination: {target}, {asv_id}; total reads: {total}; "
                f"cutoff: {cutoff}; affected samples: {affected}"
            )
    return report


def read_asv_table(
    path: Path, fieldnames: list[str] | None = None,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=fieldnames, delimiter="\t")
        columns = list(reader.fieldnames or [])
        if not columns or columns[0] != "ASV_ID" or len(set(columns)) != len(columns):
            raise ValueError(f"Expected unique columns starting with ASV_ID: {path}")
        rows = {}
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"Wrong number of fields: {path}, line {reader.line_num}")
            asv_id = row["ASV_ID"]
            if not asv_id or asv_id in rows:
                raise ValueError(f"Empty or duplicate ASV_ID {asv_id!r}: {path}")
            rows[asv_id] = row
    return columns, rows


def summarise_asvs(
    target_dirs: dict[str, Path], output_dir: Path, config_dir: Path,
) -> None:
    with (config_dir / "pipeline.yaml").open() as handle:
        config = yaml.safe_load(handle)
    summary = config["summary"]
    settings = {}
    for target in target_dirs:
        config["ampliseq"]["targets"][target]
        settings[target] = {**summary["common"], **summary["targets"].get(target, {})}
    output_dir = output_dir.resolve()
    report: list[str] = []
    ranks: list[str] = []
    asvs: list[dict[str, str | int]] = []
    assignments: list[dict[str, str | int]] = []
    for target, results_dir in target_dirs.items():
        results_dir = results_dir.resolve()
        feature_path = results_dir / "savont/pooled_feature-table.tsv"
        sample_columns, features = read_asv_table(feature_path)
        taxonomy_paths = [
            path for path in sorted((results_dir / "sintax").glob("ASV_tax_sintax.*.tsv"))
            if not path.name.endswith(".raw.tsv")
        ]
        if len(taxonomy_paths) != 1:
            raise ValueError(f"Expected one processed SINTAX table: {results_dir / 'sintax'}")
        taxonomy_path = taxonomy_paths[0]
        raw_path = taxonomy_path.with_suffix(".raw.tsv")
        report.extend([
            f"Target {target}: {results_dir}",
            f"  Counts: {feature_path.resolve()}",
            f"  Taxonomy: {taxonomy_path.resolve()}",
            f"  Raw SINTAX: {raw_path.resolve()}",
        ])
        columns, taxonomy = read_asv_table(taxonomy_path)
        _, raw = read_asv_table(raw_path, RAW_COLUMNS)
        if features.keys() != taxonomy.keys() or features.keys() != raw.keys():
            raise ValueError(f"Inconsistent ASV IDs across counts and SINTAX: {target}")
        for column in columns:
            if column not in ("ASV_ID", "sequence", "confidence") and column not in ranks:
                ranks.append(column)
        assignment_count = 0
        read_count = 0
        for asv_id, feature in features.items():
            annotation = taxonomy[asv_id]
            sequence = annotation["sequence"]
            if not sequence:
                raise ValueError(f"Empty sequence: {target}, {asv_id}")
            asvs.append({
                **{column: value for column, value in annotation.items()
                   if column != "ASV_ID"},
                **{column: raw[asv_id][column] for column in RAW_COLUMNS[1:]},
                "target": target,
                "asv_id": asv_id,
                "confidence": annotation["confidence"],
                "length": len(sequence),
            })
            for sample_id in sample_columns[1:]:
                count = int(feature[sample_id])
                if count < 0:
                    raise ValueError(f"Negative read count: {target}, {asv_id}, {sample_id}")
                if count:
                    assignment_count += 1
                    read_count += count
                    assignments.append({
                        "sample_id": sample_id, "target": target,
                        "asv_id": asv_id, "read_count": count,
                    })
        report.append(
            f"  Unique ASVs: {len(features)}; non-zero sample-ASV assignments: "
            f"{assignment_count}; assigned reads: {read_count}"
        )

    annotate_asvs(asvs, settings)
    report.extend(annotate_assignments(asvs, assignments, settings))
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, fields, rows in (
        ("asvs.tsv", ["target", "asv_id", "sequence", "length", *ranks,
                      "confidence", *RAW_COLUMNS[1:], *SEQUENCE_QC_COLUMNS,
                      "total_reads", "contam_target_cutoff"], asvs),
        ("assignments.tsv", ["sample_id", "target", "asv_id", "read_count",
                             "contam_possible", "locov"], assignments),
    ):
        with (output_dir / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        report.append(f"Output: {output_dir / filename}")
    log_path = output_dir / "summarise_asvs.log"
    report.append(f"Log: {log_path}")
    text = "\n".join(report) + "\n"
    log_path.write_text(text)
    print(text, end="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target", nargs=2, action="append", required=True,
        metavar=("TARGET", "RESULTS_DIR"),
        help="Target name and published ampliseq directory; repeat for each target",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, required=True)
    args = parser.parse_args()
    target_dirs = {}
    for target, directory in args.target:
        if target in target_dirs:
            raise ValueError(f"Duplicate target: {target}")
        target_dirs[target] = Path(directory)
    summarise_asvs(target_dirs, args.output_dir, args.config_dir)


if __name__ == "__main__":
    main()
