#!/usr/bin/env python3
"""Aggregate pooled Savont counts and SINTAX evidence without filtering."""

import argparse
import csv
from pathlib import Path


RAW_COLUMNS = [
    "ASV_ID", "sintax_raw_taxonomy", "sintax_strand", "sintax_filtered_taxonomy",
]


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


def summarise_asvs(target_dirs: dict[str, Path], output_dir: Path) -> None:
    ranks: list[str] = []
    asvs: list[dict[str, str | int]] = []
    assignments: list[dict[str, str | int]] = []
    for target, results_dir in target_dirs.items():
        sample_columns, features = read_asv_table(
            results_dir / "savont/pooled_feature-table.tsv"
        )
        taxonomy_paths = [
            path for path in sorted((results_dir / "sintax").glob("ASV_tax_sintax.*.tsv"))
            if not path.name.endswith(".raw.tsv")
        ]
        if len(taxonomy_paths) != 1:
            raise ValueError(f"Expected one processed SINTAX table: {results_dir / 'sintax'}")
        taxonomy_path = taxonomy_paths[0]
        columns, taxonomy = read_asv_table(taxonomy_path)
        _, raw = read_asv_table(
            taxonomy_path.with_suffix(".raw.tsv"), RAW_COLUMNS
        )
        if features.keys() != taxonomy.keys() or features.keys() != raw.keys():
            raise ValueError(f"Inconsistent ASV IDs across counts and SINTAX: {target}")
        for column in columns:
            if column not in ("ASV_ID", "sequence", "confidence") and column not in ranks:
                ranks.append(column)
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
                    assignments.append({
                        "sample_id": sample_id, "target": target,
                        "asv_id": asv_id, "read_count": count,
                    })

    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, fields, rows in (
        ("asvs.tsv", ["target", "asv_id", "sequence", "length", *ranks,
                      "confidence", *RAW_COLUMNS[1:]], asvs),
        ("assignments.tsv", ["sample_id", "target", "asv_id", "read_count"], assignments),
    ):
        with (output_dir / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target", nargs=2, action="append", required=True,
        metavar=("TARGET", "RESULTS_DIR"),
        help="Target name and published ampliseq directory; repeat for each target",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    target_dirs = {}
    for target, directory in args.target:
        if target in target_dirs:
            raise ValueError(f"Duplicate target: {target}")
        target_dirs[target] = Path(directory)
    summarise_asvs(target_dirs, args.output_dir)


if __name__ == "__main__":
    main()
