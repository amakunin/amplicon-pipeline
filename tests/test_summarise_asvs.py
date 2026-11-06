"""ASV aggregation from frozen outputs; no external analysis tools."""

import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from scripts.summarise_asvs import summarise_asvs


ROOT = Path(__file__).resolve().parents[1]
CAPTURED = ROOT / "tests/fixtures/small/ampliseq"
CONFIG = ROOT / "tests/fixtures/small/config"


class SummariseAsvsTests(unittest.TestCase):
    def test_captured_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "summary"
            command = [sys.executable, str(ROOT / "scripts/summarise_asvs.py")]
            for target in ("COI", "CytB", "ITS2", "P3"):
                command += ["--target", target, str(CAPTURED / target)]
            command += ["--output-dir", str(output), "--config-dir", str(CONFIG)]
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertEqual(result.stderr, "")
            self.assertEqual(result.stdout, (output / "summarise_asvs.log").read_text())
            for target, unique, assignments, reads in (
                ("COI", 0, 0, 0), ("CytB", 1, 2, 162),
                ("ITS2", 1, 2, 20), ("P3", 2, 3, 35),
            ):
                database = "coidb" if target == "COI" else "user"
                directory = CAPTURED / target
                self.assertIn(f"Target {target}: {directory.resolve()}\n", result.stdout)
                for label, path in (
                    ("Counts", directory / "savont/pooled_feature-table.tsv"),
                    ("Taxonomy", directory / f"sintax/ASV_tax_sintax.{database}.tsv"),
                    ("Raw SINTAX", directory / f"sintax/ASV_tax_sintax.{database}.raw.tsv"),
                ):
                    self.assertIn(f"  {label}: {path.resolve()}\n", result.stdout)
                self.assertIn(
                    f"  Unique ASVs: {unique}; non-zero sample-ASV assignments: "
                    f"{assignments}; assigned reads: {reads}\n",
                    result.stdout,
                )
            for filename in ("asvs.tsv", "assignments.tsv"):
                self.assertIn(f"Output: {output / filename}\n", result.stdout)
            self.assertIn(f"Log: {output / 'summarise_asvs.log'}\n", result.stdout)
            with (output / "asvs.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                asvs = list(reader)
                self.assertEqual(reader.fieldnames, [
                    "target", "asv_id", "sequence", "length", "Kingdom", "Phylum",
                    "Class", "Order", "Family", "Genus", "Species", "BOLD_bin",
                    "Domain", "Subgenus", "Group", "confidence",
                    "sintax_raw_taxonomy", "sintax_strand", "sintax_filtered_taxonomy",
                    "length_status", "stop_free_frames", "uncertain_frames", "coding_status",
                    "total_reads", "contam_target_cutoff",
                ])
            self.assertEqual(
                [(row["target"], row["asv_id"], int(row["length"])) for row in asvs],
                [("CytB", "0aece8f365778d45", 307),
                 ("ITS2", "71e9be273051fbdb", 527),
                 ("P3", "0a4b74ea4eddec52", 371),
                 ("P3", "83ea72fd937a3d48", 373)],
            )
            for target in ("CytB", "ITS2", "P3"):
                path = CAPTURED / target / "sintax/ASV_tax_sintax.user.tsv"
                with path.open(newline="") as handle:
                    source = list(csv.DictReader(handle, delimiter="\t"))
                with path.with_suffix(".raw.tsv").open(newline="") as handle:
                    raw = {row[0]: row[1:] for row in csv.reader(handle, delimiter="\t")}
                for original in source:
                    row = next(row for row in asvs if
                               (row["target"], row["asv_id"]) == (target, original["ASV_ID"]))
                    for column, value in original.items():
                        if column != "ASV_ID":
                            self.assertEqual(row[column], value)
                    self.assertEqual(
                        [row[column] for column in (
                            "sintax_raw_taxonomy", "sintax_strand",
                            "sintax_filtered_taxonomy",
                        )],
                        raw[row["asv_id"]],
                    )
            self.assertEqual(asvs[1]["Species"], "")
            self.assertIn("Anopheles coluzzii(0.59)", asvs[1]["sintax_raw_taxonomy"])
            self.assertEqual(asvs[0]["Domain"], "")
            self.assertEqual(asvs[2]["Kingdom"], "")
            with (output / "assignments.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, [
                    "sample_id", "target", "asv_id", "read_count", "contam_possible", "locov",
                ])
                observed = [
                    (row["target"], row["asv_id"], row["sample_id"], int(row["read_count"]))
                    for row in reader
                ]
            self.assertEqual(observed, [
                ("CytB", "0aece8f365778d45", "BFHK_023_A2", 150),
                ("CytB", "0aece8f365778d45", "BFHK_023_H8", 12),
                ("ITS2", "71e9be273051fbdb", "BFHK_023_A2", 9),
                ("ITS2", "71e9be273051fbdb", "BFHK_023_H8", 11),
                ("P3", "0a4b74ea4eddec52", "ANLCx_002_A12", 10),
                ("P3", "0a4b74ea4eddec52", "BFHK_023_H8", 1),
                ("P3", "83ea72fd937a3d48", "BFHK_023_A2", 24),
            ])

    def test_empty_target_writes_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "summary"
            summarise_asvs({"COI": CAPTURED / "COI"}, output, CONFIG)
            for filename in ("asvs.tsv", "assignments.tsv"):
                with (output / filename).open(newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    self.assertEqual(list(reader), [])
                    self.assertIn("target", reader.fieldnames)

    def test_inconsistent_or_duplicate_ids_do_not_replace_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            target = directory / "CytB"
            shutil.copytree(CAPTURED / "CytB", target)
            raw_path = target / "sintax/ASV_tax_sintax.user.raw.tsv"
            original = raw_path.read_text()
            output = directory / "summary"
            summarise_asvs({"CytB": target}, output, CONFIG)
            before = {path.name: path.read_bytes() for path in output.iterdir()}
            for content, message in (
                (original.replace("0aece8f365778d45", "different"), "Inconsistent ASV IDs"),
                (original + original, "duplicate ASV_ID"),
            ):
                with self.subTest(message=message):
                    raw_path.write_text(content)
                    with self.assertRaisesRegex(ValueError, message):
                        summarise_asvs({"CytB": target}, output, CONFIG)
                    self.assertEqual(
                        {path.name: path.read_bytes() for path in output.iterdir()}, before
                    )


if __name__ == "__main__":
    unittest.main()
