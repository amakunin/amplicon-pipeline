"""Fixture contracts and read-subsetting helper tests; no analysis jobs."""

import csv
import gzip
from pathlib import Path
import subprocess
import tempfile
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/coi_small"
BARCODES = ROOT / "assets/barcodes"
HELPERS = ROOT / "tests/helpers"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class FixtureContractTests(unittest.TestCase):
    def test_manifest_and_expected_barcodes(self) -> None:
        manifest = read_tsv(FIXTURE / "manifest.tsv")
        expected = read_tsv(FIXTURE / "expected_demux.tsv")
        ont = {
            row["tag_sequence"]: row["index"]
            for row in read_tsv(BARCODES / "ONT_nb24_barcode.tsv")
        }
        rear = {
            row["tag_sequence"]: row["index"]
            for row in read_tsv(BARCODES / "Sanger_custom_barcode_R.tsv")
        }
        self.assertEqual(len(manifest), 5)
        self.assertEqual(len(expected), 5)
        self.assertEqual(
            {row["sample_id"] for row in manifest},
            {row["sample_id"] for row in expected},
        )
        by_sample = {row["sample_id"]: row for row in expected}
        for sample in manifest:
            row = by_sample[sample["sample_id"]]
            self.assertEqual(ont[sample["ont_barcode"]], row["ont_barcode_id"])
            self.assertEqual(rear[sample["rear_barcode"]], row["rear_barcode_id"])
        self.assertEqual(
            {row["sample_id"] for row in expected if row["expected_reads"] == "positive"},
            {"ANLCx_001_A1", "ANLCx_001_C1", "ANLCx_001_E1"},
        )
        self.assertEqual(
            {row["sample_id"] for row in expected if row["expected_reads"] == "zero"},
            {"ANLCx_001_B1", "ANLCx_001_D1"},
        )

    def test_custom_barcode_assets_are_consistent(self) -> None:
        for end, count in (("F", 16), ("R", 96)):
            with self.subTest(end=end):
                prefix = BARCODES / f"Sanger_custom_barcode_{end}"
                lines = prefix.with_suffix(".fasta").read_text().splitlines()
                fasta = dict(zip(lines[::2], lines[1::2]))
                mapping = read_tsv(prefix.with_suffix(".tsv"))
                arrangement = tomllib.loads(prefix.with_suffix(".tmol").read_text())
                self.assertEqual(len(fasta), count)
                self.assertEqual(len(mapping), count)
                self.assertEqual(arrangement["arrangement"]["first_index"], 1)
                self.assertEqual(arrangement["arrangement"]["last_index"], count)
                for row in mapping:
                    name = f">SANGER-BC{end}{row['index']}"
                    self.assertEqual(fasta[name], row["tag_sequence"])

    def test_optional_real_fastq_integrity_and_selected_ids(self) -> None:
        fastq = FIXTURE / "run.fastq.gz"
        ids_path = FIXTURE / "read_ids.txt"
        if not fastq.exists():
            self.skipTest("Real FASTQ is not distributed with this checkout")
        self.assertTrue(ids_path.is_file(), "Record the fixture's selected read IDs")
        selected = set(ids_path.read_text().split())
        self.assertTrue(selected)
        observed: set[str] = set()
        with gzip.open(fastq, "rt") as handle:
            while header := handle.readline():
                sequence = handle.readline().rstrip("\n")
                plus = handle.readline()
                quality = handle.readline().rstrip("\n")
                self.assertTrue(header.startswith("@"))
                self.assertTrue(plus.startswith("+"))
                self.assertGreater(len(sequence), 0)
                self.assertEqual(len(sequence), len(quality))
                read_id = header.split()[0][1:]
                self.assertIn(read_id, selected)
                self.assertNotIn(read_id, observed)
                observed.add(read_id)
        self.assertEqual(observed, selected)


class ReadHelperTests(unittest.TestCase):
    def test_subsetting_plain_and_compressed_fastq(self) -> None:
        sequence = "ACATAAAGATATTGGT"
        record = f"@read1 metadata\n{sequence}\n+\n{'I' * len(sequence)}\n"
        other = "@read2\nACGT\n+\nIIII\n"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            ids = directory / "ids.txt"
            ids.write_text("read1\n")
            for compressed in (False, True):
                with self.subTest(compressed=compressed):
                    fastq = directory / ("reads.fastq.gz" if compressed else "reads.fastq")
                    if compressed:
                        with gzip.open(fastq, "wt") as handle:
                            handle.write(record + other)
                    else:
                        fastq.write_text(record + other)
                    extracted = subprocess.run(
                        ["bash", str(HELPERS / "extract_read_ids.sh"), str(fastq)],
                        capture_output=True, text=True, check=False,
                    )
                    self.assertEqual(extracted.returncode, 0, extracted.stderr)
                    self.assertEqual(extracted.stdout, "read1\n")
                    retrieved = subprocess.run(
                        ["bash", str(HELPERS / "retrieve_reads.sh"), str(ids), str(fastq)],
                        capture_output=True, text=True, check=False,
                    )
                    self.assertEqual(retrieved.returncode, 0, retrieved.stderr)
                    self.assertEqual(retrieved.stdout, record)
                    self.assertIn("Missing IDs: 0", retrieved.stderr)

    def test_empty_or_unrecoverable_ids_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fastq = directory / "reads.fastq"
            fastq.write_text("@read1\nACGT\n+\nIIII\n")
            ids = directory / "ids.txt"
            for contents in ("", "\n", "missing_read\n"):
                with self.subTest(contents=contents):
                    ids.write_text(contents)
                    result = subprocess.run(
                        ["bash", str(HELPERS / "retrieve_reads.sh"), str(ids), str(fastq)],
                        capture_output=True, text=True, check=False,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertTrue(result.stderr)
                    self.assertEqual(result.stdout, "")

    def test_usage_errors(self) -> None:
        for script in ("extract_read_ids.sh", "retrieve_reads.sh"):
            result = subprocess.run(
                ["bash", str(HELPERS / script)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("Usage:", result.stderr)


if __name__ == "__main__":
    unittest.main()
