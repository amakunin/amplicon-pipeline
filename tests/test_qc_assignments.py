"""Assignment QC annotations without external tools."""

import csv
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
import unittest

import yaml

from scripts.qc_assignments import annotate_assignments
from scripts.qc_asvs import annotate_asvs
from scripts.summarise_asvs import summarise_asvs


FIXTURE = Path(__file__).resolve().parents[1] / "tests/fixtures/small"


class AssignmentQcTests(unittest.TestCase):
    def test_boundaries_overrides_and_order_independence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = {
                "ampliseq": {"targets": {"COI": {}, "CytB": {}}},
                "summary": {
                    "common": {},
                    "targets": {"CytB": {
                        "min_contam_total_reads": 5000,
                        "max_contam_prop": 0.02,
                        "max_locov": 50,
                    }},
                },
            }
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(config))
            qc_path = directory / "asv_qc.tsv"
            with qc_path.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, delimiter="\t",
                    fieldnames=[
                        "target", "asv_id", "length", "length_status",
                        "stop_free_frames", "uncertain_frames", "coding_status",
                    ],
                )
                writer.writeheader()
                for target, asv_id, status in (
                    ("COI", "shared", "too_short"),
                    ("COI", "boundary", "pass"),
                    ("CytB", "shared", "pass"),
                ):
                    writer.writerow({
                        "target": target, "asv_id": asv_id, "length": 1,
                        "length_status": status, "coding_status": "not_checked",
                    })
            rows = [
                ("COI", "shared", "A", 9000),
                ("COI", "shared", "B", 9000),
                ("COI", "shared", "C", 1601),
                ("COI", "shared", "D", 200),
                ("COI", "shared", "E", 100),
                ("COI", "shared", "F", 99),
                ("COI", "boundary", "G", 9901),
                ("COI", "boundary", "H", 99),
                ("CytB", "shared", "I", 9901),
                ("CytB", "shared", "J", 99),
            ]
            assignments_path = directory / "assignments.tsv"
            self.write_assignments(assignments_path, rows)

            first = directory / "first"
            stdout = StringIO()
            with redirect_stdout(stdout):
                annotate_assignments(assignments_path, qc_path, directory, first)
            first_rows = self.read_assignments(first / "assignment_qc.tsv")
            expected = {
                "A": "", "B": "", "C": "", "D": "",
                "E": "contam_possible", "F": "contam_possible,locov",
                "G": "", "H": "locov", "I": "", "J": "contam_possible",
            }
            self.assertEqual(
                {row["sample_id"]: row["qc_labels"] for row in first_rows}, expected,
            )
            self.assertEqual(first_rows[0]["length_status"], "too_short")
            self.assertEqual(first_rows[0]["contam_total_reads"], "20000")
            self.assertEqual(first_rows[0]["contam_target_cutoff"], "200.0")
            report_header = (
                "target\tasv_id\tcontam_total_reads\tcontam_target_cutoff"
                "\taffected_samples\n"
            )
            report = (first / "assignment_qc.log").read_text()
            self.assertEqual(
                report, report_header
                + "COI\tshared\t20000\t200.0\t2\n"
                + "CytB\tshared\t10000\t200.0\t1\n",
            )
            self.assertEqual(
                stdout.getvalue(), report + f"Log: {first / 'assignment_qc.log'}\n",
            )

            self.write_assignments(assignments_path, list(reversed(rows)))
            second = directory / "second"
            annotate_assignments(assignments_path, qc_path, directory, second)
            self.assertEqual(
                {row["sample_id"]: row["qc_labels"] for row in
                 self.read_assignments(second / "assignment_qc.tsv")},
                expected,
            )
            self.assertEqual(
                (first / "assignment_qc.log").read_bytes(),
                (second / "assignment_qc.log").read_bytes(),
            )

    def test_fixture_annotations_preserve_raw_assignments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            summarise_asvs(
                {target: FIXTURE / "ampliseq" / target
                 for target in ("COI", "CytB", "ITS2", "P3")},
                output,
            )
            annotate_asvs(output / "asvs.tsv", FIXTURE / "config", output)
            before = (output / "assignments.tsv").read_bytes()
            annotate_assignments(
                output / "assignments.tsv", output / "asv_qc.tsv",
                FIXTURE / "config", output,
            )
            with (output / "assignment_qc.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                self.assertEqual(reader.fieldnames, [
                    "sample_id", "target", "asv_id", "read_count", "length",
                    "length_status", "stop_free_frames", "uncertain_frames",
                    "coding_status", "qc_labels", "contam_total_reads",
                    "contam_target_cutoff",
                ])
            self.assertEqual(len(rows), 7)
            self.assertEqual(
                [(int(row["read_count"]), row["qc_labels"]) for row in rows],
                [(150, ""), (12, "locov"), (9, "locov"), (11, "locov"),
                 (10, "locov"), (1, "locov"), (24, "locov")],
            )
            self.assertEqual((output / "assignments.tsv").read_bytes(), before)
            self.assertEqual(
                (output / "assignment_qc.log").read_text(),
                "target\tasv_id\tcontam_total_reads\tcontam_target_cutoff"
                "\taffected_samples\n",
            )

    def test_obsolete_settings_fail_instead_of_being_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = {
                "ampliseq": {"targets": {"COI": {}}},
                "summary": {
                    "common": {"min_contam_source_reads": 500},
                    "targets": {},
                },
            }
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(config))
            with self.assertRaisesRegex(ValueError, "min_contam_source_reads"):
                annotate_assignments(
                    directory / "assignments.tsv", directory / "asv_qc.tsv",
                    directory, directory / "summary",
                )

    @staticmethod
    def write_assignments(
        path: Path, assignments: list[tuple[str, str, str, int]],
    ) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", "target", "asv_id", "read_count"])
            writer.writerows(
                (sample, target, asv_id, count)
                for target, asv_id, sample, count in assignments
            )

    @staticmethod
    def read_assignments(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))


if __name__ == "__main__":
    unittest.main()
