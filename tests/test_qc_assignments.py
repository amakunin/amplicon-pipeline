"""Assignment contamination annotations without external tools."""

import csv
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
import unittest

import yaml

from scripts.qc_assignments import annotate_assignments


class AssignmentQcTests(unittest.TestCase):
    def test_contamination_labels_boundaries_and_order_independence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = {
                "ampliseq": {"targets": {"COI": {}}},
                "summary": {
                    "common": {
                        "min_contam_source_reads": 10000,
                        "min_contam_read_ratio": 100,
                    },
                    "targets": {},
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
                writer.writerow({
                    "target": "COI", "asv_id": "source", "length": 1,
                    "length_status": "too_short", "coding_status": "not_checked",
                })
                writer.writerow({
                    "target": "COI", "asv_id": "boundary", "length": 2,
                    "length_status": "pass", "coding_status": "not_checked",
                })
            rows = [
                ("source", "A", 2_000_000),
                ("source", "B", 2_000_000),
                ("source", "C", 15_000),
                ("source", "D", 19_999),
                ("source", "E", 20_000),
                ("source", "F", 10_000),
                ("boundary", "G", 10_000),
            ]
            assignments_path = directory / "assignments.tsv"
            self.write_assignments(assignments_path, rows)

            first = directory / "first"
            stdout = StringIO()
            with redirect_stdout(stdout):
                annotate_assignments(assignments_path, qc_path, directory, first)
            expected = {
                "A": "contam_source",
                "B": "contam_source",
                "C": "possible_contam_source,contam_target",
                "D": "possible_contam_source,contam_target",
                "E": "possible_contam_source",
                "F": "contam_target",
                "G": "",
            }
            first_rows = self.read_assignments(first / "assignment_qc.tsv")
            self.assertEqual(
                {row["sample_id"]: row["contam_labels"] for row in first_rows},
                expected,
            )
            self.assertEqual(first_rows[0]["length_status"], "too_short")
            report_header = (
                "target\tasv_id\tcontam_source_reads\tcontam_target_cutoff"
                "\taffected_samples\n"
            )
            report = (first / "assignment_qc.log").read_text()
            self.assertEqual(
                report, report_header + "COI\tsource\t2000000\t20000.0\t3\n",
            )
            self.assertEqual(
                stdout.getvalue(), report + f"Log: {first / 'assignment_qc.log'}\n",
            )

            shuffled_path = directory / "shuffled.tsv"
            self.write_assignments(shuffled_path, list(reversed(rows)))
            second = directory / "second"
            annotate_assignments(shuffled_path, qc_path, directory, second)
            second_rows = self.read_assignments(second / "assignment_qc.tsv")
            self.assertEqual(
                {row["sample_id"]: row["contam_labels"] for row in second_rows},
                expected,
            )
            self.assertEqual(
                (first / "assignment_qc.log").read_bytes(),
                (second / "assignment_qc.log").read_bytes(),
            )

            self.write_assignments(assignments_path, [rows[0], rows[-1]])
            unaffected = directory / "unaffected"
            stdout = StringIO()
            with redirect_stdout(stdout):
                annotate_assignments(assignments_path, qc_path, directory, unaffected)
            self.assertEqual(
                (unaffected / "assignment_qc.log").read_text(), report_header,
            )
            self.assertEqual(
                stdout.getvalue(),
                report_header + f"Log: {unaffected / 'assignment_qc.log'}\n",
            )
            self.assertEqual(
                [row["contam_labels"] for row in
                 self.read_assignments(unaffected / "assignment_qc.tsv")],
                ["contam_source", ""],
            )

    @staticmethod
    def write_assignments(path: Path, assignments: list[tuple[str, str, int]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", "target", "asv_id", "read_count"])
            writer.writerows(
                (sample, "COI", asv_id, count)
                for asv_id, sample, count in assignments
            )

    @staticmethod
    def read_assignments(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))


if __name__ == "__main__":
    unittest.main()
