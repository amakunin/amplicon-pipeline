"""Assignment QC flags and unique-ASV totals without external tools."""

import csv
from pathlib import Path
import tempfile
import unittest

import yaml

from scripts.summarise_asvs import annotate_assignments, summarise_asvs


FIXTURE = Path(__file__).resolve().parents[1] / "tests/fixtures/small"


class AssignmentQcTests(unittest.TestCase):
    def test_boundaries_overrides_and_order_independence(self) -> None:
        settings = {
            "COI": {},
            "CytB": {"min_contam_total_reads": 5000, "max_contam_prop": 0.02, "max_locov": 50},
        }
        asvs = [
            {"target": "COI", "asv_id": "shared", "length_status": "too_short"},
            {"target": "COI", "asv_id": "boundary"},
            {"target": "CytB", "asv_id": "shared"},
        ]
        rows = [
            ("COI", "shared", "A", 9000), ("COI", "shared", "B", 9000),
            ("COI", "shared", "C", 1601), ("COI", "shared", "D", 200),
            ("COI", "shared", "E", 100), ("COI", "shared", "F", 99),
            ("COI", "boundary", "G", 9901), ("COI", "boundary", "H", 99),
            ("CytB", "shared", "I", 9901), ("CytB", "shared", "J", 99),
        ]
        assignments = [
            {"target": target, "asv_id": asv_id, "sample_id": sample, "read_count": count}
            for target, asv_id, sample, count in rows
        ]
        report = annotate_assignments(asvs, assignments, settings)
        expected = {
            "A": ("false", "false"), "B": ("false", "false"),
            "C": ("false", "false"), "D": ("false", "false"),
            "E": ("true", "false"), "F": ("true", "true"),
            "G": ("false", "false"), "H": ("false", "true"),
            "I": ("false", "false"), "J": ("true", "false"),
        }
        self.assertEqual(
            {r["sample_id"]: (r["contam_possible"], r["locov"]) for r in assignments},
            expected,
        )
        self.assertEqual(asvs[0]["total_reads"], 20000)
        self.assertEqual(asvs[0]["contam_target_cutoff"], "200.0")
        self.assertEqual(asvs[1]["total_reads"], 10000)
        self.assertEqual(report, [
            "Contamination: COI, shared; total reads: 20000; cutoff: 200.0; affected samples: 2",
            "Contamination: CytB, shared; total reads: 10000; cutoff: 200.0; affected samples: 1",
        ])
        self.assertEqual(
            annotate_assignments(list(reversed(asvs)), list(reversed(assignments)), settings),
            report,
        )
        self.assertEqual(
            {r["sample_id"]: (r["contam_possible"], r["locov"]) for r in assignments},
            expected,
        )
        self.assertTrue(all(set(r) == {
            "target", "asv_id", "sample_id", "read_count", "contam_possible", "locov",
        } for r in assignments))

    def test_fixture_canonical_flags_and_totals(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            summarise_asvs(
                {target: FIXTURE / "ampliseq" / target
                 for target in ("COI", "CytB", "ITS2", "P3")},
                output, FIXTURE / "config",
            )
            with (output / "assignments.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(int(r["read_count"]), r["contam_possible"], r["locov"]) for r in rows],
                [(150, "false", "false"), (12, "false", "true"),
                 (9, "false", "true"), (11, "false", "true"),
                 (10, "false", "true"), (1, "false", "true"), (24, "false", "true")],
            )
            with (output / "asvs.tsv").open(newline="") as handle:
                asvs = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([int(r["total_reads"]) for r in asvs], [162, 20, 11, 24])
            self.assertTrue(all("total_reads" not in r and "length_status" not in r for r in rows))

    def test_obsolete_settings_do_not_replace_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = yaml.safe_load((FIXTURE / "config/pipeline.yaml").read_text())
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(config))
            output = directory / "summary"
            summarise_asvs({"CytB": FIXTURE / "ampliseq/CytB"}, output, directory)
            before = {p.name: p.read_bytes() for p in output.iterdir()}
            config["summary"]["common"]["min_contam_source_reads"] = 500
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(config))
            with self.assertRaisesRegex(ValueError, "min_contam_source_reads"):
                summarise_asvs({"CytB": FIXTURE / "ampliseq/CytB"}, output, directory)
            self.assertEqual({p.name: p.read_bytes() for p in output.iterdir()}, before)


if __name__ == "__main__":
    unittest.main()
