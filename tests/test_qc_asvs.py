"""Sequence QC in canonical ASV summaries without external tools."""

import csv
from pathlib import Path
import tempfile
import unittest

import yaml

from scripts.summarise_asvs import annotate_asvs, check_frames, summarise_asvs


FIXTURE = Path(__file__).resolve().parents[1] / "tests/fixtures/small"


class SequenceQcTests(unittest.TestCase):
    def test_six_frames_and_codon_handling(self) -> None:
        sequence = "CTAACTAGTTA"
        for shift in range(3):
            forward = "C" * shift + sequence
            reverse = forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            with self.subTest(shift=shift):
                self.assertEqual(check_frames(forward.lower()), ([f"+{shift}"], []))
                self.assertEqual(check_frames(reverse), ([f"-{shift}"], []))
        self.assertEqual(check_frames("TGATGATGA"), (
            ["+0", "+1", "+2", "-0", "-1", "-2"], [],
        ))
        self.assertEqual(check_frames(sequence + "A"), ([], []))
        self.assertEqual(check_frames(sequence + "N"), ([], ["+0"]))
        self.assertEqual(check_frames("TAANNN"), (
            [], ["+1", "+2", "-0", "-1", "-2"],
        ))

    def test_coding_statuses(self) -> None:
        rows = [
            {"target": "COI", "asv_id": name, "sequence": sequence}
            for name, sequence in (
                ("one", "CTAACTAGTTA"), ("none", "CTAACTAGTTAA"),
                ("uncertain", "CTAACTAGTTAN"),
            )
        ]
        annotate_asvs(rows, {"COI": {"genetic_code": 5}})
        self.assertEqual(
            [(r["coding_status"], r["stop_free_frames"], r["uncertain_frames"])
             for r in rows],
            [("one_stop_free_frame", "+0", ""),
             ("no_stop_free_frame", "", ""), ("uncertain", "", "+0")],
        )

    def test_length_boundaries_and_common_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = yaml.safe_load((FIXTURE / "config/pipeline.yaml").read_text())
            config["summary"]["common"] = {"min_length": 1}
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(config))
            cases = [
                ("COI", 639, "too_short"), ("COI", 641, "unexpected_length"),
                ("COI", 674, "too_long"),
                *[("COI", length, "pass") for length in range(640, 674, 3)],
                ("CytB", 297, "too_short"), ("CytB", 299, "unexpected_length"),
                ("CytB", 317, "too_long"),
                *[("CytB", length, "pass") for length in range(298, 317, 3)],
                ("P3", 371, "too_short"), ("P3", 372, "pass"),
                ("P3", 373, "pass"), ("P3", 380, "pass"), ("P3", 381, "too_long"),
            ]
            rows = [
                {"target": target, "asv_id": str(length), "sequence": "C" * length}
                for target, length, _ in cases
            ]
            settings = {
                target: {**config["summary"]["common"], **overrides}
                for target, overrides in config["summary"]["targets"].items()
            }
            annotate_asvs(rows, settings)
            self.assertEqual(
                [(r["target"], r["length"], r["length_status"]) for r in rows], cases,
            )
            self.assertEqual(rows[0]["coding_status"], "multiple_stop_free_frames")
            self.assertEqual(rows[-1]["coding_status"], "not_checked")
            output = directory / "summary"
            summarise_asvs({"CytB": FIXTURE / "ampliseq/CytB"}, output, directory)
            before = {p.name: p.read_bytes() for p in output.iterdir()}
            config["summary"]["targets"]["CytB"]["genetic_code"] = 1
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(config))
            with self.assertRaisesRegex(ValueError, "Only genetic code 5"):
                summarise_asvs({"CytB": FIXTURE / "ampliseq/CytB"}, output, directory)
            self.assertEqual({p.name: p.read_bytes() for p in output.iterdir()}, before)
            with self.assertRaisesRegex(ValueError, "Conflicting length bounds"):
                annotate_asvs(rows, {**settings, "COI": {"min_length": 700, "max_length": 600}})

    def test_fixture_annotations_in_canonical_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            summarise_asvs(
                {target: FIXTURE / "ampliseq" / target
                 for target in ("COI", "CytB", "ITS2", "P3")},
                output, FIXTURE / "config",
            )
            with (output / "asvs.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(r["target"], r["length_status"], r["coding_status"]) for r in rows],
                [("CytB", "pass", "multiple_stop_free_frames"),
                 ("ITS2", "not_checked", "not_checked"),
                 ("P3", "too_short", "not_checked"), ("P3", "pass", "not_checked")],
            )
            self.assertEqual(rows[0]["stop_free_frames"], "+0,+1")
            self.assertEqual(rows[0]["uncertain_frames"], "")
            self.assertTrue(all(not r["stop_free_frames"] and not r["uncertain_frames"]
                                for r in rows[1:]))
            self.assertEqual(
                sorted(p.name for p in output.iterdir()),
                ["assignments.tsv", "asvs.tsv", "summarise_asvs.log"],
            )


if __name__ == "__main__":
    unittest.main()
