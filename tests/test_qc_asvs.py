"""Sequence QC evidence without external tools or assignment filtering."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml

from scripts.qc_asvs import annotate_asvs, check_frames
from scripts.summarise_asvs import summarise_asvs


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/small"


class SequenceQcTests(unittest.TestCase):
    def test_six_frames_and_codon_handling(self) -> None:
        # Alternate frames contain TAA or TAG; +0 contains CTA, ACT and AGT.
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
        # Completing the trailing TA creates a stop, including at the terminus.
        self.assertEqual(check_frames(sequence + "A"), ([], []))
        self.assertEqual(check_frames(sequence + "N"), ([], ["+0"]))
        self.assertEqual(check_frames("TAANNN"), (
            [], ["+1", "+2", "-0", "-1", "-2"],
        ))

    def test_coding_statuses(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "asvs.tsv"
            source.write_text(
                "target\tasv_id\tsequence\n"
                "COI\tone\tCTAACTAGTTA\n"
                "COI\tnone\tCTAACTAGTTAA\n"
                "COI\tuncertain\tCTAACTAGTTAN\n"
            )
            annotate_asvs(source, FIXTURE / "config", directory)
            with (directory / "asv_qc.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(r["coding_status"], r["stop_free_frames"], r["uncertain_frames"])
                 for r in rows],
                [("one_stop_free_frame", "+0", ""),
                 ("no_stop_free_frame", "", ""), ("uncertain", "", "+0")],
            )

    def test_length_boundaries_and_common_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            settings = yaml.safe_load(
                (FIXTURE / "config/pipeline.yaml").read_text()
            )
            settings["summary"]["common"] = {"min_length": 1}
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(settings))
            source = directory / "asvs.tsv"
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
            with source.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["target", "asv_id", "sequence"])
                for target, length, _ in cases:
                    writer.writerow([target, str(length), "C" * length])
            annotate_asvs(source, directory, directory / "summary")
            with (directory / "summary/asv_qc.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(r["target"], int(r["length"]), r["length_status"]) for r in rows],
                cases,
            )
            self.assertEqual(rows[0]["coding_status"], "multiple_stop_free_frames")
            self.assertEqual(rows[-1]["coding_status"], "not_checked")
            settings["summary"]["targets"]["COI"]["genetic_code"] = 1
            (directory / "pipeline.yaml").write_text(yaml.safe_dump(settings))
            before = (directory / "summary/asv_qc.tsv").read_bytes()
            with self.assertRaisesRegex(ValueError, "Only genetic code 5"):
                annotate_asvs(source, directory, directory / "summary")
            self.assertEqual((directory / "summary/asv_qc.tsv").read_bytes(), before)

    def test_fixture_annotations_leave_raw_tables_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            summarise_asvs(
                {target: FIXTURE / "ampliseq" / target
                 for target in ("COI", "CytB", "ITS2", "P3")},
                output,
            )
            before = {name: (output / name).read_bytes()
                      for name in ("asvs.tsv", "assignments.tsv")}
            subprocess.run([
                sys.executable, str(ROOT / "scripts/qc_asvs.py"),
                "--asvs", str(output / "asvs.tsv"),
                "--config-dir", str(FIXTURE / "config"),
                "--output-dir", str(output),
            ], check=True)
            with (output / "asv_qc.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                self.assertEqual(reader.fieldnames, [
                    "target", "asv_id", "length", "length_status",
                    "stop_free_frames", "uncertain_frames", "coding_status",
                ])
            self.assertEqual(
                [(r["target"], r["length_status"], r["coding_status"]) for r in rows],
                [("CytB", "pass", "multiple_stop_free_frames"),
                 ("ITS2", "not_checked", "not_checked"),
                 ("P3", "too_short", "not_checked"), ("P3", "pass", "not_checked")],
            )
            self.assertEqual(
                [(r["asv_id"], r["length"], r["length_status"]) for r in rows[2:]],
                [("0a4b74ea4eddec52", "371", "too_short"),
                 ("83ea72fd937a3d48", "373", "pass")],
            )
            self.assertEqual(rows[0]["stop_free_frames"], "+0,+1")
            self.assertEqual(rows[0]["uncertain_frames"], "")
            self.assertTrue(all(not r["stop_free_frames"] and not r["uncertain_frames"]
                                for r in rows[1:]))
            self.assertEqual(
                {name: (output / name).read_bytes() for name in before}, before
            )


if __name__ == "__main__":
    unittest.main()
