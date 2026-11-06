"""Read preprocessing orchestration without external tools."""

import gzip
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from typing import BinaryIO

from scripts.preprocess_reads import preprocess_reads


class PreprocessReadsTests(unittest.TestCase):
    def test_preprocesses_each_nonempty_sample_once_and_records_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            demux = root / "demux"
            demux.mkdir()
            self.write_fastq(demux / "barcode02_barcode89.fastq.gz", b"raw one")
            self.write_fastq(demux / "barcode08_barcode09.fastq.gz", b"raw two")
            self.write_fastq(demux / "barcode08_barcode64.fastq.gz", b"")
            output = root / "preprocessed"
            config = {
                "threads": 3,
                "demux": {"method": "illumina"},
                "preprocess": {
                    "porechop_abi_options": "--discard_database -v 0",
                    "chopper_options": "--minlength 10",
                },
            }
            samples = [
                self.sample("ANLCx_002_A12", "02", "89"),
                self.sample("BFHK_023_H8", "08", "64"),
                self.sample("BFHK_023_A2", "08", "09"),
                self.sample("ANLCx_001_A1", "09", "01"),
            ]
            fastq = b"@read\nACGT\n+\nIIII\n"
            commands = []

            def run(
                command: list[str],
                *,
                stdout: BinaryIO | None = None,
                **kwargs: object,
            ) -> None:
                commands.append(command)
                if command[0] == "porechop_abi":
                    path = Path(command[command.index("-o") + 1])
                    path.write_bytes(b"trimmed")
                else:
                    assert stdout is not None
                    stdout.write(fastq)

            with (
                patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                patch("scripts.preprocess_reads.subprocess.run", side_effect=run),
            ):
                preprocess_reads(root / "manifest.tsv", root, demux, output)

            self.assertEqual(
                [command[0] for command in commands],
                ["porechop_abi", "chopper", "porechop_abi", "chopper"],
            )
            self.assertIn("--discard_database", commands[0])
            self.assertIn("-v", commands[0])
            self.assertIn("--minlength", commands[1])
            self.assertEqual(commands[0][commands[0].index("-t") + 1], "3")
            self.assertEqual(commands[1][commands[1].index("-t") + 1], "3")
            self.assertEqual(
                gzip.open(output / "barcode02_barcode89.fastq.gz", "rb").read(),
                fastq,
            )
            self.assertEqual(
                gzip.open(output / "barcode08_barcode09.fastq.gz", "rb").read(),
                fastq,
            )
            self.assertEqual(
                (output / "excluded_samples.tsv").read_text(),
                "sample_id\treason\n"
                "BFHK_023_H8\tempty_fastq\n"
                "ANLCx_001_A1\tno_demux_file\n",
            )

    def test_tool_failure_is_not_swallowed_and_existing_output_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            demux = root / "demux"
            demux.mkdir()
            self.write_fastq(demux / "barcode02_barcode89.fastq.gz", b"raw")
            output = root / "preprocessed"
            config = {"threads": 1, "demux": {"method": "illumina"}}
            samples = [self.sample("ANLCx_002_A12", "02", "89")]
            failure = subprocess.CalledProcessError(1, ["porechop_abi"])
            with (
                patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                patch("scripts.preprocess_reads.subprocess.run", side_effect=failure),
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    preprocess_reads(root / "manifest.tsv", root, demux, output)
            self.assertTrue(output.is_dir())

            previous = output / "keep.txt"
            previous.write_text("existing")
            with self.assertRaisesRegex(FileExistsError, "already exists"):
                with (
                    patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                    patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                ):
                    preprocess_reads(root / "manifest.tsv", root, demux, output)
            self.assertEqual(previous.read_text(), "existing")

    @staticmethod
    def sample(sample_id: str, front_index: str, rear_index: str) -> dict[str, object]:
        return {
            "sample_id": sample_id,
            "barcode_indices": {
                "front_barcode": front_index,
                "rear_barcode": rear_index,
            },
        }

    @staticmethod
    def write_fastq(path: Path, contents: bytes) -> None:
        with gzip.open(path, "wb") as handle:
            handle.write(contents)


if __name__ == "__main__":
    unittest.main()
