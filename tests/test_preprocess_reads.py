"""Read preprocessing orchestration without external tools."""

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from scripts.preprocess_reads import preprocess_reads, run_chopper


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

            def run(command: list[str], **kwargs: object) -> None:
                commands.append(command)
                path = Path(command[command.index("-o") + 1])
                scratch = Path(command[command.index("--temp_dir") + 1])
                self.assertEqual(scratch, path.parent / "tmp")
                self.assertEqual(path.parent.parent, output.resolve() / "work")
                scratch.mkdir()
                (scratch / "tool-work.txt").write_text("scratch")
                self.write_fastq(path, b"trimmed")

            def chopper(input_path, output_path, options, threads):
                self.write_fastq(output_path, fastq)

            with (
                patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                patch("scripts.preprocess_reads.subprocess.run", side_effect=run),
                patch("scripts.preprocess_reads.run_chopper", side_effect=chopper) as pipe,
            ):
                preprocess_reads(root / "manifest.tsv", root, demux, output)

            self.assertEqual(
                [command[0] for command in commands],
                ["porechop_abi", "porechop_abi"],
            )
            self.assertIn("--discard_database", commands[0])
            self.assertIn("-v", commands[0])
            self.assertEqual(commands[0][commands[0].index("-t") + 1], "3")
            self.assertEqual(pipe.call_count, 2)
            self.assertEqual(pipe.call_args.args[2:], (["--minlength", "10"], "3"))
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
            self.assertFalse((output / "work").exists())

    def test_keep_intermediates_retains_stage_work(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            demux = root / "demux"
            demux.mkdir()
            self.write_fastq(demux / "barcode02_barcode89.fastq.gz", b"raw")
            output = root / "preprocessed"
            config = {"threads": 1, "demux": {"method": "illumina"}}
            samples = [self.sample("ANLCx_002_A12", "02", "89")]

            def run(command: list[str], **kwargs: object) -> None:
                path = Path(command[command.index("-o") + 1])
                self.write_fastq(path, b"trimmed")
                scratch = Path(command[command.index("--temp_dir") + 1])
                scratch.mkdir()
                (scratch / "tool-work.txt").write_text("scratch")

            with (
                patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                patch("scripts.preprocess_reads.subprocess.run", side_effect=run),
                patch("scripts.preprocess_reads.run_chopper") as chopper,
            ):
                preprocess_reads(root / "manifest.tsv", root, demux, output, True)
                chopper.side_effect = subprocess.CalledProcessError(7, ["chopper"])
                with self.assertRaises(subprocess.CalledProcessError):
                    preprocess_reads(
                        root / "manifest.tsv", root, demux, root / "failed",
                    )
            self.assertTrue(
                (root / "failed/work/barcode02_barcode89/porechop.fastq.gz").is_file(),
            )
            sample_work = output / "work/barcode02_barcode89"
            with gzip.open(sample_work / "porechop.fastq.gz", "rb") as handle:
                self.assertEqual(handle.read(), b"trimmed")
            self.assertEqual(
                (sample_work / "tmp/tool-work.txt").read_text(), "scratch",
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
            def run(command: list[str], **kwargs: object) -> None:
                scratch = Path(command[command.index("--temp_dir") + 1])
                scratch.mkdir()
                (scratch / "failed-tool.txt").write_text("diagnostics")
                raise failure

            with (
                patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                patch("scripts.preprocess_reads.subprocess.run", side_effect=run),
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    preprocess_reads(root / "manifest.tsv", root, demux, output)
            self.assertTrue(output.is_dir())
            self.assertEqual(
                (output / "work/barcode02_barcode89/tmp/failed-tool.txt").read_text(),
                "diagnostics",
            )

            previous = output / "keep.txt"
            previous.write_text("existing")
            with self.assertRaisesRegex(FileExistsError, "already exists"):
                with (
                    patch("scripts.preprocess_reads.load_config", return_value=(config, {})),
                    patch("scripts.preprocess_reads.read_manifest", return_value=samples),
                ):
                    preprocess_reads(root / "manifest.tsv", root, demux, output)
            self.assertEqual(previous.read_text(), "existing")

    def test_compression_pipeline_checks_all_processes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq = b"@read\nACGT\n+\nIIII\n"
            source = root / "input reads.fastq.gz"
            self.write_fastq(source, fastq)
            chopper = root / "chopper"
            chopper.write_text(
                f"#!{sys.executable}\n"
                "import sys\n"
                "sys.stdout.buffer.write(sys.stdin.buffer.read())\n"
            )
            chopper.chmod(0o755)
            with (
                patch.dict(os.environ, {"PATH": f"{root}{os.pathsep}{os.environ['PATH']}"}),
            ):
                output = root / "output reads.fastq.gz"
                run_chopper(source, output, [], "1")
                with gzip.open(output, "rb") as handle:
                    self.assertEqual(handle.read(), fastq)
                chopper.write_text(f"#!{sys.executable}\nimport sys\nsys.exit(7)\n")
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    run_chopper(source, output, [], "1")
                self.assertEqual(raised.exception.returncode, 7)
                self.assertEqual(raised.exception.cmd[0], "chopper")
                chopper.write_text(
                    f"#!{sys.executable}\nimport sys\n"
                    "sys.stdout.buffer.write(sys.stdin.buffer.read())\n"
                )
                source.write_bytes(b"not gzip")
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    run_chopper(source, output, [], "1")
                self.assertEqual(raised.exception.cmd[0], "zcat")
                self.write_fastq(source, fastq)
                gzip_command = root / "gzip"
                gzip_command.write_text(f"#!{sys.executable}\nimport sys\nsys.exit(8)\n")
                gzip_command.chmod(0o755)
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    run_chopper(source, output, [], "1")
                self.assertEqual(raised.exception.returncode, 8)
                self.assertEqual(raised.exception.cmd, ["gzip"])

    def test_pipeline_startup_failure_reaps_started_children(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            real_popen = subprocess.Popen
            processes = []

            missing_command = "chopper"

            def popen(command, **kwargs):
                if command[0] == missing_command:
                    raise FileNotFoundError(missing_command)
                process = real_popen(
                    [sys.executable, "-c", "import time; time.sleep(60)"], **kwargs,
                )
                processes.append(process)
                return process

            for missing_command, expected_count in (("chopper", 1), ("gzip", 2)):
                processes.clear()
                with patch("scripts.preprocess_reads.subprocess.Popen", side_effect=popen):
                    with self.assertRaises(FileNotFoundError):
                        run_chopper(root / "input.fastq.gz", root / "output.fastq.gz", [], "1")
                self.assertEqual(len(processes), expected_count)
                for process in processes:
                    self.assertIsNotNone(process.returncode)
                    self.assertTrue(process.stdout.closed)

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
