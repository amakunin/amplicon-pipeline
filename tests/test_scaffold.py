"""Smoke tests for the scaffold, not for analysis behaviour."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/small"
sys.path.insert(0, str(ROOT))

import run_pipeline
from run_pipeline import demux_command, load_config, read_manifest


class ScaffoldTests(unittest.TestCase):
    def test_load_fixture_config(self) -> None:
        config, targets = load_config(FIXTURE / "config")
        self.assertEqual(config["demux"]["method"], "hybrid")
        self.assertEqual(
            config["demux"]["reverse_barcode_mapping"],
            ROOT / "assets/barcodes/Sanger_custom_barcode_R.tsv",
        )
        self.assertEqual(config["demux"]["barcode_dir"], ROOT / "assets/barcodes")
        self.assertNotIn("seqtk_exec", config["demux"])
        self.assertEqual(
            config["targets"]["COI"], FIXTURE / "config/targets/COI.yaml"
        )
        self.assertEqual(targets["COI"]["asv_calling"], "savont")

    def test_load_config_resolves_relative_executable_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config_dir = Path(temporary)
            targets_dir = config_dir / "targets"
            targets_dir.mkdir()
            (config_dir / "pipeline.yaml").write_text(
                "targets:\n"
                "  COI: targets/COI.yaml\n"
                "demux:\n"
                "  method: hybrid\n"
                "  dorado_exec: bin/dorado\n"
                "  reverse_barcode_mapping: null\n"
            )
            (targets_dir / "COI.yaml").write_text("asv_calling: savont\n")

            config, _ = load_config(config_dir)

            self.assertEqual(
                config["demux"]["dorado_exec"],
                (config_dir / "bin/dorado").resolve(),
            )
            self.assertIsNone(config["demux"]["reverse_barcode_mapping"])

    def test_read_fixture_manifest(self) -> None:
        config, _ = load_config(FIXTURE / "config")

        samples = read_manifest(FIXTURE / "manifest.tsv", config)

        self.assertEqual(len(samples), 5)
        self.assertEqual(
            samples[0]["barcode_indices"],
            {"ont_barcode": "10", "rear_barcode": "89"},
        )
        self.assertEqual(samples[0]["well_id"], "A12")

    def test_read_manifest_reports_unknown_barcode_with_sample(self) -> None:
        config, _ = load_config(FIXTURE / "config")
        manifest = (FIXTURE / "manifest.tsv").read_text().replace(
            "AACGTGAT", "UNKNOWN"
        )
        with tempfile.TemporaryDirectory() as temporary:
            manifest_path = Path(temporary) / "manifest.tsv"
            manifest_path.write_text(manifest)

            with self.assertRaisesRegex(
                ValueError, "Unknown rear_barcode .*ANLCx_001_A1"
            ):
                read_manifest(manifest_path, config)

    def test_help(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_pipeline.py"), "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for argument in (
            "run_fastq_gz",
            "manifest_tsv",
            "config_dir",
            "output_dir",
            "--profile",
            "--dry-run",
        ):
            self.assertIn(argument, result.stdout)
        self.assertIn("Only demultiplexing runs", result.stdout)

    def test_missing_arguments(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_pipeline.py")],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("required", result.stderr)

    def test_demux_command(self) -> None:
        config, _ = load_config(FIXTURE / "config")
        command = demux_command(config, Path("run.fastq.gz"), Path("out"))
        self.assertEqual(command[:2], ["bash", str(ROOT / "scripts/demux_hybrid.sh")])
        self.assertEqual(
            command[2:],
            [
                "--input", str(Path("run.fastq.gz").resolve()),
                "--output-dir", str(Path("out").resolve() / "demux"),
                "--barcode-dir", str(ROOT / "assets/barcodes"),
                "--threads", "4",
                "--dorado", "/software/sciops/external/ont-dorado/0.8.2/bin/dorado",
            ],
        )

        config["demux"]["method"] = "illumina"
        config["demux"]["dorado_exec"] = None
        command = demux_command(config, Path("run.fastq.gz"), Path("out"))
        self.assertTrue(command[1].endswith("scripts/demux_illumina.sh"))
        self.assertNotIn("--dorado", command)

    def test_execution_runs_demux_then_stops(self) -> None:
        argv = [
            "run_pipeline.py",
            "run.fastq.gz",
            str(FIXTURE / "manifest.tsv"),
            str(FIXTURE / "config"),
            "out",
        ]
        with mock.patch.object(sys, "argv", argv), mock.patch.object(
            run_pipeline.subprocess, "run"
        ) as run:
            with self.assertRaises(SystemExit) as raised:
                run_pipeline.main()
        self.assertEqual(raised.exception.code, 1)
        run.assert_called_once()
        self.assertIn("demux_hybrid.sh", run.call_args.args[0][1])

    def test_dry_run_prints_plan_without_creating_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output_dir = directory / "output"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "run_pipeline.py"),
                    "run.fastq.gz",
                    str(FIXTURE / "manifest.tsv"),
                    str(FIXTURE / "config"),
                    str(output_dir),
                    "--dry-run",
                ],
                cwd=directory,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            plan, demux_line = result.stdout.rsplit("Demux command: ", 1)
            self.assertIn("demux_hybrid.sh", demux_line)
            self.assertEqual(
                plan,
                "Method: hybrid\n"
                "Targets: COI\n"
                "Samples:\n"
                "  ANLCx_002_A12: ont_barcode=10, rear_barcode=89\n"
                "  BFHK_023_H8: ont_barcode=16, rear_barcode=64\n"
                "  BFHK_023_A2: ont_barcode=16, rear_barcode=09\n"
                "  ANLCx_001_A1: ont_barcode=09, rear_barcode=01\n"
                "  ANLCx_001_B1: ont_barcode=09, rear_barcode=02\n",
            )
            self.assertEqual(result.stderr, "")
            self.assertEqual(list(directory.iterdir()), [])

    def test_stage_placeholders_fail_explicitly(self) -> None:
        commands = [
            [sys.executable, str(ROOT / "scripts" / name)]
            for name in (
                "summarise_asvs.py",
                "summarise_samples.py",
                "plot_summary.py",
            )
        ]
        with tempfile.TemporaryDirectory() as temporary:
            for command in commands:
                with self.subTest(script=command[-1]):
                    result = subprocess.run(
                        command,
                        cwd=temporary,
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    self.assertEqual(result.returncode, 1)
                    self.assertIn("not implemented", result.stderr)
                    self.assertEqual(list(Path(temporary).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
