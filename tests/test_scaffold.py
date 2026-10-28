"""Smoke tests for the scaffold, not for analysis behaviour."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/coi_small"
sys.path.insert(0, str(ROOT))

from run_pipeline import load_config


class ScaffoldTests(unittest.TestCase):
    def test_load_fixture_config(self) -> None:
        config, targets = load_config(FIXTURE / "config")
        self.assertEqual(config["demux"]["method"], "hybrid")
        self.assertEqual(
            config["demux"]["reverse_barcode_mapping"],
            ROOT / "assets/barcodes/Sanger_custom_barcode_R.tsv",
        )
        self.assertEqual(
            config["targets"]["COI"], FIXTURE / "config/targets/COI.yaml"
        )
        self.assertEqual(targets["COI"]["asv_calling"], "savont")

    def test_help(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_pipeline.py"), "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for argument in (
            "run_fastq_gz", "manifest_tsv", "config_dir", "output_dir", "--profile"
        ):
            self.assertIn(argument, result.stdout)
        self.assertIn("not implemented", result.stdout)

    def test_missing_arguments(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_pipeline.py")],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("required", result.stderr)

    def test_execution_fails_without_creating_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for profile in ("local", "lsf"):
                with self.subTest(profile=profile):
                    result = subprocess.run(
                        [
                            sys.executable,
                            str(ROOT / "run_pipeline.py"),
                            "run.fastq.gz",
                            "manifest.tsv",
                            str(FIXTURE / "config"),
                            str(directory / "output"),
                            "--profile",
                            profile,
                        ],
                        cwd=directory,
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    self.assertEqual(result.returncode, 1)
                    self.assertIn("not implemented", result.stderr)
                    self.assertEqual(list(directory.iterdir()), [])

    def test_stage_placeholders_fail_explicitly(self) -> None:
        commands = [
            [sys.executable, str(ROOT / "scripts" / name)]
            for name in (
                "prepare_inputs.py",
                "summarise_asvs.py",
                "summarise_samples.py",
                "plot_summary.py",
            )
        ]
        commands.append(["bash", str(ROOT / "scripts/demux_hybrid.sh")])
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
