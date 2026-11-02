"""Driver configuration, command/restart and CLI tests without external tools."""

import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/small"
sys.path.insert(0, str(ROOT))

import run_pipeline
from run_pipeline import (
    ampliseq_command,
    check_nextflow,
    demux_command,
    load_config,
    read_manifest,
    run_ampliseq,
    run_demux,
    write_resource_config,
)


class ConfigManifestTests(unittest.TestCase):
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
            config["targets"]["CytB"], FIXTURE / "config/targets/CytB.yaml"
        )
        self.assertEqual(targets["CytB"]["asv_calling"], "savont")

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
            (targets_dir / "COI.yaml").write_text(
                "asv_calling: savont\n"
                "sintax_ref: refs/legacy.fasta\n"
                "sintax_ref_tax_custom: ../culicidae_ref/CytB.fasta\n"
            )

            config, targets = load_config(config_dir)

            self.assertEqual(
                config["demux"]["dorado_exec"],
                (config_dir / "bin/dorado").resolve(),
            )
            self.assertIsNone(config["demux"]["reverse_barcode_mapping"])
            self.assertEqual(
                targets["COI"]["sintax_ref"],
                str((config_dir / "refs/legacy.fasta").resolve()),
            )
            self.assertEqual(
                targets["COI"]["sintax_ref_tax_custom"],
                str((config_dir.parent / "culicidae_ref/CytB.fasta").resolve()),
            )

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


class CommandRestartTests(unittest.TestCase):
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

    def test_ampliseq_command_and_resource_config(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            target_dir = Path(temporary)
            write_resource_config(target_dir, 4)
            self.assertEqual(
                (target_dir / "resources.config").read_text(),
                "executor.cpus = 4\nprocess.resourceLimits = [cpus: 4]\n",
            )
            self.assertEqual(
                ampliseq_command(FIXTURE / "config", target_dir, "local"),
                [
                    "nextflow", "run", "nf-core/ampliseq", "-r", "3.0.0",
                    "-resume",
                    "-profile", "local",
                    "-c", str(FIXTURE / "config/ampliseq.config"),
                    "-c", str(target_dir / "resources.config"),
                    "-params-file", str(target_dir / "params.yaml"),
                ],
            )

    def test_ampliseq_completion_record_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            demux_dir = output_dir / "demux"
            demux_dir.mkdir(parents=True)
            with gzip.open(
                demux_dir / "ont_barcode10_barcode89.fastq.gz", "wb"
            ) as handle:
                handle.write(b"@r\nACGT\n+\nIIII\n")
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            target_dir = output_dir / "ampliseq/CytB"
            reference = Path(temporary) / "reference.fasta"
            reference.write_text(">ref\nACGT\n")
            args = SimpleNamespace(
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                profile="local",
                run_fastq_gz=run_fastq,
            )
            nextflow_commands = []
            fail_nextflow = [True]

            def fake_run(command, **kwargs):
                if "prepare_inputs.py" in command[1]:
                    target_dir.mkdir(parents=True, exist_ok=True)
                    reads = demux_dir / "ont_barcode10_barcode89.fastq.gz"
                    (target_dir / "samplesheet.tsv").write_text(
                        f"sample\tfastq_1\nsample\t{reads}\n"
                    )
                    (target_dir / "excluded_samples.tsv").write_text(
                        "sample_id\treason\n"
                    )
                    (target_dir / "params.yaml").write_text(
                        run_pipeline.yaml.safe_dump(
                            {
                                "input": str(target_dir / "samplesheet.tsv"),
                                "outdir": str(target_dir / "results"),
                                "sintax_ref_tax_custom": str(reference),
                            }
                        )
                    )
                    return subprocess.CompletedProcess(command, 0)
                if command[0] != "nextflow":
                    raise AssertionError(f"Unexpected command: {command}")
                nextflow_commands.append(command)
                if fail_nextflow[0]:
                    raise subprocess.CalledProcessError(1, command)
                results_dir = Path(kwargs["cwd"]) / "results"
                results_dir.mkdir(exist_ok=True)
                (results_dir / "overall_summary.tsv").write_text("result\n")
                return subprocess.CompletedProcess(command, 0)

            record_path = output_dir / "ampliseq/CytB/.completion.json"
            with mock.patch.object(
                run_pipeline.subprocess, "run", side_effect=fake_run
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    run_ampliseq(config, args)
                self.assertFalse(record_path.exists())

                fail_nextflow[0] = False
                run_ampliseq(config, args)
                self.assertTrue(record_path.exists())
                run_ampliseq(config, args)
                self.assertEqual(len(nextflow_commands), 2)

                (target_dir / "results/overall_summary.tsv").unlink()
                run_ampliseq(config, args)
                self.assertEqual(len(nextflow_commands), 3)

                args.profile = "lsf"
                run_ampliseq(config, args)
                self.assertEqual(len(nextflow_commands), 4)

                config["threads"] += 1
                run_ampliseq(config, args)
                self.assertEqual(len(nextflow_commands), 5)

                with gzip.open(
                    demux_dir / "ont_barcode10_barcode89.fastq.gz", "wb"
                ) as handle:
                    handle.write(b"@r2\nTGCA\n+\nIIII\n")
                run_ampliseq(config, args)
                self.assertEqual(len(nextflow_commands), 6)

                reference.write_text(">ref\nTGCA\n")
                run_ampliseq(config, args)
                self.assertEqual(len(nextflow_commands), 7)

            self.assertTrue(all("-resume" in command for command in nextflow_commands))

    def test_demux_completion_record_and_partial_output_guard(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            output_dir.mkdir()
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                run_fastq_gz=run_fastq,
            )
            calls = []

            def fake_run(command, **kwargs):
                calls.append(command)
                demux_dir = Path(command[command.index("--output-dir") + 1])
                demux_dir.mkdir()
                (demux_dir / "reads.fastq.gz").write_bytes(b"reads")
                (demux_dir / "demux_stats.txt").write_text("reads\t1\n")
                return subprocess.CompletedProcess(command, 0)

            with mock.patch.object(
                run_pipeline.subprocess, "run", side_effect=fake_run
            ):
                run_demux(config, args)
                record_path = output_dir / ".demux-completion.json"
                self.assertTrue(record_path.exists())
                run_demux(config, args)
                self.assertEqual(len(calls), 1)

                (output_dir / "demux/reads.fastq.gz").unlink()
                with self.assertRaisesRegex(SystemExit, "output files are missing"):
                    run_demux(config, args)
                self.assertEqual(len(calls), 1)
                self.assertTrue((output_dir / "demux/demux_stats.txt").exists())

                (output_dir / "demux/reads.fastq.gz").write_bytes(b"reads")
                run_fastq.write_bytes(b"changed raw run")
                with self.assertRaisesRegex(SystemExit, "inputs or configuration"):
                    run_demux(config, args)
                self.assertEqual(len(calls), 1)

                record_path.unlink()
                with self.assertRaisesRegex(
                    SystemExit, "without a completion record"
                ):
                    run_demux(config, args)
                self.assertEqual(len(calls), 1)

    def test_demux_failure_does_not_write_completion_record(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            output_dir.mkdir()
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                run_fastq_gz=run_fastq,
            )

            def fail_after_partial_output(command, **kwargs):
                demux_dir = Path(command[command.index("--output-dir") + 1])
                (demux_dir / ".work").mkdir(parents=True)
                raise subprocess.CalledProcessError(1, command)

            with mock.patch.object(
                run_pipeline.subprocess,
                "run",
                side_effect=fail_after_partial_output,
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    run_demux(config, args)

            self.assertFalse((output_dir / ".demux-completion.json").exists())
            with self.assertRaisesRegex(SystemExit, "without a completion record"):
                run_demux(config, args)

    def test_check_nextflow_requires_version_26(self) -> None:
        def completed(version: str) -> subprocess.CompletedProcess:
            return subprocess.CompletedProcess([], 0, stdout=f"\n  version {version} build 1\n")

        with mock.patch.object(
            run_pipeline.subprocess, "run", return_value=completed("26.04.0")
        ):
            check_nextflow()
        with mock.patch.object(
            run_pipeline.subprocess, "run", return_value=completed("25.10.4")
        ):
            with self.assertRaisesRegex(SystemExit, "26 or higher"):
                check_nextflow()
        with mock.patch.object(
            run_pipeline.subprocess, "run", side_effect=FileNotFoundError
        ):
            with self.assertRaisesRegex(SystemExit, "not found"):
                check_nextflow()


class CliTests(unittest.TestCase):
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
        self.assertIn("later", result.stdout)

    def test_missing_arguments(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_pipeline.py")],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("required", result.stderr)

    def test_execution_runs_demux_then_ampliseq_then_stops(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            argv = [
                "run_pipeline.py",
                str(Path(temporary) / "run.fastq.gz"),
                str(FIXTURE / "manifest.tsv"),
                str(FIXTURE / "config"),
                str(Path(temporary) / "out"),
            ]
            Path(argv[1]).write_text("raw run")
            target_dir = Path(temporary) / "out/ampliseq/CytB"

            def fake_run(command, **kwargs):
                if "demux_" in command[1]:
                    demux_dir = Path(
                        command[command.index("--output-dir") + 1]
                    )
                    demux_dir.mkdir(parents=True)
                    (demux_dir / "demux_stats.txt").write_text("reads\t1\n")
                    return subprocess.CompletedProcess(command, 0)
                if "prepare_inputs.py" in command[1]:
                    target_dir.mkdir(parents=True, exist_ok=True)
                    reads = Path(temporary) / "input.fastq.gz"
                    reads.write_bytes(b"reads")
                    (target_dir / "samplesheet.tsv").write_text(
                        f"sample\tfastq_1\nsample\t{reads}\n"
                    )
                    (target_dir / "excluded_samples.tsv").write_text(
                        "sample_id\treason\n"
                    )
                    (target_dir / "params.yaml").write_text("input: samplesheet.tsv\n")
                    return subprocess.CompletedProcess(command, 0)
                if command[0] == "nextflow":
                    results_dir = target_dir / "results"
                    results_dir.mkdir(parents=True, exist_ok=True)
                    (results_dir / "overall_summary.tsv").write_text("result\n")
                    return subprocess.CompletedProcess(command, 0)

            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                run_pipeline, "check_nextflow"
            ), mock.patch.object(
                run_pipeline.subprocess, "run", side_effect=fake_run
            ) as run:
                with self.assertRaises(SystemExit) as raised:
                    run_pipeline.main()

            self.assertEqual(raised.exception.code, 1)
            commands = [call.args[0] for call in run.call_args_list]
            self.assertIn("demux_hybrid.sh", commands[0][1])
            self.assertIn("prepare_inputs.py", commands[1][1])
            self.assertEqual(commands[2][:3], ["nextflow", "run", "nf-core/ampliseq"])
            self.assertEqual(run.call_args_list[2].kwargs["cwd"], target_dir)

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
            plan, commands = result.stdout.split("Demux command: ", 1)
            self.assertIn("demux_hybrid.sh", commands)
            self.assertIn("Ampliseq command (CytB): nextflow run", commands)
            self.assertEqual(
                plan,
                "Method: hybrid\n"
                "Targets: CytB\n"
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
