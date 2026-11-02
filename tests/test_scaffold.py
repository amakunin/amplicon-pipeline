"""Driver configuration, command/restart and CLI tests without external tools."""

import gzip
from pathlib import Path
import shlex
import shutil
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
        self.assertEqual(config["demux"]["method"], "illumina")
        self.assertEqual(
            config["demux"]["reverse_barcode_mapping"],
            ROOT / "assets/barcodes/Sanger_custom_barcode_R.tsv",
        )
        self.assertEqual(config["demux"]["barcode_dir"], ROOT / "assets/barcodes")
        self.assertNotIn("seqtk_exec", config["demux"])
        self.assertEqual(config["ampliseq"]["targets"], targets)
        self.assertEqual(list(targets), ["COI", "CytB", "ITS2", "P3"])
        self.assertEqual(targets["CytB"]["asv_calling"], "savont")

    def test_load_config_merges_parameters_and_resolves_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config_dir = Path(temporary)
            (config_dir / "pipeline.yaml").write_text(
                "ampliseq:\n"
                "  common:\n"
                "    asv_calling: savont\n"
                "    skip_fastqc: true\n"
                "    savont_options: '--single-strand'\n"
                "    sintax_ref: refs/legacy.fasta\n"
                "  targets:\n"
                "    COI:\n"
                "      skip_fastqc: false\n"
                "      savont_options: '--min-cluster-size 4'\n"
                "      sintax_ref_tax_custom: ../culicidae_ref/CytB.fasta\n"
                "    CytB: {}\n"
                "summary:\n"
                "  common: {}\n"
                "  targets: {}\n"
                "demux:\n"
                "  method: hybrid\n"
                "  dorado_exec: bin/dorado\n"
                "  reverse_barcode_mapping: null\n"
            )
            config, targets = load_config(config_dir)

            self.assertEqual(targets["COI"]["asv_calling"], "savont")
            self.assertFalse(targets["COI"]["skip_fastqc"])
            self.assertEqual(targets["COI"]["savont_options"], "--min-cluster-size 4")
            self.assertTrue(targets["CytB"]["skip_fastqc"])
            self.assertEqual(targets["CytB"]["savont_options"], "--single-strand")
            self.assertEqual(config["ampliseq"]["common"]["sintax_ref"], "refs/legacy.fasta")
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
            {"front_barcode": "02", "rear_barcode": "89"},
        )
        self.assertEqual(samples[0]["well_id"], "A12")
        config["demux"]["method"] = "hybrid"
        config["demux"]["ont_barcode_mapping"] = (
            ROOT / "assets/barcodes/ONT_nb24_barcode.tsv"
        )
        self.assertEqual(
            read_manifest(FIXTURE / "manifest.tsv", config)[0]["barcode_indices"],
            {"ont_barcode": "10", "rear_barcode": "89"},
        )

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
    def test_stage_configuration_fingerprints_are_isolated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config_dir = directory / "config"
            shutil.copytree(FIXTURE / "config", config_dir)
            barcodes = directory / "barcodes"
            shutil.copytree(ROOT / "assets/barcodes", barcodes)
            pipeline_path = config_dir / "pipeline.yaml"
            settings = run_pipeline.yaml.safe_load(pipeline_path.read_text())
            settings["ampliseq"]["targets"] = {
                name: settings["ampliseq"]["targets"][name] for name in ("COI", "CytB")
            }
            settings["demux"]["barcode_dir"] = str(barcodes)
            for key in ("front_barcode_mapping", "reverse_barcode_mapping"):
                settings["demux"][key] = str(
                    barcodes / Path(settings["demux"][key]).name
                )
            for target, parameters in settings["ampliseq"]["targets"].items():
                reference = directory / f"{target}.fasta"
                reference.write_text(">ref\nACGT\n")
                parameters["sintax_ref_tax_custom"] = str(reference)
            pipeline_path.write_text(run_pipeline.yaml.safe_dump(settings))
            args = SimpleNamespace(
                config_dir=config_dir,
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=directory / "output",
                run_fastq_gz=directory / "run.fastq.gz",
                profile="local",
            )
            args.run_fastq_gz.write_bytes(b"raw run")
            demux_dir = args.output_dir / "demux"
            demux_dir.mkdir(parents=True)

            def fingerprints():
                config, _ = load_config(config_dir)
                targets = {}
                for target in config["ampliseq"]["targets"]:
                    subprocess.run([
                        sys.executable, str(ROOT / "scripts/prepare_inputs.py"),
                        "--manifest", str(args.manifest_tsv),
                        "--config-dir", str(config_dir),
                        "--demux-dir", str(demux_dir),
                        "--target", target,
                        "--output-dir", str(args.output_dir / "ampliseq"),
                    ], check=True)
                    target_dir = args.output_dir / "ampliseq" / target
                    write_resource_config(target_dir, config["threads"])
                    targets[target] = run_pipeline.target_fingerprints(
                        args, config, target, target_dir
                    )
                demux = run_pipeline.demux_fingerprints(
                    config, args, demux_command(config, args.run_fastq_gz, args.output_dir)
                )
                return demux, targets

            baseline_demux, baseline_targets = fingerprints()
            settings["summary"]["targets"] = {"COI": {}}
            pipeline_path.write_text(run_pipeline.yaml.safe_dump(settings))
            self.assertEqual(fingerprints(), (baseline_demux, baseline_targets))

            parameters = settings["ampliseq"]["targets"]["CytB"]
            original = parameters.get("sintax_cutoff")
            parameters["sintax_cutoff"] = 0.9
            pipeline_path.write_text(run_pipeline.yaml.safe_dump(settings))
            demux, targets = fingerprints()
            self.assertEqual(demux, baseline_demux)
            self.assertEqual(targets["COI"], baseline_targets["COI"])
            self.assertNotEqual(targets["CytB"], baseline_targets["CytB"])
            if original is None:
                del parameters["sintax_cutoff"]
            else:
                parameters["sintax_cutoff"] = original

            settings["threads"] += 1
            pipeline_path.write_text(run_pipeline.yaml.safe_dump(settings))
            demux, targets = fingerprints()
            self.assertNotEqual(demux, baseline_demux)
            for target in targets:
                self.assertNotEqual(targets[target], baseline_targets[target])
            settings["threads"] -= 1
            pipeline_path.write_text(run_pipeline.yaml.safe_dump(settings))

            mapping = barcodes / "Sanger_custom_barcode_F.tsv"
            mapping.write_text(mapping.read_text() + "\n")
            demux, targets = fingerprints()
            self.assertNotEqual(demux, baseline_demux)
            self.assertEqual(targets, baseline_targets)

            (directory / "CytB.fasta").write_text(">ref\nTGCA\n")
            _, targets = fingerprints()
            self.assertEqual(targets["COI"], baseline_targets["COI"])
            self.assertNotEqual(targets["CytB"], baseline_targets["CytB"])

    def test_demux_command(self) -> None:
        config, _ = load_config(FIXTURE / "config")
        command = demux_command(config, Path("run.fastq.gz"), Path("out"))
        self.assertEqual(command[:2], ["bash", str(ROOT / "scripts/demux_illumina.sh")])
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

        config["demux"]["method"] = "hybrid"
        config["demux"]["dorado_exec"] = None
        command = demux_command(config, Path("run.fastq.gz"), Path("out"))
        self.assertTrue(command[1].endswith("scripts/demux_hybrid.sh"))
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

    def test_multiple_targets_stop_on_failure_and_resume_independently(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            run_fastq = output_dir / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                profile="local",
                run_fastq_gz=run_fastq,
            )
            launches = []
            fail_target = "CytB"

            def fake_run(command, **kwargs):
                if command[0] == "nextflow":
                    target_dir = kwargs["cwd"]
                    launches.append(target_dir.name)
                    if target_dir.name == fail_target:
                        raise subprocess.CalledProcessError(1, command)
                    results = target_dir / "results"
                    results.mkdir()
                    (results / "overall_summary.tsv").write_text(
                        f"{target_dir.name}\n"
                    )
                else:
                    target = command[command.index("--target") + 1]
                    target_dir = output_dir / "ampliseq" / target
                    target_dir.mkdir(parents=True, exist_ok=True)
                    (target_dir / "samplesheet.tsv").write_text("sample\tfastq_1\n")
                    (target_dir / "excluded_samples.tsv").write_text(
                        "sample_id\treason\n"
                    )
                    (target_dir / "params.yaml").write_text(
                        f"input: {target_dir / 'samplesheet.tsv'}\n"
                        f"outdir: {target_dir / 'results'}\n"
                    )
                return subprocess.CompletedProcess(command, 0)

            with mock.patch.object(
                run_pipeline.subprocess, "run", side_effect=fake_run
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    run_ampliseq(config, args)
                self.assertEqual(launches, ["COI", "CytB"])
                self.assertTrue(
                    (output_dir / "ampliseq/COI/.completion.json").exists()
                )
                self.assertFalse(
                    (output_dir / "ampliseq/CytB/.completion.json").exists()
                )
                self.assertFalse((output_dir / "ampliseq/ITS2").exists())

                fail_target = None
                launches.clear()
                run_ampliseq(config, args)
                self.assertEqual(launches, ["CytB", "ITS2", "P3"])
                for target in config["ampliseq"]["targets"]:
                    target_dir = output_dir / "ampliseq" / target
                    self.assertTrue((target_dir / ".completion.json").exists())
                    self.assertEqual(
                        (target_dir / "results/overall_summary.tsv").read_text(),
                        f"{target}\n",
                    )

    def test_ampliseq_completion_record_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            demux_dir = output_dir / "demux"
            demux_dir.mkdir(parents=True)
            with gzip.open(
                demux_dir / "barcode02_barcode89.fastq.gz", "wb"
            ) as handle:
                handle.write(b"@r\nACGT\n+\nIIII\n")
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            config["ampliseq"]["targets"] = {"CytB": config["ampliseq"]["targets"]["CytB"]}
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
                    reads = demux_dir / "barcode02_barcode89.fastq.gz"
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
                    demux_dir / "barcode02_barcode89.fastq.gz", "wb"
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

    def test_execution_aggregates_after_ampliseq_and_after_restart(self) -> None:
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
            config, targets = load_config(FIXTURE / "config")
            config["ampliseq"]["targets"] = {"CytB": config["ampliseq"]["targets"]["CytB"]}
            real_run = subprocess.run

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
                    shutil.copytree(FIXTURE / "ampliseq/CytB", results_dir)
                    return subprocess.CompletedProcess(command, 0)
                return real_run(command, **kwargs)

            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                run_pipeline, "check_nextflow"
            ), mock.patch.object(
                run_pipeline, "load_config", return_value=(config, targets)
            ), mock.patch.object(
                run_pipeline.subprocess, "run", side_effect=fake_run
            ) as run:
                with self.assertRaises(SystemExit) as raised:
                    run_pipeline.main()
                self.assertEqual(raised.exception.code, 1)
                commands = [call.args[0] for call in run.call_args_list]
                self.assertIn("demux_illumina.sh", commands[0][1])
                self.assertIn("prepare_inputs.py", commands[1][1])
                self.assertEqual(commands[2][:3], ["nextflow", "run", "nf-core/ampliseq"])
                self.assertEqual(run.call_args_list[2].kwargs["cwd"], target_dir)
                expected_summary = [
                    sys.executable, str(ROOT / "scripts/summarise_asvs.py"),
                    "--target", "CytB", str(target_dir / "results"),
                    "--output-dir", str(Path(argv[-1]) / "summary"),
                ]
                self.assertEqual(commands[3], expected_summary)
                self.assertEqual(run.call_args_list[3].kwargs, {"check": True})
                summary = Path(argv[-1]) / "summary"
                before = {
                    path.name: path.read_bytes() for path in summary.iterdir()
                }
                self.assertEqual(
                    set(before), {"asvs.tsv", "assignments.tsv", "summarise_asvs.log"}
                )
                self.assertIn(b"0aece8f365778d45", before["asvs.tsv"])
                (summary / "assignments.tsv").unlink()
                run.reset_mock()
                with self.assertRaises(SystemExit) as raised:
                    run_pipeline.main()
                self.assertEqual(raised.exception.code, 1)
                self.assertEqual(
                    [call.args[0] for call in run.call_args_list][1:],
                    [expected_summary],
                )
                self.assertEqual(
                    {path.name: path.read_bytes() for path in summary.iterdir()},
                    before,
                )

    def test_asv_summary_failure_propagates(self) -> None:
        argv = [
            "run_pipeline.py", "run.fastq.gz", str(FIXTURE / "manifest.tsv"),
            str(FIXTURE / "config"), "out",
        ]
        failure = subprocess.CalledProcessError(2, ["summarise_asvs.py"])
        with mock.patch.object(sys, "argv", argv), mock.patch.object(
            run_pipeline, "check_nextflow"
        ), mock.patch.object(run_pipeline, "run_demux"), mock.patch.object(
            run_pipeline, "run_ampliseq"
        ), mock.patch.object(
            run_pipeline.subprocess, "run", side_effect=failure
        ) as run:
            with self.assertRaises(subprocess.CalledProcessError) as raised:
                run_pipeline.main()
        self.assertIs(raised.exception, failure)
        run.assert_called_once()
        self.assertIn("summarise_asvs.py", run.call_args.args[0][1])
        self.assertEqual(run.call_args.kwargs, {"check": True})

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
            self.assertIn("demux_illumina.sh", commands)
            self.assertIn("Ampliseq command (CytB): nextflow run", commands)
            summary_command = [
                sys.executable, str(ROOT / "scripts/summarise_asvs.py"),
            ]
            for target in ("COI", "CytB", "ITS2", "P3"):
                summary_command += [
                    "--target", target,
                    str(output_dir / "ampliseq" / target / "results"),
                ]
            summary_command += ["--output-dir", str(output_dir / "summary")]
            self.assertTrue(commands.endswith(
                f"ASV summary command: {shlex.join(summary_command)}\n"
            ))
            self.assertEqual(
                plan,
                "Method: illumina\n"
                "Targets: COI, CytB, ITS2, P3\n"
                "Samples:\n"
                "  ANLCx_002_A12: front_barcode=02, rear_barcode=89\n"
                "  BFHK_023_H8: front_barcode=08, rear_barcode=64\n"
                "  BFHK_023_A2: front_barcode=08, rear_barcode=09\n"
                "  ANLCx_001_A1: front_barcode=01, rear_barcode=01\n"
                "  ANLCx_001_B1: front_barcode=01, rear_barcode=02\n",
            )
            self.assertEqual(result.stderr, "")
            self.assertEqual(list(directory.iterdir()), [])

    def test_stage_placeholders_fail_explicitly(self) -> None:
        commands = [
            [sys.executable, str(ROOT / "scripts" / name)]
            for name in (
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
