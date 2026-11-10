"""Driver configuration, command/restart and CLI tests without external tools."""

import gzip
import hashlib
import json
import os
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
    PipelinePaths,
    TargetPaths,
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
    def test_legacy_barcode_mapping_settings_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = run_pipeline.yaml.safe_load(
                (FIXTURE / "config/pipeline.yaml").read_text(),
            )
            for key in (
                "ont_barcode_mapping", "front_barcode_mapping", "reverse_barcode_mapping",
            ):
                with self.subTest(key=key):
                    config["demux"][key] = None
                    (directory / "pipeline.yaml").write_text(
                        run_pipeline.yaml.safe_dump(config),
                    )
                    with self.assertRaisesRegex(ValueError, "place mapping TSVs in demux.barcode_dir"):
                        load_config(directory)
                    del config["demux"][key]

    def test_load_fixture_config(self) -> None:
        config, targets = load_config(FIXTURE / "config")
        self.assertEqual(config["demux"]["method"], "illumina")
        self.assertEqual(
            run_pipeline.barcode_mapping_paths(config)["rear_barcode"],
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
                "  barcode_dir: barcodes\n"
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
            self.assertEqual(config["demux"]["barcode_dir"], config_dir / "barcodes")
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
    def test_shared_paths_relocate_commands_and_fingerprint_dependencies(self) -> None:
        config, _ = load_config(FIXTURE / "config")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(Path, "open", side_effect=AssertionError("read file")):
                paths = PipelinePaths(
                    root / "raw.fastq.gz", root / "manifest.tsv",
                    root / "settings", root / "analysis", root / "driver/scripts",
                )
                target = TargetPaths(paths.ampliseq / "CytB", paths.executor_config)
                demux = demux_command(config, paths, True)
                preprocess = run_pipeline.preprocess_command(paths, True)
                prepare = run_pipeline.prepare_inputs_command(paths, "CytB")
                ampliseq = ampliseq_command(target, "lsf", True)
                summary = run_pipeline.asv_summary_command(config, paths)
                counts = run_pipeline.read_count_command(paths)
                self.assertEqual(demux[1], str(paths.demux_script("illumina")))
                self.assertEqual(demux[demux.index("--input") + 1], str(paths.run_fastq))
                self.assertEqual(preprocess[1], str(paths.preprocess_script))
                self.assertEqual(prepare[1], str(paths.prepare_inputs_script))
                self.assertEqual(summary[1], str(paths.asv_summary_script))
                self.assertEqual(counts[1], str(paths.read_count_script))
                self.assertEqual(ampliseq[ampliseq.index("-c") + 1], str(target.executor_config))
                self.assertEqual(ampliseq[-1], str(target.parameters))
                self.assertIn(str(target.results), summary)
                self.assertEqual(paths.demux_record.parent, paths.demux.parent)
                self.assertEqual(paths.preprocess_record.parent, paths.preprocessed.parent)
                self.assertEqual(target.record.parent, target.directory)
                self.assertEqual(target.work, target.directory / "work")
                self.assertEqual(target.cache, target.directory / ".nextflow/cache")

            paths.demux.mkdir(parents=True)
            target.directory.mkdir(parents=True)
            fastq = paths.preprocessed / "sample.fastq.gz"
            reference = root / "reference.fasta"
            target.samplesheet.write_text(f"sample\tfastq_1\ns1\t{fastq}\n")
            target.parameters.write_text(run_pipeline.yaml.safe_dump({
                "sintax_ref_tax_custom": str(reference),
            }))
            with (
                mock.patch.object(
                    run_pipeline, "file_fingerprints", return_value={},
                ) as files,
                mock.patch.object(
                    run_pipeline, "output_fingerprints",
                    return_value={"exists": True, "files": {}},
                ) as outputs,
            ):
                run_pipeline.demux_fingerprints(config, paths, demux)
                self.assertEqual(set(files.call_args.args[0]), {
                    paths.run_fastq, paths.manifest, paths.demux_script("illumina"),
                    *run_pipeline.barcode_mapping_paths(config).values(),
                    *(
                        config["demux"]["barcode_dir"] / f"Sanger_custom_barcode_{end}.{ext}"
                        for end in ("F", "R") for ext in ("tmol", "fasta")
                    ),
                })
                run_pipeline.preprocess_fingerprints(config, paths, preprocess)
                self.assertEqual(set(files.call_args.args[0]), {
                    paths.manifest, paths.preprocess_script, paths.driver,
                    paths.prepare_inputs_script,
                    *run_pipeline.barcode_mapping_paths(config).values(),
                })
                outputs.assert_called_once_with(paths.demux, exclude_work=True)
                hashes = run_pipeline.target_fingerprints(paths, target, ampliseq, "lsf")
                self.assertEqual(set(files.call_args.args[0]), {
                    paths.run_fastq, paths.manifest, target.executor_config,
                    target.samplesheet, target.exclusions, target.resources, fastq, reference,
                })
                self.assertEqual(
                    hashes["runtime:ampliseq-command"], run_pipeline.json_fingerprint(ampliseq),
                )
            self.assertFalse(paths.run_fastq.exists())
            self.assertFalse(paths.preprocessed.exists())

    def test_recorded_stage_records_and_cleans_only_after_success(self) -> None:
        for keep in (False, True):
            with self.subTest(keep=keep), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                output = root / "results"
                record = root / ".completion.json"
                command = ["tool", "--output", str(output)]
                fingerprints = mock.Mock(return_value={"source": "digest"})
                cleanup = mock.Mock()
                failure = subprocess.CalledProcessError(7, command)

                def execute(command, **kwargs):
                    output.mkdir(exist_ok=True)
                    (output / "published.tsv").write_text("result")
                    return subprocess.CompletedProcess(command, 0)

                with mock.patch.object(
                    run_pipeline.subprocess, "run", side_effect=failure,
                ) as run:
                    with self.assertRaises(subprocess.CalledProcessError) as raised:
                        run_pipeline.run_recorded_stage(
                            "test", command, output, record, fingerprints,
                            keep_intermediates=keep, check_integrity=False,
                            refuse_existing=False, cwd=root, cleanup=cleanup,
                        )
                    self.assertIs(raised.exception, failure)
                    cleanup.assert_not_called()
                    self.assertFalse(record.exists())
                    run.side_effect = execute
                    run_pipeline.run_recorded_stage(
                        "test", command, output, record, fingerprints,
                        keep_intermediates=keep, check_integrity=False,
                        refuse_existing=False, cwd=root, cleanup=cleanup,
                    )
                    run.assert_called_with(command, cwd=root, check=True)
                self.assertEqual(record.exists(), keep)
                self.assertEqual(cleanup.call_count, int(not keep))
                self.assertEqual(fingerprints.call_count, 2 if keep else 0)
                self.assertEqual((output / "published.tsv").read_text(), "result")

    def test_stage_completion_matches_selects_fast_or_strict_checks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output = directory / "results"
            output.mkdir()
            reads = output / "reads.fastq.gz"
            reads.write_bytes(b"reads")
            record_path = directory / ".completion.json"
            fingerprints = {"source": "digest"}
            self.assertFalse(run_pipeline.stage_completion_matches(
                record_path, output, None,
            ))
            with mock.patch.object(
                run_pipeline, "sha256_file", side_effect=AssertionError("hashed"),
            ):
                self.assertFalse(run_pipeline.stage_completion_matches(
                    record_path, output, fingerprints,
                ))
            run_pipeline.write_completion_record(
                record_path, run_pipeline.completion_record(fingerprints, output),
            )
            self.assertTrue(run_pipeline.stage_completion_matches(
                record_path, output, fingerprints,
            ))
            self.assertFalse(run_pipeline.stage_completion_matches(
                record_path, output, {"source": "changed"},
            ))
            reads.write_bytes(b"modified")
            self.assertFalse(run_pipeline.stage_completion_matches(
                record_path, output, fingerprints,
            ))
            with mock.patch.object(
                run_pipeline, "sha256_file", side_effect=AssertionError("hashed"),
            ):
                self.assertTrue(run_pipeline.stage_completion_matches(
                    record_path, output, None,
                ))
                scratch = output / "work"
                scratch.mkdir()
                (scratch / "intermediate").write_text("scratch")
                self.assertTrue(run_pipeline.stage_completion_matches(
                    record_path, output, None, exclude_work=True,
                ))
                self.assertFalse(run_pipeline.stage_completion_matches(
                    record_path, output, None,
                ))
                reads.unlink()
                self.assertFalse(run_pipeline.stage_completion_matches(
                    record_path, output, None, exclude_work=True,
                ))

    def test_runtime_fingerprints_preserve_serialization(self) -> None:
        config, _ = load_config(FIXTURE / "config")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "demux").mkdir()
            target_dir = output / "ampliseq/COI"
            target_dir.mkdir(parents=True)
            (target_dir / "samplesheet.tsv").write_text("sample\tfastq_1\n")
            parameters = {"outdir": str(target_dir / "results"), "input": "samplesheet.tsv"}
            (target_dir / "params.yaml").write_text(run_pipeline.yaml.safe_dump(parameters))
            paths = PipelinePaths(
                output / "run.fastq.gz", FIXTURE / "manifest.tsv",
                FIXTURE / "config", output,
            )
            target_paths = TargetPaths(target_dir, paths.executor_config)
            demux = demux_command(config, paths, True)
            preprocess = run_pipeline.preprocess_command(paths, True)
            ampliseq = ampliseq_command(target_paths, "local", True)
            published = {"exists": True, "files": {"reads.fastq.gz": "digest"}}
            with (
                mock.patch.object(
                    run_pipeline, "file_fingerprints",
                    side_effect=lambda paths: {"source": "digest"},
                ),
                mock.patch.object(run_pipeline, "output_fingerprints", return_value=published),
            ):
                for actual, expected in (
                    (
                        run_pipeline.demux_fingerprints(config, paths, demux),
                        {
                            "runtime:demux-config": {
                                "method": config["demux"]["method"],
                                "barcode_dir": str(config["demux"]["barcode_dir"]),
                            },
                            "runtime:demux-command": demux,
                        },
                    ),
                    (
                        run_pipeline.preprocess_fingerprints(config, paths, preprocess),
                        {
                            "runtime:preprocess-config": {
                                "method": config["demux"]["method"],
                                "threads": config["threads"],
                                "porechop_abi_options": [], "chopper_options": ["-q", "16"],
                            },
                            "runtime:preprocess-command": preprocess,
                            "runtime:preprocess-demux": published,
                        },
                    ),
                    (
                        run_pipeline.target_fingerprints(
                            paths, target_paths, ampliseq, "local",
                        ),
                        {
                            "runtime:ampliseq-params": parameters,
                            "runtime:ampliseq-command": ampliseq,
                        },
                    ),
                ):
                    legacy = {
                        key: hashlib.sha256(
                            json.dumps(value, sort_keys=isinstance(value, dict)).encode(),
                        ).hexdigest()
                        for key, value in expected.items()
                    }
                    if "runtime:ampliseq-params" in expected:
                        legacy["runtime:profile"] = hashlib.sha256(b"local").hexdigest()
                    self.assertEqual(actual, {"source": "digest", **legacy})

    def test_completion_outputs_match_checks_inventory_without_hashing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "reads.fastq.gz").write_bytes(b"reads")
            record = run_pipeline.completion_record({"source": "digest"}, output)
            (output / "reads.fastq.gz").write_bytes(b"modified")

            with mock.patch.object(
                run_pipeline, "sha256_file", side_effect=AssertionError("hashed"),
            ):
                self.assertTrue(run_pipeline.completion_outputs_match(record, output))
                (output / "unexpected.txt").write_text("extra")
                self.assertFalse(run_pipeline.completion_outputs_match(record, output))
                (output / "unexpected.txt").unlink()
                (output / "reads.fastq.gz").unlink()
                self.assertFalse(run_pipeline.completion_outputs_match(record, output))

    def test_output_fingerprints_exclude_only_top_level_work(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            for relative in ("reads.fastq.gz", "work/scratch", "results/work/published"):
                path = output / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"data")
            filtered = run_pipeline.output_fingerprints(output, exclude_work=True)
            self.assertEqual(
                set(filtered["files"]), {"reads.fastq.gz", "results/work/published"},
            )
            self.assertIn("work/scratch", run_pipeline.output_fingerprints(output)["files"])
            (output / "results/work/published").write_bytes(b"changed")
            self.assertNotEqual(
                run_pipeline.output_fingerprints(output, exclude_work=True), filtered,
            )

    def test_demux_work_cleanup_retention_and_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tools = root / "tools"
            tools.mkdir()
            dorado = tools / "dorado"
            dorado.write_text(
                f"#!{sys.executable}\n"
                "from pathlib import Path\n"
                "import sys\n"
                "output = Path(sys.argv[sys.argv.index('--output-dir') + 1])\n"
                "output.mkdir(parents=True, exist_ok=True)\n"
                "(output / 'mock_barcode01.fastq').write_bytes(Path(sys.argv[-1]).read_bytes())\n"
                "(output / 'mock_unclassified.fastq').write_text('')\n"
                f"if Path({str(root / 'fail')!r}).exists(): sys.exit(7)\n"
            )
            dorado.chmod(0o755)
            seqtk = tools / "seqtk"
            seqtk.write_text(
                f"#!{sys.executable}\n"
                "from pathlib import Path\n"
                "import sys\n"
                "sys.stdout.buffer.write(Path(sys.argv[-1]).read_bytes())\n"
            )
            seqtk.chmod(0o755)
            source = root / "input.fastq"
            source.write_text("@read\nACGT\n+\nIIII\n")
            for method in ("illumina", "hybrid"):
                for mode in ("cleanup", "keep", "failure"):
                    with self.subTest(method=method, mode=mode):
                        output = root / f"{method}-{mode}"
                        command = [
                            "bash", str(ROOT / f"scripts/demux_{method}.sh"),
                            "--input", str(source), "--output-dir", str(output),
                            "--barcode-dir", str(ROOT / "assets/barcodes"),
                            "--dorado", str(dorado),
                        ]
                        if mode == "keep":
                            command.append("--keep-intermediates")
                        if mode == "failure":
                            (root / "fail").touch()
                        with mock.patch.dict(
                            "os.environ",
                            {"PATH": f"{tools}:{os.environ['PATH']}"},
                        ):
                            if mode == "failure":
                                with self.assertRaises(subprocess.CalledProcessError):
                                    subprocess.run(command, check=True)
                            else:
                                subprocess.run(command, check=True)
                        self.assertEqual(
                            (output / "work").exists(), mode != "cleanup",
                        )
                        self.assertFalse((output / ".work").exists())
                        if mode == "failure":
                            self.assertTrue(list((output / "work").rglob("*.fastq")))
                            (root / "fail").unlink()
                        else:
                            fastqs = list(output.glob("*.fastq.gz"))
                            self.assertTrue(fastqs)
                            self.assertTrue((output / "demux_stats.txt").is_file())
                            for fastq in fastqs:
                                with gzip.open(fastq, "rb") as handle:
                                    self.assertIn(handle.read(), (b"", source.read_bytes()))

    def test_preprocess_restart_and_configuration_isolation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                keep_intermediates=True,
                check_integrity=True,
                config_dir=FIXTURE / "config", manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=directory, profile="local",
            )
            paths = PipelinePaths(
                directory / "run.fastq.gz", args.manifest_tsv, args.config_dir, directory,
            )
            demux_dir = directory / "demux"
            demux_dir.mkdir()
            reads = demux_dir / "reads.fastq.gz"
            reads.write_bytes(b"reads")
            command = run_pipeline.preprocess_command(paths, args.keep_intermediates)
            baseline = run_pipeline.preprocess_fingerprints(config, paths, command)
            scratch = demux_dir / "work"
            scratch.mkdir()
            (scratch / "intermediate.fastq").write_bytes(b"scratch")
            self.assertEqual(
                run_pipeline.preprocess_fingerprints(config, paths, command), baseline,
            )
            shutil.rmtree(scratch)
            args.keep_intermediates = False
            self.assertEqual(
                run_pipeline.preprocess_command(paths, args.keep_intermediates), command[:-1],
            )
            args.keep_intermediates = True
            config["summary"]["common"]["max_locov"] += 1
            config["ampliseq"]["targets"]["CytB"]["sintax_cutoff"] = 0.9
            args.profile = "lsf"
            self.assertEqual(
                run_pipeline.preprocess_fingerprints(config, paths, command), baseline,
            )
            config["preprocess"] = {"chopper_options": "--quality 10"}
            self.assertNotEqual(
                run_pipeline.preprocess_fingerprints(config, paths, command), baseline,
            )
            config["preprocess"] = {}
            config["threads"] += 1
            self.assertNotEqual(
                run_pipeline.preprocess_fingerprints(config, paths, command), baseline,
            )
            config["threads"] -= 1
            output = directory / "preprocessed"
            record = directory / ".preprocess-completion.json"

            def fake_run(command, **kwargs):
                output.mkdir()
                (output / "reads.fastq.gz").write_bytes(b"processed")
                (output / "excluded_samples.tsv").write_text("sample_id\treason\n")
                return subprocess.CompletedProcess(command, 0)

            with mock.patch.object(run_pipeline.subprocess, "run", side_effect=fake_run) as run:
                run_pipeline.run_preprocess(config, args, paths)
                scratch = output / "work"
                scratch.mkdir()
                (scratch / "intermediate.fastq").write_bytes(b"scratch")
                run_pipeline.run_preprocess(config, args, paths)
                shutil.rmtree(scratch)
                run_pipeline.run_preprocess(config, args, paths)
                run.assert_called_once_with(command, check=True)
                args.check_integrity = False
                with mock.patch.object(
                    run_pipeline, "sha256_file",
                    side_effect=AssertionError("fast restart hashed a file"),
                ):
                    run_pipeline.run_preprocess(config, args, paths)
                args.check_integrity = True
                before = run_pipeline.output_fingerprints(output)
                record_bytes = record.read_bytes()
                reads.write_bytes(b"changed")
                with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
                    run_pipeline.run_preprocess(config, args, paths)
                self.assertEqual(run_pipeline.output_fingerprints(output), before)
                self.assertEqual(record.read_bytes(), record_bytes)
                reads.write_bytes(b"reads")
                config["preprocess"] = {"chopper_options": "--quality 10"}
                with self.assertRaises(SystemExit):
                    run_pipeline.run_preprocess(config, args, paths)
                config["preprocess"] = {}
                (output / "reads.fastq.gz").unlink()
                with self.assertRaises(SystemExit):
                    run_pipeline.run_preprocess(config, args, paths)
                (output / "reads.fastq.gz").write_bytes(b"modified")
                with self.assertRaises(SystemExit):
                    run_pipeline.run_preprocess(config, args, paths)
                record.unlink()
                with self.assertRaises(SystemExit):
                    run_pipeline.run_preprocess(config, args, paths)
                run.assert_called_once()

            args.output_dir = directory / "failed"
            paths = PipelinePaths(
                paths.run_fastq, paths.manifest, paths.config_dir, args.output_dir,
            )
            (args.output_dir / "demux").mkdir(parents=True)
            failure = subprocess.CalledProcessError(7, command)

            def fail_run(command, **kwargs):
                (args.output_dir / "preprocessed").mkdir()
                raise failure

            with mock.patch.object(run_pipeline.subprocess, "run", side_effect=fail_run):
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    run_pipeline.run_preprocess(config, args, paths)
                self.assertIs(raised.exception, failure)
                self.assertFalse((args.output_dir / ".preprocess-completion.json").exists())
                with self.assertRaises(SystemExit):
                    run_pipeline.run_preprocess(config, args, paths)

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
            for target, parameters in settings["ampliseq"]["targets"].items():
                reference = directory / f"{target}.fasta"
                reference.write_text(">ref\nACGT\n")
                parameters["sintax_ref_tax_custom"] = str(reference)
            pipeline_path.write_text(run_pipeline.yaml.safe_dump(settings))
            args = SimpleNamespace(
                keep_intermediates=False,
                config_dir=config_dir,
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=directory / "output",
                run_fastq_gz=directory / "run.fastq.gz",
                profile="local",
            )
            args.run_fastq_gz.write_bytes(b"raw run")
            paths = PipelinePaths(
                args.run_fastq_gz, args.manifest_tsv, args.config_dir, args.output_dir,
            )
            demux_dir = args.output_dir / "demux"
            demux_dir.mkdir(parents=True)
            preprocessed_dir = args.output_dir / "preprocessed"
            preprocessed_dir.mkdir()
            (preprocessed_dir / "excluded_samples.tsv").write_text(
                "sample_id\treason\n"
                + "".join(
                    f"{sample['sample_id']}\tno_demux_file\n"
                    for sample in read_manifest(args.manifest_tsv, load_config(config_dir)[0])
                )
            )

            def fingerprints():
                config, _ = load_config(config_dir)
                targets = {}
                for target in config["ampliseq"]["targets"]:
                    subprocess.run([
                        sys.executable, str(ROOT / "scripts/prepare_inputs.py"),
                        "--manifest", str(args.manifest_tsv),
                        "--config-dir", str(config_dir),
                        "--preprocessed-dir", str(preprocessed_dir),
                        "--target", target,
                        "--output-dir", str(args.output_dir / "ampliseq"),
                    ], check=True)
                    target_dir = args.output_dir / "ampliseq" / target
                    target_paths = TargetPaths(target_dir, paths.executor_config)
                    write_resource_config(target_paths, config["threads"])
                    targets[target] = run_pipeline.target_fingerprints(
                        paths, target_paths,
                        ampliseq_command(target_paths, args.profile), args.profile,
                    )
                demux = run_pipeline.demux_fingerprints(
                    config, paths, demux_command(config, paths),
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
        paths = PipelinePaths(
            Path("run.fastq.gz"), FIXTURE / "manifest.tsv", FIXTURE / "config", Path("out"),
        )
        command = demux_command(config, paths)
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
        command = demux_command(config, paths)
        self.assertTrue(command[1].endswith("scripts/demux_hybrid.sh"))
        self.assertNotIn("--dorado", command)
        self.assertEqual(
            demux_command(config, paths, True),
            command + ["--keep-intermediates"],
        )

    def test_ampliseq_command_and_resource_config(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            target_dir = Path(temporary)
            paths = TargetPaths(target_dir, FIXTURE / "config/ampliseq.config")
            write_resource_config(paths, 4)
            self.assertEqual(
                (target_dir / "resources.config").read_text(),
                "executor.cpus = 4\nprocess.resourceLimits = [cpus: 4]\n",
            )
            self.assertEqual(
                ampliseq_command(paths, "local"),
                [
                    "nextflow", "run", "nf-core/ampliseq", "-r", "3.0.0",
                    "-profile", "local",
                    "-work-dir", str(target_dir / "work"),
                    "-c", str(FIXTURE / "config/ampliseq.config"),
                    "-c", str(target_dir / "resources.config"),
                    "-params-file", str(target_dir / "params.yaml"),
                ],
            )
            self.assertIn(
                "-resume", ampliseq_command(paths, "local", True),
            )

    def test_prepare_inputs_command_preserves_arguments(self) -> None:
        paths = PipelinePaths(
            Path("run.fastq.gz"), Path("manifest.tsv"), Path("config"), Path("out"),
        )
        self.assertEqual(
            run_pipeline.prepare_inputs_command(paths, "CytB"),
            [
                sys.executable, str(ROOT / "scripts/prepare_inputs.py"),
                "--manifest", "manifest.tsv", "--config-dir", "config",
                "--preprocessed-dir", str(Path("out").resolve() / "preprocessed"),
                "--target", "CytB",
                "--output-dir", str(Path("out").resolve() / "ampliseq"),
            ],
        )

    def test_multiple_targets_stop_on_failure_and_resume_independently(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            run_fastq = output_dir / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                keep_intermediates=True,
                check_integrity=True,
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                profile="local",
                run_fastq_gz=run_fastq,
            )
            paths = PipelinePaths(
                run_fastq, args.manifest_tsv, args.config_dir, output_dir,
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
                    self.assertEqual(
                        command[command.index("--preprocessed-dir") + 1],
                        str(output_dir / "preprocessed"),
                    )
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
                    run_ampliseq(config, args, paths)
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
                run_ampliseq(config, args, paths)
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
            demux_dir = output_dir / "preprocessed"
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
                keep_intermediates=False,
                check_integrity=True,
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                profile="local",
                run_fastq_gz=run_fastq,
            )
            paths = PipelinePaths(
                run_fastq, args.manifest_tsv, args.config_dir, output_dir,
            )
            nextflow_commands = []
            fail_nextflow = [True]

            def fake_run(command, **kwargs):
                if "prepare_inputs.py" in command[1]:
                    self.assertEqual(
                        command[command.index("--preprocessed-dir") + 1],
                        str(demux_dir),
                    )
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
                (target_dir / "work").mkdir(exist_ok=True)
                (target_dir / "work/task.txt").write_text("task")
                (target_dir / ".nextflow/cache").mkdir(parents=True, exist_ok=True)
                (target_dir / ".nextflow/cache/record").write_text("cache")
                (target_dir / ".nextflow.log").write_text("log")
                (target_dir / ".nextflow/history").write_text("metadata")
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
                    run_ampliseq(config, args, paths)
                self.assertFalse(record_path.exists())
                self.assertTrue((target_dir / "work/task.txt").exists())
                self.assertTrue((target_dir / ".nextflow/cache/record").exists())

                fail_nextflow[0] = False
                args.keep_intermediates = True
                run_ampliseq(config, args, paths)
                self.assertTrue(record_path.exists())
                self.assertTrue((target_dir / "work/task.txt").exists())
                self.assertTrue((target_dir / ".nextflow/cache/record").exists())
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 2)
                args.check_integrity = False
                with mock.patch.object(
                    run_pipeline, "sha256_file",
                    side_effect=AssertionError("fast restart hashed a file"),
                ):
                    run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 2)
                args.check_integrity = True

                (target_dir / "results/overall_summary.tsv").unlink()
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 3)

                args.profile = "lsf"
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 4)

                config["threads"] += 1
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 5)

                with gzip.open(
                    demux_dir / "barcode02_barcode89.fastq.gz", "wb"
                ) as handle:
                    handle.write(b"@r2\nTGCA\n+\nIIII\n")
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 6)

                reference.write_text(">ref\nTGCA\n")
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 7)
                args.keep_intermediates = False
                run_ampliseq(config, args, paths)
                self.assertEqual(len(nextflow_commands), 8)
                self.assertFalse((target_dir / "work").exists())
                self.assertFalse((target_dir / ".nextflow/cache").exists())
                self.assertEqual((target_dir / ".nextflow.log").read_text(), "log")
                self.assertEqual((target_dir / ".nextflow/history").read_text(), "metadata")
            self.assertNotIn("-resume", nextflow_commands[0])
            self.assertTrue(all("-resume" in command for command in nextflow_commands[1:-1]))
            self.assertNotIn("-resume", nextflow_commands[-1])

    def test_ampliseq_cleanup_refuses_links_before_deleting_work(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target"
            paths = TargetPaths(target, FIXTURE / "config/ampliseq.config")
            work = target / "work"
            cache = target / ".nextflow/cache"
            work.mkdir(parents=True)
            cache.mkdir(parents=True)
            (work / "asvs.tsv").write_text("published")
            results = target / "results"
            results.mkdir()
            output = results / "asvs.tsv"
            output.symlink_to(work / "asvs.tsv")
            with self.assertRaisesRegex(ValueError, "Published output depends"):
                run_pipeline.cleanup_ampliseq_work(paths)
            self.assertEqual(output.read_text(), "published")
            self.assertTrue(cache.exists())
            output.unlink()
            output.write_text("published copy")
            shutil.rmtree(cache)
            external = root / "external"
            external.mkdir()
            (external / "keep.txt").write_text("external cache")
            cache.symlink_to(external, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "Refusing cleanup"):
                run_pipeline.cleanup_ampliseq_work(paths)
            self.assertTrue(work.exists())
            self.assertEqual((external / "keep.txt").read_text(), "external cache")
            cache.unlink()
            run_pipeline.cleanup_ampliseq_work(paths)
            self.assertFalse(work.exists())
            self.assertEqual(output.read_text(), "published copy")

    def test_demux_completion_record_and_partial_output_guard(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            output_dir.mkdir()
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                keep_intermediates=True,
                check_integrity=True,
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                run_fastq_gz=run_fastq,
            )
            paths = PipelinePaths(
                run_fastq, args.manifest_tsv, args.config_dir, output_dir,
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
                run_demux(config, args, paths)
                record_path = output_dir / ".demux-completion.json"
                self.assertTrue(record_path.exists())
                scratch = output_dir / "demux/work"
                scratch.mkdir()
                (scratch / "intermediate.fastq").write_bytes(b"scratch")
                run_demux(config, args, paths)
                shutil.rmtree(scratch)
                run_demux(config, args, paths)
                self.assertEqual(len(calls), 1)

                (output_dir / "demux/reads.fastq.gz").unlink()
                with self.assertRaisesRegex(SystemExit, "output files are missing"):
                    run_demux(config, args, paths)
                self.assertEqual(len(calls), 1)
                self.assertTrue((output_dir / "demux/demux_stats.txt").exists())

                (output_dir / "demux/reads.fastq.gz").write_bytes(b"reads")
                run_fastq.write_bytes(b"changed raw run")
                with self.assertRaisesRegex(SystemExit, "inputs or configuration"):
                    run_demux(config, args, paths)
                self.assertEqual(len(calls), 1)

                record_path.unlink()
                with self.assertRaisesRegex(
                    SystemExit, "without a completion record"
                ):
                    run_demux(config, args, paths)
                self.assertEqual(len(calls), 1)

    def test_fast_demux_restart_uses_output_presence_without_hashing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            output_dir.mkdir()
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                keep_intermediates=True,
                check_integrity=False,
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                run_fastq_gz=run_fastq,
            )
            paths = PipelinePaths(
                run_fastq, args.manifest_tsv, args.config_dir, output_dir,
            )

            def fake_run(command, **kwargs):
                demux_dir = Path(command[command.index("--output-dir") + 1])
                demux_dir.mkdir()
                (demux_dir / "reads.fastq.gz").write_bytes(b"reads")
                (demux_dir / "demux_stats.txt").write_text("reads\t1\n")
                return subprocess.CompletedProcess(command, 0)

            with mock.patch.object(
                run_pipeline.subprocess, "run", side_effect=fake_run,
            ) as run:
                run_demux(config, args, paths)
                run_fastq.write_bytes(b"changed input")
                (output_dir / "demux/reads.fastq.gz").write_bytes(b"changed output")
                with mock.patch.object(
                    run_pipeline, "sha256_file",
                    side_effect=AssertionError("fast restart hashed a file"),
                ):
                    run_demux(config, args, paths)
                    self.assertEqual(run.call_count, 1)

                    (output_dir / "demux/unexpected.txt").write_text("extra")
                    with self.assertRaisesRegex(SystemExit, "output files are missing or unexpected"):
                        run_demux(config, args, paths)

    def test_demux_failure_does_not_write_completion_record(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "output"
            output_dir.mkdir()
            run_fastq = Path(temporary) / "run.fastq.gz"
            run_fastq.write_bytes(b"raw run")
            config, _ = load_config(FIXTURE / "config")
            args = SimpleNamespace(
                keep_intermediates=True,
                check_integrity=True,
                config_dir=FIXTURE / "config",
                manifest_tsv=FIXTURE / "manifest.tsv",
                output_dir=output_dir,
                run_fastq_gz=run_fastq,
            )
            paths = PipelinePaths(
                run_fastq, args.manifest_tsv, args.config_dir, output_dir,
            )

            def fail_after_partial_output(command, **kwargs):
                demux_dir = Path(command[command.index("--output-dir") + 1])
                (demux_dir / "work").mkdir(parents=True)
                raise subprocess.CalledProcessError(1, command)

            with mock.patch.object(
                run_pipeline.subprocess,
                "run",
                side_effect=fail_after_partial_output,
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    run_demux(config, args, paths)

            self.assertFalse((output_dir / ".demux-completion.json").exists())
            with self.assertRaisesRegex(SystemExit, "without a completion record"):
                run_demux(config, args, paths)

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
    def test_fresh_mode_never_fingerprints_or_writes_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            argv = [
                "run_pipeline.py", str(FIXTURE / "run.fastq.gz"),
                str(FIXTURE / "manifest.tsv"), str(FIXTURE / "config"), str(output),
            ]
            def run(command, **kwargs):
                if "demux_" in command[1] or "preprocess_reads.py" in command[1]:
                    Path(command[command.index("--output-dir") + 1]).mkdir(parents=True)
                elif "prepare_inputs.py" in command[1]:
                    target = command[command.index("--target") + 1]
                    (output / "ampliseq" / target).mkdir(parents=True)
                elif command[0] == "nextflow":
                    (kwargs["cwd"] / "results").mkdir()
                    self.assertNotIn("-resume", command)
                return subprocess.CompletedProcess(command, 0)

            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(run_pipeline, "check_nextflow"),
                mock.patch.object(run_pipeline.subprocess, "run", side_effect=run),
                mock.patch.object(run_pipeline, "sha256_file", side_effect=AssertionError("hashing")),
                mock.patch.object(run_pipeline, "output_fingerprints", side_effect=AssertionError("hashing")),
                mock.patch.object(run_pipeline, "demux_fingerprints", side_effect=AssertionError("hashing")),
                mock.patch.object(run_pipeline, "preprocess_fingerprints", side_effect=AssertionError("hashing")),
                mock.patch.object(run_pipeline, "target_fingerprints", side_effect=AssertionError("hashing")),
                mock.patch.object(run_pipeline, "write_completion_record") as write,
                mock.patch("builtins.print") as progress,
            ):
                with self.assertRaises(SystemExit) as raised:
                    run_pipeline.main()
                self.assertEqual(raised.exception.code, 1)
                write.assert_not_called()
                progress.assert_any_call("Validating configuration and manifest", flush=True)
                progress.assert_any_call("Input validation complete: 5 samples", flush=True)
                self.assertFalse(list(output.rglob("*completion.json")))
                with mock.patch.object(run_pipeline.subprocess, "run") as commands:
                    with self.assertRaises(SystemExit):
                        run_pipeline.main()
                    commands.assert_not_called()

    def test_hashing_progress_identifies_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input"
            path.write_bytes(b"abc")
            with mock.patch("builtins.print") as progress:
                self.assertEqual(
                    run_pipeline.sha256_file(path),
                    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
                )
            progress.assert_any_call(f"Fingerprinting: {path} (3 bytes)", flush=True)
            progress.assert_any_call(f"Fingerprint complete: {path}", flush=True)

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
                "--keep-intermediates",
                str(Path(temporary) / "run.fastq.gz"),
                str(FIXTURE / "manifest.tsv"),
                str(FIXTURE / "config"),
                str(Path(temporary) / "out"),
            ]
            shutil.copyfile(FIXTURE / "run.fastq.gz", argv[2])
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
                    (demux_dir / "demux_stats.txt").write_text("barcode08_barcode09\t931\n")
                    return subprocess.CompletedProcess(command, 0)
                if "preprocess_reads.py" in command[1]:
                    preprocessed_dir = Path(command[command.index("--output-dir") + 1])
                    preprocessed_dir.mkdir()
                    (preprocessed_dir / "reads.fastq.gz").write_bytes(b"reads")
                    (preprocessed_dir / "excluded_samples.tsv").write_text(
                        "sample_id\treason\n"
                    )
                    (preprocessed_dir / "read_counts.tsv").write_text(
                        "sample_id\tinput_reads\tporechop_reads\tchopper_reads\n"
                        "BFHK_023_A2\t931\t342\t314\n"
                        "BFHK_023_H8\t\t\t\nANLCx_002_A12\t\t\t\n"
                        "ANLCx_001_A1\t\t\t\nANLCx_001_B1\t\t\t\n"
                    )
                    return subprocess.CompletedProcess(command, 0)
                if "prepare_inputs.py" in command[1]:
                    target_dir.mkdir(parents=True, exist_ok=True)
                    reads = Path(argv[-1]) / "preprocessed/reads.fastq.gz"
                    self.assertEqual(
                        command[command.index("--preprocessed-dir") + 1],
                        str(reads.parent),
                    )
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
                if command[1].endswith("count_reads.py"):
                    from scripts.count_reads import count_reads
                    with mock.patch(
                        "scripts.count_reads.load_config",
                        return_value=(config, config["ampliseq"]["targets"]),
                    ):
                        count_reads(Path(argv[2]), FIXTURE / "manifest.tsv",
                                    FIXTURE / "config", Path(argv[-1]))
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
                self.assertIn("preprocess_reads.py", commands[1][1])
                self.assertIn("prepare_inputs.py", commands[2][1])
                self.assertEqual(commands[3][:3], ["nextflow", "run", "nf-core/ampliseq"])
                self.assertEqual(run.call_args_list[3].kwargs["cwd"], target_dir)
                expected_summary = [
                    sys.executable, str(ROOT / "scripts/summarise_asvs.py"),
                    "--target", "CytB", str(target_dir / "results"),
                    "--output-dir", str(Path(argv[-1]) / "summary"),
                    "--config-dir", str(FIXTURE / "config"),
                ]
                self.assertEqual(commands[4], expected_summary)
                expected_counts = run_pipeline.read_count_command(PipelinePaths(
                    Path(argv[2]), FIXTURE / "manifest.tsv",
                    FIXTURE / "config", Path(argv[-1]),
                ))
                self.assertEqual(commands[5], expected_counts)
                self.assertEqual(len(commands), 6)
                self.assertEqual(run.call_args_list[4].kwargs, {"check": True})
                summary = Path(argv[-1]) / "summary"
                before = {
                    path.name: path.read_bytes() for path in summary.iterdir()
                }
                self.assertEqual(
                    set(before), {"asvs.tsv", "assignments.tsv", "summarise_asvs.log", "read_flow.tsv"}
                )
                self.assertIn(b"0aece8f365778d45", before["asvs.tsv"])
                self.assertIn(b"length_status", before["asvs.tsv"])
                self.assertIn(b"contam_possible\tlocov", before["assignments.tsv"])
                self.assertIn(b"150\tfalse\tfalse", before["assignments.tsv"])
                self.assertIn(b"12\tfalse\ttrue", before["assignments.tsv"])
                (summary / "assignments.tsv").unlink()
                (summary / "read_flow.tsv").unlink()
                run.reset_mock()
                with self.assertRaises(SystemExit) as raised:
                    run_pipeline.main()
                self.assertEqual(raised.exception.code, 1)
                self.assertEqual(
                    [call.args[0] for call in run.call_args_list][1:],
                    [expected_summary, expected_counts],
                )
                self.assertEqual(
                    {path.name: path.read_bytes() for path in summary.iterdir()},
                    before,
                )

    def test_asv_summary_failure_propagates(self) -> None:
        argv = [
            "run_pipeline.py", str(FIXTURE / "run.fastq.gz"), str(FIXTURE / "manifest.tsv"),
            str(FIXTURE / "config"), "out",
        ]
        failure = subprocess.CalledProcessError(2, ["summarise_asvs.py"])
        with mock.patch.object(sys, "argv", argv), mock.patch.object(
            run_pipeline, "check_nextflow"
        ), mock.patch.object(run_pipeline, "run_demux"), mock.patch.object(
            run_pipeline, "run_preprocess"
        ), mock.patch.object(
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

    def test_read_count_failure_propagates(self) -> None:
        failure = subprocess.CalledProcessError(2, ["count_reads.py"])
        argv = [
            "run_pipeline.py", str(FIXTURE / "run.fastq.gz"), str(FIXTURE / "manifest.tsv"),
            str(FIXTURE / "config"), "out",
        ]
        with (
            mock.patch.object(sys, "argv", argv),
            mock.patch.object(run_pipeline, "check_nextflow"),
            mock.patch.object(run_pipeline, "run_demux"),
            mock.patch.object(run_pipeline, "run_preprocess"),
            mock.patch.object(run_pipeline, "run_ampliseq"),
            mock.patch.object(run_pipeline.subprocess, "run", side_effect=[
                subprocess.CompletedProcess([], 0), failure,
            ]) as run,
        ):
            with self.assertRaises(subprocess.CalledProcessError) as raised:
                run_pipeline.main()
        self.assertIs(raised.exception, failure)
        self.assertEqual(run.call_count, 2)
        self.assertTrue(run.call_args.args[0][1].endswith("count_reads.py"))

    def test_preprocess_failure_stops_before_ampliseq(self) -> None:
        argv = [
            "run_pipeline.py", str(FIXTURE / "run.fastq.gz"), str(FIXTURE / "manifest.tsv"),
            str(FIXTURE / "config"), "out",
        ]
        failure = subprocess.CalledProcessError(7, ["preprocess_reads.py"])
        with (
            mock.patch.object(sys, "argv", argv),
            mock.patch.object(run_pipeline, "check_nextflow"),
            mock.patch.object(run_pipeline, "run_demux"),
            mock.patch.object(run_pipeline, "run_preprocess", side_effect=failure),
            mock.patch.object(run_pipeline, "run_ampliseq") as ampliseq,
        ):
            with self.assertRaises(subprocess.CalledProcessError) as raised:
                run_pipeline.main()
            self.assertIs(raised.exception, failure)
            ampliseq.assert_not_called()

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
            preprocess = run_pipeline.preprocess_command(PipelinePaths(
                directory / "run.fastq.gz", FIXTURE / "manifest.tsv",
                FIXTURE / "config", output_dir,
            ))
            self.assertIn(f"Preprocess command: {shlex.join(preprocess)}", commands)
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
            summary_command += ["--config-dir", str(FIXTURE / "config")]
            self.assertIn(
                f"ASV summary command: {shlex.join(summary_command)}\n", commands,
            )
            expected_counts = run_pipeline.read_count_command(PipelinePaths(
                directory / "run.fastq.gz", FIXTURE / "manifest.tsv",
                FIXTURE / "config", output_dir,
            ))
            self.assertTrue(commands.endswith(
                f"Read count command: {shlex.join(expected_counts)}\n"
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

    def test_dry_run_forwards_keep_intermediates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            result = subprocess.run(
                [
                    sys.executable, str(ROOT / "run_pipeline.py"), "run.fastq.gz",
                    str(FIXTURE / "manifest.tsv"), str(FIXTURE / "config"),
                    str(output), "--dry-run", "--keep-intermediates",
                ],
                capture_output=True, text=True, check=True,
            )
            commands = result.stdout.splitlines()
            for prefix in ("Demux command:", "Preprocess command:"):
                command = next(line for line in commands if line.startswith(prefix))
                self.assertIn("--keep-intermediates", shlex.split(command))
            ampliseq = [line for line in commands if line.startswith("Ampliseq command")]
            self.assertEqual(len(ampliseq), 4)
            self.assertTrue(all("-resume" in shlex.split(line) for line in ampliseq))
            self.assertFalse(output.exists())

    def test_dry_run_check_integrity_implies_keep_intermediates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            result = subprocess.run(
                [
                    sys.executable, str(ROOT / "run_pipeline.py"), "run.fastq.gz",
                    str(FIXTURE / "manifest.tsv"), str(FIXTURE / "config"),
                    str(output), "--dry-run", "--check-integrity",
                ],
                capture_output=True, text=True, check=True,
            )
            commands = result.stdout.splitlines()
            for prefix in ("Demux command:", "Preprocess command:"):
                command = next(line for line in commands if line.startswith(prefix))
                self.assertIn("--keep-intermediates", shlex.split(command))
            ampliseq = [line for line in commands if line.startswith("Ampliseq command")]
            self.assertTrue(all("-resume" in shlex.split(line) for line in ampliseq))
            self.assertFalse(output.exists())

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
