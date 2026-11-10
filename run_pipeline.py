#!/usr/bin/env python3
"""Entrypoint scaffold for the fixed-order amplicon pipeline."""

import argparse
from collections.abc import Callable
import csv
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml


def load_config(config_dir: Path) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    config_dir = config_dir.resolve()
    with (config_dir / "pipeline.yaml").open() as handle:
        config = yaml.safe_load(handle)

    demux = config["demux"]
    if demux["method"] not in ("hybrid", "illumina"):
        raise ValueError("demux.method must be 'hybrid' or 'illumina'")
    obsolete = {
        "ont_barcode_mapping", "front_barcode_mapping", "reverse_barcode_mapping",
    } & demux.keys()
    if obsolete:
        raise ValueError(
            f"Obsolete demux settings: {', '.join(sorted(obsolete))}; "
            "remove them and place mapping TSVs in demux.barcode_dir"
        )
    for key in ("dorado_exec", "barcode_dir"):
        if key in demux and demux[key] is not None:
            demux[key] = (config_dir / demux[key]).resolve()

    targets = {}
    ampliseq = config["ampliseq"]
    for name, overrides in ampliseq["targets"].items():
        parameters = {**ampliseq["common"], **overrides}
        for key in ("sintax_ref", "sintax_ref_tax_custom"):
            if parameters.get(key) is not None:
                parameters[key] = str((config_dir / parameters[key]).resolve())
        targets[name] = parameters
    config["ampliseq"]["targets"] = targets
    return config, targets


def barcode_mapping_paths(config: dict[str, Any]) -> dict[str, Path]:
    directory = config["demux"]["barcode_dir"]
    first = (
        ("ont_barcode", "ONT_nb24_barcode.tsv")
        if config["demux"]["method"] == "hybrid"
        else ("front_barcode", "Sanger_custom_barcode_F.tsv")
    )
    return {
        first[0]: directory / first[1],
        "rear_barcode": directory / "Sanger_custom_barcode_R.tsv",
    }


def load_barcode_mapping(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {
            row["tag_sequence"].casefold(): row["index"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_manifest(
    manifest_path: Path, config: dict[str, Any]
) -> list[dict[str, Any]]:
    barcode_mappings = {
        column: load_barcode_mapping(path)
        for column, path in barcode_mapping_paths(config).items()
    }
    required_columns = ("sample_id", *barcode_mappings)
    samples = []
    sample_ids = set()
    barcode_pairs = set()

    with manifest_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for column in required_columns:
                row[column]
            sample_id = row["sample_id"]
            if sample_id in sample_ids:
                raise ValueError(f"Duplicate sample_id: {sample_id}")
            sample_ids.add(sample_id)

            indices = {}
            for column, mapping in barcode_mappings.items():
                barcode = row[column]
                try:
                    indices[column] = mapping[barcode.casefold()]
                except KeyError:
                    raise ValueError(
                        f"Unknown {column} {barcode!r} for sample {sample_id}"
                    ) from None

            barcode_pair = tuple(indices.values())
            if barcode_pair in barcode_pairs:
                raise ValueError(
                    f"Duplicate barcode pair for sample {sample_id}: {barcode_pair}"
                )
            barcode_pairs.add(barcode_pair)
            row["barcode_indices"] = indices
            samples.append(row)

    return samples


def demux_command(
    config: dict[str, Any], run_fastq: Path, output_dir: Path,
    keep_intermediates: bool = False,
) -> list[str]:
    demux = config["demux"]
    command = [
        "bash",
        str(Path(__file__).resolve().parent / f"scripts/demux_{demux['method']}.sh"),
        "--input", str(run_fastq.resolve()),
        "--output-dir", str(output_dir.resolve() / "demux"),
        "--barcode-dir", str(demux["barcode_dir"]),
        "--threads", str(config["threads"]),
    ]
    if demux.get("dorado_exec"):
        command += ["--dorado", str(demux["dorado_exec"])]
    if keep_intermediates:
        command.append("--keep-intermediates")
    return command


def preprocess_command(args: argparse.Namespace) -> list[str]:
    return [
        sys.executable,
        str(Path(__file__).resolve().parent / "scripts/preprocess_reads.py"),
        "--manifest", str(args.manifest_tsv.resolve()),
        "--config-dir", str(args.config_dir.resolve()),
        "--demux-dir", str(args.output_dir.resolve() / "demux"),
        "--output-dir", str(args.output_dir.resolve() / "preprocessed"),
        *(["--keep-intermediates"] if args.keep_intermediates else []),
    ]


def prepare_inputs_command(args: argparse.Namespace, target: str) -> list[str]:
    return [
        sys.executable,
        str(Path(__file__).resolve().parent / "scripts/prepare_inputs.py"),
        "--manifest", str(args.manifest_tsv),
        "--config-dir", str(args.config_dir),
        "--preprocessed-dir", str(args.output_dir.resolve() / "preprocessed"),
        "--target", target,
        "--output-dir", str(args.output_dir.resolve() / "ampliseq"),
    ]


def ampliseq_command(
    config_dir: Path, target_dir: Path, profile: str,
    keep_intermediates: bool = False,
) -> list[str]:
    return [
        "nextflow", "run", "nf-core/ampliseq", "-r", "3.0.0",
        *(["-resume"] if keep_intermediates else []),
        "-profile", profile,
        "-work-dir", str(target_dir.resolve() / "work"),
        "-c", str(config_dir.resolve() / "ampliseq.config"),
        "-c", str(target_dir / "resources.config"),
        "-params-file", str(target_dir / "params.yaml"),
    ]


def asv_summary_command(
    config: dict[str, Any], output_dir: Path, config_dir: Path,
) -> list[str]:
    output_dir = output_dir.resolve()
    command = [
        sys.executable,
        str(Path(__file__).resolve().parent / "scripts/summarise_asvs.py"),
    ]
    for target in config["ampliseq"]["targets"]:
        command += [
            "--target", target, str(output_dir / "ampliseq" / target / "results"),
        ]
    return command + [
        "--output-dir", str(output_dir / "summary"),
        "--config-dir", str(config_dir.resolve()),
    ]


def read_count_command(args: argparse.Namespace) -> list[str]:
    return [
        sys.executable,
        str(Path(__file__).resolve().parent / "scripts/count_reads.py"),
        "--run-fastq", str(args.run_fastq_gz.resolve()),
        "--manifest", str(args.manifest_tsv.resolve()),
        "--config-dir", str(args.config_dir.resolve()),
        "--output-dir", str(args.output_dir.resolve()),
    ]


def json_fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def sha256_file(path: Path) -> str:
    print(f"Fingerprinting: {path} ({path.stat().st_size} bytes)", flush=True)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    print(f"Fingerprint complete: {path}", flush=True)
    return digest.hexdigest()


def file_fingerprints(paths: list[Path]) -> dict[str, str]:
    return {
        str(path.resolve()): sha256_file(path)
        for path in sorted(set(paths))
    }


def output_fingerprints(
    directory: Path, exclude_work: bool = False,
) -> dict[str, Any]:
    exists = directory.is_dir()
    files = {}
    if exists:
        files = {
            str(path.relative_to(directory)): sha256_file(path)
            for path in sorted(directory.rglob("*"))
            if not (exclude_work and path.relative_to(directory).parts[0] == "work")
            and path.is_file()
        }
    return {"exists": exists, "files": files}


def write_completion_record(path: Path, record: dict[str, Any]) -> None:
    temporary_path = path.with_name(f"{path.name}.tmp")
    temporary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    temporary_path.replace(path)


def completion_record(
    fingerprints: dict[str, str], output_dir: Path, exclude_work: bool = False,
) -> dict[str, Any]:
    return {
        "fingerprints": fingerprints,
        "outputs": output_fingerprints(output_dir, exclude_work),
        "status": "complete",
    }


def completion_outputs_match(
    record: dict[str, Any], output_dir: Path, exclude_work: bool = False,
) -> bool:
    fingerprints = record.get("fingerprints")
    outputs = record.get("outputs")
    if (
        record.get("status") != "complete"
        or not isinstance(fingerprints, dict)
        or not fingerprints
        or not isinstance(outputs, dict)
        or outputs.get("exists") is not True
        or not isinstance(outputs.get("files"), dict)
        or not outputs["files"]
        or not output_dir.is_dir()
    ):
        return False

    actual_files = {
        str(path.relative_to(output_dir))
        for path in output_dir.rglob("*")
        if not (
            exclude_work and path.relative_to(output_dir).parts[0] == "work"
        )
        and path.is_file()
    }
    return actual_files == set(outputs["files"])


def stage_completion_matches(
    record_path: Path, output_dir: Path, fingerprints: dict[str, str] | None,
    *, exclude_work: bool = False,
) -> bool:
    if not record_path.exists():
        return False
    record = json.loads(record_path.read_text())
    if fingerprints is None:
        return completion_outputs_match(record, output_dir, exclude_work)
    return record == completion_record(fingerprints, output_dir, exclude_work)


def demux_fingerprints(
    config: dict[str, Any],
    args: argparse.Namespace,
    command: list[str],
) -> dict[str, str]:
    demux = config["demux"]
    method = demux["method"]
    barcode_ends = ("R",) if method == "hybrid" else ("F", "R")
    resources = list(barcode_mapping_paths(config).values())
    resources.extend(
        demux["barcode_dir"] / f"Sanger_custom_barcode_{end}.{extension}"
        for end in barcode_ends
        for extension in ("tmol", "fasta")
    )
    paths = [
        args.run_fastq_gz,
        args.manifest_tsv,
        Path(command[1]),
        *resources,
    ]
    fingerprints = file_fingerprints(paths)
    fingerprints["runtime:demux-config"] = json_fingerprint(
        {"method": method, "barcode_dir": str(demux["barcode_dir"])},
    )
    fingerprints["runtime:demux-command"] = json_fingerprint(command)
    return fingerprints


def preprocess_fingerprints(
    config: dict[str, Any], args: argparse.Namespace, command: list[str],
) -> dict[str, str]:
    demux = config["demux"]
    fingerprints = file_fingerprints([
        args.manifest_tsv, Path(command[1]),
        Path(__file__).resolve(),
        Path(__file__).resolve().parent / "scripts/prepare_inputs.py",
        *barcode_mapping_paths(config).values(),
    ])
    demux_dir = args.output_dir.resolve() / "demux"
    if not demux_dir.is_dir():
        raise FileNotFoundError(demux_dir)
    options = config.get("preprocess", {})
    effective_config = {
        "method": demux["method"], "threads": config["threads"],
        **{
            key: shlex.split(options.get(key, ""))
            for key in ("porechop_abi_options", "chopper_options")
        },
    }
    for key, value in (
        ("config", effective_config),
        ("command", command),
        ("demux", output_fingerprints(demux_dir, exclude_work=True)),
    ):
        fingerprints[f"runtime:preprocess-{key}"] = json_fingerprint(value)
    return fingerprints


def target_fingerprints(
    args: argparse.Namespace,
    config: dict[str, Any],
    target: str,
    target_dir: Path,
) -> dict[str, str]:
    samplesheet = target_dir / "samplesheet.tsv"
    paths = [
        args.run_fastq_gz,
        args.manifest_tsv,
        args.config_dir / "ampliseq.config",
        samplesheet,
        target_dir / "excluded_samples.tsv",
        target_dir / "resources.config",
    ]
    parameters = yaml.safe_load((target_dir / "params.yaml").read_text())
    for key in ("sintax_ref", "sintax_ref_tax_custom"):
        value = parameters.get(key)
        if value and value != "coidb":
            paths.append(Path(value))
    with samplesheet.open(newline="") as handle:
        paths.extend(
            Path(row["fastq_1"])
            for row in csv.DictReader(handle, delimiter="\t")
        )
    fingerprints = file_fingerprints(paths)
    fingerprints["runtime:ampliseq-params"] = json_fingerprint(parameters)
    fingerprints["runtime:ampliseq-command"] = json_fingerprint(
        ampliseq_command(
            args.config_dir, target_dir, args.profile, args.keep_intermediates,
        ),
    )
    fingerprints["runtime:profile"] = hashlib.sha256(
        args.profile.encode()
    ).hexdigest()
    return fingerprints


def check_nextflow() -> None:
    try:
        output = subprocess.run(
            ["nextflow", "-version"], capture_output=True, text=True, check=True
        ).stdout
    except FileNotFoundError:
        raise SystemExit("nextflow not found on PATH (version 26 or higher required)")
    version = re.search(r"version (\d+)", output)
    if version is None or int(version.group(1)) < 26:
        raise SystemExit(f"nextflow 26 or higher is required by ampliseq 3.0.0:{output}")


def write_resource_config(target_dir: Path, threads: int) -> Path:
    path = target_dir / "resources.config"
    path.write_text(
        f"executor.cpus = {threads}\n"
        f"process.resourceLimits = [cpus: {threads}]\n"
    )
    return path


def run_recorded_stage(
    name: str, command: list[str], output_dir: Path, record_path: Path,
    fingerprint_inputs: Callable[[], dict[str, str]], *,
    keep_intermediates: bool, check_integrity: bool, refuse_existing: bool,
    exclude_work: bool = False, cwd: Path | None = None,
    cleanup: Callable[[], None] | None = None,
) -> None:
    print(f"Starting {name}", flush=True)
    fingerprints = None
    if keep_intermediates:
        if check_integrity:
            print(f"Checking {name} fingerprints", flush=True)
            fingerprints = fingerprint_inputs()
        if stage_completion_matches(
            record_path, output_dir, fingerprints, exclude_work=exclude_work,
        ):
            reason = (
                "integrity check passed" if check_integrity
                else "recorded outputs are present"
            )
            print(f"Skipping {name}: {reason}", flush=True)
            return

    if refuse_existing and output_dir.exists():
        if not keep_intermediates:
            raise SystemExit(f"Refusing existing {name} output: {output_dir}")
        if not record_path.exists():
            reason = "without a completion record"
        elif check_integrity:
            record = json.loads(record_path.read_text())
            if record.get("fingerprints") != fingerprints:
                reason = "inputs or configuration changed"
            else:
                reason = "output files are missing, changed, or unexpected"
        else:
            reason = "output files are missing or unexpected, or completion data is invalid"
        raise SystemExit(
            f"{name} output is incomplete or stale ({reason}); "
            f"refusing to overwrite: {output_dir}"
        )

    if keep_intermediates and fingerprints is None:
        print(f"Fingerprinting {name} inputs", flush=True)
        fingerprints = fingerprint_inputs()
    if cwd is None:
        subprocess.run(command, check=True)
    else:
        subprocess.run(command, cwd=cwd, check=True)
    if not output_dir.is_dir():
        raise FileNotFoundError(output_dir)
    if not keep_intermediates and cleanup is not None:
        cleanup()
    if fingerprints is not None:
        write_completion_record(
            record_path, completion_record(fingerprints, output_dir, exclude_work),
        )
    print(f"{name} complete", flush=True)


def run_demux(config: dict[str, Any], args: argparse.Namespace) -> None:
    command = demux_command(
        config, args.run_fastq_gz, args.output_dir, args.keep_intermediates,
    )
    run_recorded_stage(
        "demultiplexing", command, args.output_dir.resolve() / "demux",
        args.output_dir.resolve() / ".demux-completion.json",
        lambda: demux_fingerprints(config, args, command),
        keep_intermediates=args.keep_intermediates,
        check_integrity=getattr(args, "check_integrity", False),
        refuse_existing=True, exclude_work=True,
    )


def run_preprocess(config: dict[str, Any], args: argparse.Namespace) -> None:
    command = preprocess_command(args)
    run_recorded_stage(
        "preprocessing", command, args.output_dir.resolve() / "preprocessed",
        args.output_dir.resolve() / ".preprocess-completion.json",
        lambda: preprocess_fingerprints(config, args, command),
        keep_intermediates=args.keep_intermediates,
        check_integrity=getattr(args, "check_integrity", False),
        refuse_existing=True, exclude_work=True,
    )


def cleanup_ampliseq_work(target_dir: Path) -> None:
    target_dir = target_dir.resolve()
    directories = [target_dir / "work", target_dir / ".nextflow/cache"]
    for directory in directories:
        if directory.is_symlink() or directory.resolve() != directory:
            raise ValueError(f"Refusing cleanup through a symlink: {directory}")
    results = target_dir / "results"
    for path in (results, *results.rglob("*")):
        if path.is_symlink() and any(
            path.resolve().is_relative_to(directory) for directory in directories
        ):
            raise ValueError(
                f"Published output depends on intermediates: {path}; "
                "use copy publication or --keep-intermediates"
            )
    for directory in directories:
        if directory.exists():
            shutil.rmtree(directory)


def run_ampliseq(
    config: dict[str, Any],
    args: argparse.Namespace,
) -> None:
    ampliseq_dir = args.output_dir.resolve() / "ampliseq"
    for target in config["ampliseq"]["targets"]:
        print(f"Preparing ampliseq inputs: {target}", flush=True)
        subprocess.run(
            prepare_inputs_command(args, target),
            check=True,
        )
        target_dir = ampliseq_dir / target
        write_resource_config(target_dir, config["threads"])
        run_recorded_stage(
            f"ampliseq for {target}",
            ampliseq_command(
                args.config_dir, target_dir, args.profile, args.keep_intermediates,
            ),
            target_dir / "results", target_dir / ".completion.json",
            lambda: target_fingerprints(args, config, target, target_dir),
            keep_intermediates=args.keep_intermediates,
            check_integrity=getattr(args, "check_integrity", False),
            refuse_existing=False, cwd=target_dir,
            cleanup=lambda: cleanup_ampliseq_work(target_dir),
        )


def print_plan(config: dict[str, Any], samples: list[dict[str, Any]]) -> None:
    print(f"Method: {config['demux']['method']}")
    print(f"Targets: {', '.join(config['ampliseq']['targets'])}")
    print("Samples:")
    for sample in samples:
        indices = ", ".join(
            f"{name}={index}"
            for name, index in sample["barcode_indices"].items()
        )
        print(f"  {sample['sample_id']}: {indices}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "ONT amplicon pipeline scaffold. Demultiplexing, preprocessing, ampliseq and ASV "
            "aggregation and QC run; later stages are not implemented."
        )
    )
    parser.add_argument("run_fastq_gz", type=Path, help="Raw run FASTQ (.fastq.gz)")
    parser.add_argument("manifest_tsv", type=Path, help="Sample manifest (.tsv)")
    parser.add_argument("config_dir", type=Path, help="Run configuration directory")
    parser.add_argument("output_dir", type=Path, help="Analysis output directory")
    parser.add_argument(
        "--profile",
        choices=("local", "lsf"),
        default="local",
        help="Nextflow profile from ampliseq.config (default: local)",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Retain stage work/cache and enable stage skips and Nextflow resume",
    )
    parser.add_argument(
        "--check-integrity", action="store_true",
        help="Enable full input/output fingerprint checks (also enables --keep-intermediates)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the loaded plan without creating outputs or running analysis",
    )
    args = parser.parse_args()
    if args.check_integrity:
        args.keep_intermediates = True
    if not args.dry_run:
        mode = (
            "integrity-checked restart" if args.check_integrity
            else "fast retained restart" if args.keep_intermediates
            else "fresh, no fingerprinting"
        )
        print(
            f"Starting pipeline: input={args.run_fastq_gz.resolve()}; "
            f"manifest={args.manifest_tsv.resolve()}; config={args.config_dir.resolve()}; "
            f"output={args.output_dir.resolve()}; profile={args.profile}; "
            f"mode={mode}",
            flush=True,
        )
        if not args.keep_intermediates and args.output_dir.exists():
            parser.exit(status=1, message=f"Output directory exists; refusing to overwrite: {args.output_dir}\n")
        print("Validating configuration and manifest", flush=True)
    config, _ = load_config(args.config_dir)
    samples = read_manifest(args.manifest_tsv, config)
    if not args.dry_run:
        for path in (args.run_fastq_gz, args.config_dir / "ampliseq.config"):
            if not path.is_file():
                raise FileNotFoundError(path)
        print(f"Input validation complete: {len(samples)} samples", flush=True)
    command = demux_command(
        config, args.run_fastq_gz, args.output_dir, args.keep_intermediates,
    )
    summary_command = asv_summary_command(config, args.output_dir, args.config_dir)
    if args.dry_run:
        print_plan(config, samples)
        print(f"Demux command: {shlex.join(command)}")
        print(f"Preprocess command: {shlex.join(preprocess_command(args))}")
        for target in config["ampliseq"]["targets"]:
            target_dir = args.output_dir.resolve() / "ampliseq" / target
            ampliseq = ampliseq_command(
                args.config_dir, target_dir, args.profile, args.keep_intermediates,
            )
            print(f"Ampliseq command ({target}): {shlex.join(ampliseq)}")
        print(f"ASV summary command: {shlex.join(summary_command)}")
        print(f"Read count command: {shlex.join(read_count_command(args))}")
        return
    print("Checking Nextflow", flush=True)
    check_nextflow()
    print("Nextflow check complete", flush=True)
    run_demux(config, args)
    run_preprocess(config, args)
    run_ampliseq(config, args)
    print("Starting ASV aggregation and QC", flush=True)
    subprocess.run(summary_command, check=True)
    print("Starting read-flow accounting", flush=True)
    subprocess.run(read_count_command(args), check=True)
    parser.exit(
        status=1,
        message="Stages after ASV aggregation are not implemented.\n",
    )


if __name__ == "__main__":
    main()
