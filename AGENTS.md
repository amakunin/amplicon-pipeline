# Development conventions

For humans and coding agents. User documentation is in [README.md](README.md);
the development plan is in [MILESTONES.md](MILESTONES.md).

## Working rules

- Work toward one top-level milestone at a time. Use lettered submilestones
  (for example, M3a and M3b) for distinct, reviewable implementation steps.
- **Propose before coding.** Before writing code, post a short proposal in
  chat: the files to touch, function names and signatures, the data flow, what
  is deliberately not handled, and how it will be verified. Wait for explicit
  approval. If the step does not fit the step budget in MILESTONES.md, split
  it and say so.
- **No unapproved structure.** Do not add files, modules, classes, CLI options,
  config keys, dependencies or output formats that the approved proposal did
  not name. Code goes where the milestone says; if unclear, ask.
- **Do not invent rules.** Naming patterns, thresholds, defaults, schemas and
  scientific definitions come from the developer. Ask instead of choosing.
- **Prefer readable code over exhaustive code.** No custom exception classes,
  defensive checks or error aggregation unless requested. Let standard Python
  errors surface (`KeyError`, `FileNotFoundError`, ...). Add a clear message
  only where the standard error would be confusing, and say where.
- **Pass external tool errors through.** Run Dorado, seqtk, Nextflow and other
  tools without capturing or rewrapping their output, so the user sees the
  tool's own error message. Use `subprocess.run(check=True)` and let a
  non-zero exit stop the run. Do not swallow or summarise tool errors.
- **Report each step briefly:** files changed with approximate line counts,
  what was verified and how to rerun it, and open questions. Then stop at the
  review gate: do not start the next submilestone or top-level milestone, chain
  steps or run long unattended without developer approval.
- Keep changes small and reviewable. Do not refactor unrelated code.
- Update the README, this file or the milestone checklist only where a change
  makes them wrong. Do not expand documentation speculatively.
- **Be brief.** Keep chat replies, proposals and any text written into
  repository files as short as the content allows: a few sentences or bullets,
  no restating of context, and no lists of every possibility. Offer detail
  only when asked.

## Repository structure

| Path | Purpose |
|---|---|
| `run_pipeline.py` | Fixed-order orchestration entrypoint, config loading, manifest handling and the dry-run |
| `scripts/demux_illumina.sh` | Dorado Illumina F & R demultiplexing |
| `scripts/demux_hybrid.sh` | Dorado hybrid ONT and R Illumina demultiplexing |
| `scripts/prepare_inputs.py` | Barcode/sample mapping and per-target ampliseq samplesheets |
| `scripts/summarise_asvs.py` | Savont counts, sequences and taxonomy into stable tables |
| `scripts/summarise_samples.py` | Manifest-driven sample annotations and QC outcomes |
| `scripts/plot_summary.py` | Plots from summary tables |
| `config/example/` | Template configuration for users to copy |
| `assets/barcodes/` | Bundled barcode arrangements, sequences and mappings |
| `tests/` | Unit tests, helpers and fixtures |
| `tests/fixtures/small/` | Small mixed-target fixture; see its README |
| `tests/helpers/` | Fixture-preparation helpers, not pipeline stages |

The demux scripts are implemented (contract in `scripts/README.md`); the other
stage scripts are currently placeholders.

## Design decisions

- **Fixed pipeline.** Stages run in a fixed order:
  `validate -> demultiplex -> prepare inputs -> ampliseq per target -> summarise ASVs -> summarise samples -> plot`.
  No module registry, configurable DAG or plugin framework.
- **Tool provisioning.** The conda environment provides every tool except
  Dorado (including `seqtk`, taken from `PATH`; no config key). Dorado is
  set by `demux.dorado_exec`, else taken from `PATH`.
- **Ampliseq** is called directly by the driver, once per target, pinned to
  version 3.0.0, from a stable per-target execution directory. Both the
  Nextflow cache and `work/` are kept for `-resume`.
- **Configuration.** Relative resource paths resolve against the user's config
  directory. `pipeline.yaml` lists targets explicitly. Target YAML files use
  native ampliseq parameters; the driver owns `input` and `outdir` and writes
  resolved parameter files. `ampliseq.config` is passed with `-c` and the
  `local` or `lsf` profile is selected by the driver.
- **Savont, not DADA2.** Ampliseq 3.0.0's long-format counts summary is
  DADA2-only, so ASV summaries must read Savont outputs. Check the real
  published filenames on a small run before fixing adapter names.
- **Sample summaries** start from every expected manifest sample/target pair.
  Distinguish no reads, preprocessing loss, no ASVs, unresolved taxonomy, QC
  rejection and processing failure. Keep raw counts and classifications next
  to acceptance flags and reasons.
- **Restart handling** validates completion records and relevant inputs and
  config, not just whether output files exist.
- **Never delete previous results** automatically, including on failed or
  changed runs.
- **Missing samples are reported**, never silently treated as success.
- Target routing (one or several targets per demultiplexed FASTQ) and the
  manifest schema are undecided. See the milestones.
- The raw run FASTQ is never passed directly to ampliseq.
- **Demultiplexing method** is `demux.method` (hybrid or illumina); the
  manifest requirements follow it (see MILESTONES.md M2).

## Coding conventions

- Python 3.12, standard library plus PyYAML. Other dependencies need approval.
  Use type hints and `pathlib`.
- Scripts take explicit arguments, fail with a non-zero exit code and a clear
  message, and do not read hidden state or environment (tools on `PATH` excepted).
- Shell scripts start with `set -euo pipefail`, quote variables and must pass
  `bash -n`.
- No hard-coded site paths, usernames or unconditional `module load`. Anything
  site-specific goes in configuration.
- Prefer clear, boring code over abstractions. Comment only what needs it.
- Keep large or external data out of the repository. The SINTAX reference
  (about 1.6 GB) and real run data are referenced by path, not committed.

## Testing

```bash
python -m unittest discover -s tests -v
bash -n scripts/demux_hybrid.sh scripts/demux_illumina.sh
```

- Keep tests few and readable: one test file per feature, the happy path on the
  fixture, and a failure case only where a mistake would be silent or costly.
  No test-helper modules or exhaustive error matrices unless requested.
- Tests must not require Dorado, Nextflow, containers or network access.
  Integration runs against the fixture are manual and recorded in the
  milestone notes.
- Do not regenerate or edit `tests/fixtures/small/` without approval. Verify
  it with `sha256sum -c checksums.sha256` from that directory.
- Check actual output against expected behaviour. Do not weaken a test to make
  it pass.

## Git

The developer owns the repository history. Agents use git read-only.