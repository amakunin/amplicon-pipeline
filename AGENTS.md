# Development conventions

For humans and coding agents. User documentation is in [README.md](README.md);
the development plan is in [MILESTONES.md](MILESTONES.md).

## Working rules

- Work on one milestone at a time, in order. Stop at its review gate and
  summarise what changed, what was verified and any open questions.
- Do not start the next milestone, or add features beyond the current one,
  without developer approval.
- Ask rather than guess when a scientific definition, threshold, schema or
  site-specific detail is undecided. Do not invent reporting thresholds.
- Keep changes small and reviewable. Do not refactor unrelated code.
- Update the README, this file or the milestone checklist when a change makes
  them wrong.


## Repository structure

| Path | Purpose |
|---|---|
| `run_pipeline.py` | Fixed-order orchestration entrypoint |
| `scripts/demux_hybrid.sh` | Dorado hybrid demultiplexing |
| `scripts/prepare_inputs.py` | Barcode/sample mapping and per-target ampliseq samplesheets |
| `scripts/summarise_asvs.py` | Savont counts, sequences and taxonomy into stable tables |
| `scripts/summarise_samples.py` | Manifest-driven sample annotations and QC outcomes |
| `scripts/plot_summary.py` | Plots from summary tables |
| `config/example/` | Template configuration for users to copy |
| `assets/barcodes/` | Bundled barcode arrangements, sequences and mappings |
| `tests/` | Unit tests, helpers and fixtures |
| `tests/fixtures/coi_small/` | Small COI fixture; see its README |
| `tests/helpers/` | Fixture-preparation helpers, not pipeline stages |

Stage scripts other than `run_pipeline.py` are currently placeholders.

## Design decisions

- **Fixed pipeline.** Stages run in a fixed order:
  `validate -> demultiplex -> prepare inputs -> ampliseq per target -> summarise ASVs -> summarise samples -> plot`.
  No module registry, configurable DAG or plugin framework.
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

## Coding conventions

- Python 3.12, standard library only unless a milestone approves a dependency.
  Use type hints and `pathlib`.
- Scripts take explicit arguments, fail with a non-zero exit code and a clear
  message, and do not read hidden state or environment.
- Shell scripts start with `set -euo pipefail`, quote variables and must pass
  `bash -n`.
- No hard-coded site paths, usernames or unconditional `module load`. Anything
  site-specific goes in configuration.
- Prefer clear, boring code over abstractions. Comment only what needs it.
- Keep large or external data out of the repository. The SINTAX reference
  (about 1.6 GB) and real run data are referenced by path, not committed.

## Testing

Use the development environment defined in [env.yaml](env.yaml). It provides
Python 3.12 and seqtk; Nextflow/Java, the container runtime and Dorado remain
external. Do not assume terminal activation persists between agent shell calls.
When the developer supplies the environment prefix and Conda executable path,
use that executable with `conda run --no-capture-output -p <prefix>` for each
environment-dependent command. Never commit machine-specific environment paths.

```bash
python -m unittest discover -s tests -v
bash -n scripts/demux_hybrid.sh
```

- Add or update tests with each change. Use `unittest`.
- Tests must not require Dorado, Nextflow, containers or network access.
  Integration runs against the fixture are manual and recorded in the
  milestone notes.
- Do not regenerate or edit `tests/fixtures/coi_small/` without approval. Verify
  it with `sha256sum -c checksums.sha256` from that directory.
- Check actual output against expected behaviour. Do not weaken a test to make
  it pass.

## Git

- Small commits with imperative messages; one milestone per branch or PR.
- Do not commit generated outputs, `work/`, `.nextflow*`, logs or credentials.
