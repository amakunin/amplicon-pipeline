# Milestones

Each milestone is implemented by an agent as one unit of work, following
[AGENTS.md](AGENTS.md). After each one the developer reviews the code, records
decisions and adjusts the later milestones before the next starts.

Status: `[ ]` not started, `[~]` in progress, `[x]` done and reviewed.
Milestones after M3 are deliberately less detailed; refine them at the review
gate before them.

## M0 - Scaffold (done)

- [x] Entrypoint interface, example configuration, barcode assets, COI fixture
  and scaffold tests.
- [x] README, AGENTS and MILESTONES split; MIT licence.

## M1 - Public repository readiness

[~]

**Goal:** the repository can be published without leaking site details.

**Tasks**
- Replace site-specific absolute paths in the fixture configuration and barcode
  documentation with placeholders or descriptions.
- Confirm that bundled assets and fixture data may be public.
- Extend `.gitignore` (`*.log`, `.ipynb_checkpoints/`, `.vscode/`).
- Add a CI workflow running the Tests commands from AGENTS.md.

**Done when:** a repository-wide scan finds no site-specific filesystem paths,
tests pass and CI is green.

**Review gate:** developer confirms the data may be public and creates the repo.

## M2 - Configuration and manifest validation

[ ]

**Goal:** `run_pipeline.py` validates its inputs and stops, without running
any analysis.

**Tasks**
- Load `pipeline.yaml` and target YAML files. YAML needs a decision first:
  add PyYAML, or write a minimal parser.
- Propose and document the manifest schema. Start from the columns in
  `tests/fixtures/coi_small/manifest.tsv`.
- Validate paths, required fields, duplicate samples, and barcode agreement
  between the manifest and `assets/barcodes/`.
- Add `--dry-run` that prints the resolved plan.

**Done when:** the fixture config and manifest validate and bad inputs give
clear errors, covered by unit tests.

**Review gate:** developer approves the manifest schema and the YAML decision.

## M3 - Demultiplexing

[ ]

**Goal:** `scripts/demux_hybrid.sh` demultiplexes the fixture.

**Tasks**
- Adapt the prototype's hybrid Dorado logic, with executables and resources
  from configuration and no `module load`.
- Produce per-sample FASTQs, a barcode mapping and read counts in `demux/`.
- Compare counts with `tests/fixtures/coi_small/expected_demux.tsv`.
- Record the manual run, including commands and versions.

**Done when:** the fixture gives positive counts for A1, C1 and E1 and zero
for B1 and D1, and the script fails clearly on bad input.

**Review gate:** developer checks the counts and the target-routing evidence
below.

## M4 - Target routing investigation

[ ]

**Goal:** find out whether a demultiplexed FASTQ holds one target or several,
before connecting it to ampliseq.

**Tasks**
- Report only: summarise primer matches per sample on the fixture and on any
  mixed-target data the developer provides.
- Propose a routing policy with its limits. Do not implement it.

**Review gate:** developer decides the routing policy. This shapes M5 onward.

## M5 - Single-target vertical slice (COI)

[ ]

**Goal:** one target goes from raw FASTQ to published ampliseq results.

**Tasks**
- Implement `scripts/prepare_inputs.py`: target samplesheets and resolved
  parameter files with absolute `input` and `outdir`.
- Implement the ampliseq launch in the driver (version 3.0.0, `local` and
  `lsf` profiles, stable execution directory, retained cache and `work/`).
- Handle empty samples, which cannot be submitted as empty FASTQs, without
  losing track of them.
- Run the fixture end to end on LSF and record the outputs, including the
  actual Savont filenames.

**Done when:** `output_dir` has the layout in the README for COI and failures
stop the run with a clear message.

**Review gate:** developer inspects the real Savont outputs. Do not start
summaries before then.

## M6 - ASV summaries

[ ]

**Goal:** stable ASV count, sequence and taxonomy tables from Savont outputs.

**Tasks:** implement `scripts/summarise_asvs.py` against the filenames found in
M5, with tests on small captured output.

**Review gate:** developer approves the table schemas.

## M7 - Sample summaries

[ ]

**Goal:** a per-sample/target table that explains every outcome.

**Tasks**
- Implement `scripts/summarise_samples.py` over all expected manifest pairs,
  keeping raw counts and classifications. Classifications: no reads,
  preprocessing loss, no ASVs, unresolved taxonomy, QC rejection, processing
  failure.
- Reporting thresholds stay unset until the developer defines them.

**Review gate:** developer agrees the QC definitions and denominators.

## M8 - Plots

[ ]

**Goal:** `scripts/plot_summary.py` plots from the summary tables.

**Tasks:** agree the plots and the plotting dependency first, then implement.

**Review gate:** developer reviews the plots on real output.

## M9 - Robustness

[ ]

**Goal:** safe reruns and clear failures.

**Tasks**
- Completion records, with checks of inputs and config on restart. Resume
  must not rely on the existence of output files alone.
- Dependency checks (Nextflow, container runtime, Dorado, seqtk).
- Never remove previous results automatically.
- Multiple targets in one run, if the routing policy needs them.

**Review gate:** developer tests an interrupted and a changed-config rerun.

## M10 - First release

[ ]

**Tasks:** update the README from the real behaviour, document site setup,
add `CITATION.cff`, tag `v0.1.0`.
