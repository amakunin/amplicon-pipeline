# Milestones

For overall development guidelines, see [AGENTS.md](AGENTS.md).

Top-level milestones define feature goals, overall completion criteria, and major decisions. Add lettered submilestones only for meaningful implementation chunks that can be reviewed independently. Each submilestone states its scope and review gate; steps follow in order, and work stops at each gate for developer approval.

**Step budget:** a step touches at most about 3 files and adds about 100-150
lines of code, excluding the fixture and docs. If it does not fit, split it
before starting. Each step begins with a structure proposal (see AGENTS.md).

Status: `[ ]` not started, `[~]` in progress, `[x]` done and reviewed.
Future milestones remain at feature level; agree their submilestones at the
review gate before starting them.

## M0 - Scaffold (done)

- [x] Entrypoint interface, example configuration, barcode assets, COI fixture
  and scaffold tests.
- [x] README, AGENTS and MILESTONES split; MIT licence.

## M1 - Public repository readiness

[x]

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

**Goal:** `run_pipeline.py` loads and checks its inputs and stops, without
running any analysis. All code stays in `run_pipeline.py` for now.

**Decisions made:** PyYAML is the YAML library. `demux.method` in
`pipeline.yaml` is `hybrid` or `illumina`, with no default. Manifest columns:
`sample_id` always; `ont_barcode` and `rear_barcode` for `hybrid`;
`front_barcode` and `rear_barcode` for `illumina`; other columns are carried
through unvalidated. Barcodes are matched case-insensitively against the
mapping files named in `pipeline.yaml`. Sample IDs must be unique, and so must
each forward/rear barcode pair. Targets come from configuration, not the
manifest; every configured target runs against the demultiplexed samples.

### M2a - Load configuration [x]
- Load `pipeline.yaml` and the target YAML files into plain dicts. Relative
  paths resolve against the config directory.
- Check only that `demux.method` is a known value; other problems surface as
  Python errors. Path and executable checks happen when a path is first used.
- One small test: the fixture config loads.
- **Review gate:** developer approves the loading structure.

### M2b - Read the manifest [x]
- Read the manifest for the chosen method, look up barcodes in the mappings
  and report an unknown barcode with its sample name. Also report duplicate
  sample IDs and duplicate barcode pairs.
- One small test on the fixture manifest, plus one failing case.
- **Review gate:** developer approves the manifest schema as implemented.

### M2c - Dry run [x]
- Add `--dry-run`: print what was loaded (method, targets, samples with
  barcode indices) and exit without creating anything. Without it, stop after
  loading with "not implemented".
- **Done when:** the fixture config and manifest load and the dry-run prints
  the plan.
- **Review gate:** developer checks the dry-run output.

## M3 - Mixed-target integration fixture

[ ]

**Goal:** one small, realistic mixed-amplicon dataset supports integration
testing from demultiplexing through plots.

**Tasks**
- Developer prepares raw reads containing multiple amplicons, with a manifest
  and configurations for both hybrid and Illumina demultiplexing.
- Record provenance, redistribution approval, checksums and known expectations;
  distinguish these from results awaiting validation at later milestones.
- Establish the fixture as the primary integration dataset. Keep small
  synthetic tests for errors and edge cases; scaling is covered in M10.

**Done when:** the dataset and configuration load, its integrity is verified,
and expected coverage of both methods and multiple targets is documented.

**Review gate:** developer approves the fixture and expectations. Retiring the
existing COI fixture is not required.

## M4 - Hybrid and Illumina demultiplexing

[ ]

**Goal:** both demultiplexing methods produce consistent downstream inputs
and account for all expected samples.

**Tasks**
- Adapt the developer-provided hybrid and Illumina scripts, with executables
  and resources from configuration and no unconditional site setup.
- Produce per-sample FASTQs, barcode mappings and read counts in `demux/`.
- Validate both methods against the mixed-target fixture and record commands,
  versions and counts.
- Check dependencies when used, fail clearly and preserve previous results;
  establish safe rerun handling for demultiplexing.

**Done when:** both methods meet the approved fixture expectations, including
samples with no reads, and fail clearly on bad input.

**Review gate:** developer approves counts and downstream inputs for both
methods.

## M5 - Single-target ampliseq execution

[ ]

**Goal:** one target goes from raw FASTQ to published ampliseq results.

**Tasks**
- Implement `scripts/prepare_inputs.py`: target samplesheets and resolved
  parameter files with absolute `input` and `outdir`.
- Implement the ampliseq launch in the driver (version 3.0.0, `local` and
  `lsf` profiles, stable execution directory, retained cache and `work/`).
- Handle empty samples, which cannot be submitted as empty FASTQs, without
  losing track of them.
- Run one configured target on the mixed-target fixture and record the
  outputs, including the actual Savont filenames.
- Check dependencies when used and validate completion against relevant
  inputs and configuration. Test interrupted and changed-config reruns;
  never remove previous results automatically.

**Done when:** one target has the output layout in the README, safe resume is
verified, and failures stop the run with a clear message.

**Review gate:** developer approves single-target results and rerun behavior.

## M6 - Multi-target ampliseq execution

[ ]

**Goal:** every target listed in configuration runs against the demultiplexed
samples, with independent outputs and execution state.

**Tasks**
- Extend the single-target workflow to all configured targets, without
  preassigning samples or reads to targets in the driver.
- Validate the mixed-target fixture across targets, retaining sample/target
  accounting and safe rerun behavior.

**Done when:** all configured targets produce distinct results and execution
state, with empty or failed sample/target outcomes accounted for.

**Review gate:** developer inspects multi-target results and published Savont
outputs before summary schemas are defined.

## M7 - ASV summaries

[ ]

**Goal:** stable ASV count, sequence and taxonomy tables from Savont outputs.

**Tasks:** implement `scripts/summarise_asvs.py` against the outputs verified
in M5-M6, with tests on small captured output.

**Review gate:** developer approves the table schemas.

## M8 - Sample summaries

[ ]

**Goal:** a per-sample/target table that explains every outcome.

**Tasks**
- Implement `scripts/summarise_samples.py` over all expected manifest pairs,
  keeping raw counts and classifications. Classifications: no reads,
  preprocessing loss, no ASVs, unresolved taxonomy, QC rejection, processing
  failure.
- Reporting thresholds stay unset until the developer defines them.

**Review gate:** developer agrees the QC definitions and denominators.

## M9 - Summary plots

[ ]

**Goal:** `scripts/plot_summary.py` plots from the summary tables.

**Tasks:** agree the plots and the plotting dependency first, then implement.

**Review gate:** developer reviews the plots on real output.

## M10 - Full-scale validation

[ ]

**Goal:** validate the complete pipeline on a real ONT run at production scale.

**Tasks**
- Developer runs the pipeline and reviews results, resource use and
  operational behavior.
- Feed observed failures and scaling issues back into focused code fixes.
- Verify end-to-end sample/target accounting, reruns and result preservation.

**Done when:** the full run completes, results are reviewed and blocking
issues are resolved.

**Review gate:** developer approves production-scale behavior for release.

## M11 - First release

[ ]

**Tasks:** update the README from the real behaviour, document site setup,
add `CITATION.cff`, tag `v0.1.0`.
