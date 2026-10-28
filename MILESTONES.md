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

[x]

**Goal:** one small, realistic mixed-amplicon dataset supports integration
testing from demultiplexing through plots.

**Tasks**
- Developer prepares raw reads containing multiple amplicons and a manifest
  suitable for both hybrid and Illumina demultiplexing.
- Record provenance, redistribution approval, checksums and known expectations;
  distinguish these from results awaiting validation at later milestones.
- Establish the fixture as the primary integration dataset. Keep small
  synthetic tests for errors and edge cases; scaling is covered in M10.

**Done when:** the dataset and manifest load, integrity is verified, and
expected coverage of both methods and multiple targets is documented.
Demultiplexing configuration follows in M4; additional target configuration
follows in M6.

**Review gate:** developer approves the fixture and expectations. Retiring the
existing COI fixture is not required.

### M3a - Prepare mixed-target read selection [x]
- Agree which samples and targets the fixture should cover and how reads
  will be selected.
- Adapt the read-ID extraction helper for the agreed mixed-target selection;
  retain raw-read recovery with the retrieval helper.
- Update helper tests and preparation instructions.
- **Review gate:** developer approves the selection behavior before preparing
  the new dataset.

### M3b - Establish the small fixture [x]
- Developer supplies recovered raw reads and an updated manifest with real
  barcodes; the extraction ID list need not be retained.
- Rename the primary fixture directory to `tests/fixtures/small/`, retaining
  existing configuration until M4 and M6.
- Record provenance, redistribution approval, checksums and developer-approved
  expectations, including expected samples with no reads. Mark outcomes not
  yet verified by analysis.
- **Review gate:** developer approves the fixture contents and expectations.

### M3c - Align tests and documentation [x]
- Update fixture references and assertions to the approved samples, barcode
  indices and configured targets.
- Verify FASTQ integrity, unique read IDs, the approved record count,
  checksums, and loading with the existing configuration and new manifest.
- Verify dry-run output and no output creation; retain synthetic helper and
  failure tests.
- Update fixture documentation and repository references affected by the rename.
- **Review gate:** tests pass and developer approves the fixture as the
  integration baseline. Actual demultiplexing and ampliseq outcomes are
  validated in M4-M6.

## M4 - Demultiplexing

[x]

**Goal:** both demultiplexing scripts are robust, standalone, config-driven,
and ready to be called from the pipeline on the small mixed-target fixture.

**Tool provisioning**
- Conda environment is the single source of all tools except Dorado,
  including `seqtk`. Scripts call `seqtk` from `PATH`; there is no `seqtk`
  config key or CLI option.
- Dorado is not in conda, so it is the only tool located by path: a
  `--dorado` script argument, set from a `dorado` key in the config. If
  omitted, scripts use `dorado` from `PATH`.
- Scripts check at start that `dorado` and `seqtk` are executable, with a
  clear message naming the missing tool. No `module load`.
- The conda environment file lists `seqtk` and is updated where it is
  currently defined.

**Tasks**
- Remove site-specific assumptions from both scripts and make the standalone
  CLI predictable.
- Standalone correctness first, then integration with the config, driver and
  tests.
- Validate both methods on the small fixture before M5/M6.
- No target-routing or ampliseq work here.

### M4a - Refactor hybrid demultiplexing script [x]
- Convert `scripts/demux_hybrid.sh` into a standalone script with explicit
  arguments: input FASTQ, output directory, barcode directory, `--dorado`,
  `--threads`.
- Remove hard-coded paths, `module load` and implicit state; `seqtk` from
  `PATH` (see Tool provisioning).
- Keep stages explicit: ONT demux, Illumina R trimming, flatten/merge,
  counts, gzip, cleanup.
- Clear failures for missing inputs or tools, existing output directory, bad
  barcode paths and failed Dorado runs. `set -euo pipefail`, quoted variables,
  must pass `bash -n`.
- **Review gate:** developer approves the hybrid script interface and output
  layout.

### M4b - Refactor Illumina demultiplexing script [x]
- Apply the same conventions to `scripts/demux_illumina.sh`.
- Keep the logic explicit: F pass, reverse-complement unclassified reads,
  second F pass, per-F R demux, flatten, counts.
- Replace fragile filename parsing in the flatten step; validate barcode
  files and empty results.
- **Review gate:** developer approves the Illumina interface and its
  output naming and flattening rules.

### M4c - Shared demux contract [x]
- Define one output contract for both scripts: directory layout, per-sample
  FASTQ naming, counts file format, unclassified handling, cleanup and exit
  codes.
- Scripts stay independently runnable; shared conventions only, no shared
  library.
- Check both scripts against the small fixture where it applies.
- Contract is in `scripts/README.md`. Hybrid files carry an `ont_` prefix;
  stats are sorted by name.
- Manual fixture runs (real Dorado): hybrid recovers all 931 reads in the
  three expected pairs; Illumina recovers 718 of 931 in the three expected
  pairs (`barcode02_barcode89`, `barcode08_barcode09`, `barcode08_barcode64`),
  mostly lost at F matching.
- **Review gate:** developer approves the common contract and what stays
  method-specific.

### M4d - Pipeline integration and tests [x]
- Driver calls the script chosen by `demux.method`, passing the config's
  `dorado` path, barcode paths and threads; it does not look for `seqtk`.
- Tests (no Dorado needed in CI) cover command construction and the
  validation failures. Real Dorado runs on the small fixture are manual and
  recorded in the milestone notes.
- Update the README and AGENTS.md for the conda and Dorado provisioning rule.
- Config: `demux.barcode_dir` replaces per-file barcode keys; `seqtk_exec`
  and `ont_kit` are removed; `dorado_exec` is optional. The driver runs only
  the demux stage into `<output_dir>/demux`, then stops.
- **Review gate:** developer approves the integration contract before M5.

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

### M5a - Prepare ampliseq inputs [x]
- Implement `scripts/prepare_inputs.py`: per-target samplesheet from the
  manifest and the demux output, using the barcode mapping from M2.
- Samples with no demux file or an empty FASTQ are left out of the samplesheet
  and written to a separate list, so they are not lost (M7 reports them).
- Write the resolved parameter file with absolute `input` and `outdir`.
- **Review gate:** developer approves the samplesheet, the empty-sample list
  and the parameter file layout.

### M5b - Launch ampliseq per target [x]
- Driver builds the `nextflow run nf-core/ampliseq -r 3.0.0` command with
  `-profile`, `-c ampliseq.config`, a generated resource config and
  `-params-file`, run from `<output_dir>/ampliseq/<target>/`.
- Top-level `threads` shared by demux and ampliseq; the generated resource
  config sets `executor.cpus` and a CPU limit; memory limits stay in the
  user's `ampliseq.config`. `--profile local|lsf` selects the profile; there
  is no profile config key.
- Fixture `local` profile: 1 CPU and 1 GB per task as a starting point.
- Targets run one after another; a failure stops the loop.
- Before demux, the driver checks that `nextflow` is on `PATH` and is
  version 26 or higher (required by ampliseq 3.0.0).
- Tests check command construction and generated files; Nextflow is run
  manually.
- **Review gate:** developer approves the command, directory layout and
  config keys.

### M5c - Fixture run and outputs [x]
- Run one target on the small fixture with the local profile. Record the
  actual Savont output filenames and the output layout in the README.
- **Review gate:** developer approves the recorded outputs.
- Manual local/Singularity run (2026-10-28), CytB: 26 successful tasks,
  no failed or ignored tasks. Output layout and filenames recorded in README.

  | Sample | Chopper input | Chopper output | Primer-trimmed | Savont counts |
  |---|---:|---:|---:|---:|
  | ANLCx_002_A12 | 202 | 184 | 38 | 29 |
  | BFHK_023_A2 | 398 | 356 | 187 | 177 |
  | BFHK_023_H8 | 309 | 269 | 137 | 131 |

- Two ASVs: `0aece8f365778d45` (177 + 131 reads; SINTAX Homo sapiens,
  confidence 0.90) and `9e07fe1354e353c7` (29 reads; Bos indicus, 0.97).
  Both have `ampliseq_accept=TRUE`; assignments are not independently verified.
  IDs agree across counts, sequences and taxonomy. Two absent samples remain
  in `excluded_samples.tsv`.
- Future work: Savont automatically enabled low-polymorphism mode after an
  88.51% no-SNPmer warning. Earlier COI run completed but produced zero ASVs,
  an empty FASTA, header-only counts and no published taxonomy; its export
  log contained an empty-FASTA error despite workflow success. Summaries must
  distinguish these outcomes rather than equating execution success with ASVs.
- The old-pipeline counts remain sanity references, not exact acceptance
  targets. No thresholds or analysis parameters were changed during inspection.

### M5d - Restart and completion records [ ]
- Completion record per target that stores the input and config checksums.
- Rerun after an interrupt and after a config change; nothing is deleted
  automatically.
- **Review gate:** developer approves rerun behaviour.

### M5e - Review demux resume logic [ ]
- Review whether and how the demux step is skipped or validated on rerun
  (currently the demux script refuses an existing output directory).
- Decide the completion record for demux together with M5d.
- **Review gate:** developer approves demux rerun behaviour.

## M6 - Multi-target ampliseq execution

[ ]

**Goal:** every target listed in configuration runs against the demultiplexed
samples, with independent outputs and execution state.

**Tasks**
- Extend the single-target workflow to all configured targets, without
  preassigning samples or reads to targets in the driver.
- Add developer-approved fixture configurations for ITS2, CytB and P3
  alongside COI.
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
