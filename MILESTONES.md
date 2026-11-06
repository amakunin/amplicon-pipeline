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
  synthetic tests for errors and edge cases; scaling is covered in M12.

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

[x]

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
- Historical run evidence, counts and empty-COI caveats are recorded in the
  [fixture README](tests/fixtures/small/README.md#historical-m5c-evidence).

### M5d - Restart and completion records [x]
- Per-target `.completion.json` stores SHA-256 fingerprints of the raw run,
  manifest, target/config files, configured reference resources and generated
  ampliseq inputs, plus the selected Nextflow profile and published result
  files.
- A matching record skips the target; a missing or stale record runs Nextflow
  with `-resume`. The record is replaced only after success, and no files are
  deleted automatically. Missing or changed published outputs trigger a rerun.
- `prepare_inputs.py` can rewrite its generated inputs in an existing target
  directory.
- **Review gate:** developer approves rerun behaviour.

### M5e - Review demux resume logic [x]
- Driver-owned `.demux-completion.json` stores fingerprints for the run,
  manifest, demux config, selected barcode mappings/resources and demux script,
  plus the generated output-file fingerprints.
- Skip demux only when both inputs/config and all generated output files match.
  If the output directory exists but is incomplete, stale or has no record,
  fail clearly without changing or deleting it. If the output directory is
  absent, run demux; write the record only after success.
- Refactoring complete: shared completion-record construction; driver tests
  grouped into configuration/manifest, command/restart and CLI classes in
  one file. Restart policy and test architecture documented in AGENTS.md.
  All 24 tests pass under Python 3.12; test bodies unchanged.
- Manual review passed (2026-11-02): developer removed selected demux and
  ampliseq outputs and confirmed the expected skip, safe-failure and resume
  behaviour.

## M6 - Multi-target ampliseq execution

[x]

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

### M6a - Align target configurations [x]
- Update the existing small fixture config to enable COI, CytB, ITS2 and P3
  in alphabetical order; use developer-supplied P3 parameters.
- Arrange example and fixture target YAML sections consistently: input,
  primers, ASV, SINTAX, skips, ignores. Add brief section annotations; the
  input section documents driver-owned `input`/`outdir` and sets
  `sequencing_type`.
- Preserve scientific parameter values; remove ITS2's explicit `input` and
  duplicate skip keys.
- Implementation ready for review: all five target YAMLs use the six sections;
  P3's driver-owned `input`/`outdir` removed and four fixture targets enabled.
  Parsed parameter values preserved; duplicate-key/section-order checks and
  all 24 Python 3.12 tests pass. Existing single-target mocks remain isolated.
- **Review gate:** developer approves the four target configs and annotations.

### M6b - Verify sequential multi-target execution [x]
- Keep target runs sequential for the local executor, in configured order,
  with separate execution directories and completion records.
- Extend `CommandRestartTests` to check target order, independent state and
  stopping before the next target after failure. Update fixture expectations
  for the four configured targets; no parallel execution or new options.
- Existing synchronous driver loop retained. One four-target test checks
  failure stopping, separate results/records and skipping a completed target
  while resuming the remaining targets.
- **Review gate:** developer approves orchestration and tests.

### M6c - Validate the four-target fixture [x]
- Developer runs COI, CytB, ITS2 and P3 and records published outputs,
  per-target exclusions and rerun behaviour.
- COI is expected to yield no ASVs; compare CytB with M5c. At least one ITS2
  sequence is hoped for, not an acceptance threshold. Agree P3 expectations
  with its configuration; do not tune parameters merely to force success.
- Four-target run captured after corrected P3 taxonomic levels: only overall
  summary, Savont feature table and processed/raw SINTAX tables retained.
  Observed outcomes and historical M5c evidence are in the
  [fixture README](tests/fixtures/small/README.md#captured-multi-target-outputs).
  Capture integrity and ID consistency are checked without running analysis.
- Developer approved the captured fixtures and observed outcomes (2026-11-02);
  all 26 Python 3.12 tests and fixture checksum checks pass.
- Primary capture replaced with the approved Illumina run `results/m7a_test`,
  including demux statistics, to cover partial/unexpected barcode assignments
  and low-confidence raw taxonomy. Historical hybrid evidence remains in the
  fixture README; fixture-specific Savont settings are investigated in M12.
- **Review gate:** developer approves observed outcomes and Savont output
  layouts before M7 summary schemas are defined.

## M7 - ASV summaries

[x]

**Goal:** stable ASV count, sequence and taxonomy tables from Savont outputs.

**Tasks:** implement `scripts/summarise_asvs.py` against the outputs verified
in M5-M6, with tests on small captured output.

**Review gate:** developer approves the table schemas.

### M7a - Standalone aggregation

[x]

- Approved two TSV tables: unique target-ASVs with sequences, lengths,
  unchanged reported taxonomy/confidence and all three raw SINTAX fields;
  non-zero sample-target-ASV read assignments.
- Implemented `scripts/summarise_asvs.py` with repeated
  `--target TARGET RESULTS_DIR` and `--output-dir`. Taxonomy columns are the
  union of reported ranks; absent ranks are blank. No QC or interpretation.
- Test frozen four-target counts, sequences, taxonomy/raw calls, empty targets
  and rejection of inconsistent/duplicate IDs without replacing existing output.
- All 29 Python 3.12 offline tests pass, including three aggregation tests;
  syntax and editor diagnostics report no errors in the new implementation.
- Developer approved exceeding the usual review-size budget for this step.
- **Review gate:** developer reviews the implementation and output tables.

### M7b - Driver integration

[x]

- Added `asv_summary_command` in `run_pipeline.py`, shared by execution and
  dry-run. Aggregate every configured target's published results into
  `summary/asvs.tsv` and `summary/assignments.tsv` after ampliseq succeeds.
- Always regenerate summaries, even when ampliseq skips completed targets;
  no new completion records. Pass aggregation failures through unchanged.
- Test command paths/order, actual aggregation of frozen outputs after mocked
  external tools, summary regeneration on restart, dry-run and failure propagation.
- Still exit explicitly for unimplemented stages after ASV aggregation.
- All 30 Python 3.12 offline tests pass; syntax and editor diagnostics report
  no errors in the changed Python files.
- Developer completed the end-to-end rerun in `results/m7b_test` (2026-11-02).
  Sanity checks confirmed four unique target-ASVs and seven non-zero assignments:
  COI 0/0, CytB 1/162, ITS2 1/20 and P3 2/35 (ASVs/assigned reads).
  All summary sequences, taxonomy/raw fields and counts match published inputs;
  feature counts match the frozen baseline. COI uses the local COIDB copy.
  The provenance log records all input/output paths and matching target totals.
- **Review gate:** developer approves integrated execution.

## M8 - Sample summaries

[ ]

**Goal:** a per-sample/target table that explains every outcome.

**Tasks**
- Implement `scripts/summarise_samples.py` over all expected manifest pairs,
  keeping raw counts and classifications. Classifications: no reads,
  preprocessing loss, no ASVs, unresolved taxonomy, QC rejection, processing
  failure.
- Annotate all ASVs and assignments before sample-level filtering; retain
  original rows, counts and every applicable QC reason. Length and coding
  checks apply to target-ASV sequences; contamination, minimum-read and
  read-count-ratio labels apply to sample-target-ASV assignments.
- Apply sequence QC, then contamination and target-specific minimum-read
  labels. Calculate cross-sample contamination comparisons once, not
  sequentially as samples acquire labels.
- Within each sample-target, select the highest-read-count ASV without
  disqualifying labels as primary. Keep that comparator fixed while labelling
  low-read-count secondary ASVs by their ratio to the primary. Secondary does
  not itself mean problematic; if none qualify, there is no primary.
- Count ASVs without disqualifying labels, then compare each eligible
  secondary's taxonomy with the primary. Retain each comparison for mixtures
  of several ASVs; distinguish supported disagreement from missing or
  low-confidence resolution. Eligibility is not biological certainty.
- Before implementation, agree thresholds, disqualifying labels, tied-primary
  handling, contamination comparison scope and whether sequence-QC failures
  can supply its cross-sample maximum. Define target-specific length rules
  and coding checks, including genetic code, orientation and reading frame.
- Interpret sample-target outcomes before comparing signals across targets.
- Resolve taxonomy-reporting inconsistencies across amplicons, still observed
  after the M6 fixture rerun with corrected P3 taxonomic levels. Agree the
  reporting rules with the developer before implementation.

**Review gate:** developer agrees the QC definitions and denominators.

### M8a - Stage-specific restart fingerprints

[x]

- Fingerprint effective demux settings/command and per-target resolved
  ampliseq parameters/command, alongside inputs, resources, executor and
  outputs. Preserve skip, safe-failure and resume behavior; summary-only
  changes do not invalidate upstream stages, and target changes stay isolated.
- Older completion records are incompatible and are never rewritten; validate
  migrations in a fresh output directory. All 31 offline tests and syntax/editor
  checks passed.
- **Review gate:** developer approves fingerprint isolation before M8b.

### M8b - Consolidated configuration

[x]

- Consolidated configuration sections: `demux`, `ampliseq.common`,
  `ampliseq.targets`, `summary.common` and `summary.targets`, with shared
  `threads` and separate Nextflow executor configuration. Target values replace
  common values by parameter name; driver-owned `input`/`outdir` win. Resolved
  parameters remain under `ampliseq/<target>/params.yaml`; no QC thresholds
  were added.
- Migrated consumers, examples and tests without changing effective scientific
  parameters. Exact parameter equivalence passed for the example and all four
  fixture targets; offline tests, fixture checksums, syntax, diagnostics and
  whitespace checks passed. The M8b integration rerun matched M7b summary TSVs
  byte-for-byte (four ASVs, seven assignments) and verified current
  configuration, fingerprints and output inventories.
- **Review gate:** agree exact interfaces before implementation, then verify
  parameter equivalence and stage-specific invalidation.

### M8c - Sequence-level QC annotations

[x]

- `scripts/qc_asvs.py` annotates each target-ASV to `asv_qc.tsv`; raw M7 tables
  are unchanged and driver wiring is deferred. Configured checks only: COI
  lengths 640-673 nt in 3-nt steps; CytB lengths 298-316 nt in 3-nt steps
  and both use genetic code 5; P3 fixture-only range 372-380 nt; ITS2 unchecked.
- Coding evidence checks all six frames (both strands, shifts 0-2), ignores an
  incomplete terminal codon, and reports confirmed stop-free and uncertain
  frames separately. Length statuses are descriptive (`too_short`,
  `too_long`, `unexpected_length`, `pass`, `not_checked`); zero/multiple
  stop-free frames are not acceptance decisions.
- All 35 offline tests, fixture checksums, syntax, diagnostics and whitespace
  checks passed. Boundary/codon cases are covered; the P3 fixture has one
  371-nt failure and one 373-nt pass. Its 307-nt CytB ASV has two confirmed
  stop-free frames (`+0,+1`); the fixture has no COI ASV.

### M8d - Assignment-level QC annotations

[~]

- Standalone `scripts/qc_assignments.py` joins assignments to sequence QC by
  `(target, asv_id)` and writes `assignment_qc.tsv`, retaining all rows, counts,
  sequence evidence, assignment-level `qc_labels` and sequence-level
  `sequence_qc_labels`. Driver wiring is deferred.
- Compare each target-ASV across the run, including sequence-QC failures.
  If total reads exceed `min_contam_total_reads` (default 10000), assignments
  below `max_contam_prop` (default 0.01) of that total are `contam_possible`.
  Assignments below `max_locov` (default 100) are `locov`; labels may overlap.
  Settings can be overridden per target.
- Sequence labels mark length statuses `too_short`, `too_long` and
  `unexpected_length`, plus coding statuses `no_stop_free_frame` and
  `uncertain`. `pass`, `not_checked`, `one_stop_free_frame` and
  `multiple_stop_free_frames` produce no sequence label.
- Stdout and `assignment_qc.log` list only ASVs with affected samples:
  target, ASV ID, total reads, target cutoff and affected-sample count.
- Old contamination configuration keys fail explicitly rather than being
  silently ignored. Focused QC tests and fixture checksum checks pass; the
  sequence-label addition remains to be implemented and reviewed.
- **Review gate:** developer reviews the revised rules and output schema.

## M9 - Read preprocessing

[x]

**Goal:** run Porechop ABI and Chopper once per expected sample FASTQ before
ampliseq, using one Conda environment for pipeline and command-line tools.

Developer approved closure after Illumina concordance, hybrid smoke validation
and retained/default-cleanup end-to-end runs with unchanged-command repeats.

### M9a - Conda tools and configuration

[x]

- Single environment pins preprocessing tools and Nextflow; optional
  `preprocess` option strings and shared threads configure tools.
  Dorado remains separate.
- Conda dry-run resolves Python 3.12.16, Porechop ABI 0.5.0post1, Chopper
  0.12.0b and Nextflow 26.04.7 without conflicts. The recreated environment
  reports Porechop ABI 0.5.0 and Chopper 0.12.0.

### M9b - Standalone preprocessing script

[x]

- Standalone script runs Porechop ABI then Chopper once per expected nonempty
  demux sample. Missing/empty inputs remain explicit exclusions; external tool
  failures propagate. Developer approved interface and output contract.

### M9c - Pipeline integration and restart

[x]

- Porechop writes a compressed intermediate, then
  `zcat`, Chopper and `gzip` are directly connected with `subprocess.Popen`.
  Each exit code is checked; children are terminated/reaped on startup failure.
- Example and fixture `ampliseq.common` disable Porechop ABI and Chopper.
- Input preparation now requires `--preprocessed-dir`; upstream exclusions
  are preserved and empty processed FASTQs remain in the samplesheet for
  ampliseq to handle. Missing nonexcluded files fail before replacing outputs.
- Driver execution and dry-run include preprocessing; all targets consume its
  shared directory. `.preprocess-completion.json` is written atomically after
  success, skips matching runs and refuses stale/incomplete outputs unchanged.
- All 44 offline tests and raw/captured fixture checksum checks pass.
  Real Chopper output decompresses exactly, including with spaces in paths.
  Pylance and whitespace checks pass. Developer approved integration.

### M9d - Fixture validation

[x]

- Developer ran Dorado and ampliseq; agents inspected outputs.
  Initial option differences were corrected: no forced `-abi`, Chopper `-q 16`.
- Illumina `results/m9d_illumina_concordance_retry` matches historical
  preprocessing records exactly apart from order (119/314/208 reads).
  All shared target read counts, ASV sequences, sample assignments and
  processed/raw taxonomy match: COI 0, CytB 162, ITS2 20, P3 35 assigned reads.
  FASTQ integrity, shared target inputs and completion/output checks pass;
  the exit after aggregation is expected.
- Hybrid `results/m9d_hybrid` accounts for all 931 raw read IDs exactly once.
  Demux -> processed counts: ANLCx_002_A12 212 -> 184,
  BFHK_023_A2 403 -> 356, BFHK_023_H8 316 -> 269.
  FASTQs are valid; both absent samples are explicitly excluded.
  Demux work was retained; preprocessing work was cleaned after success.
  No hybrid ampliseq run was required.
- Developer approved closure. Frozen captures remain unchanged; M9e verifies
  unchanged-command skips.

### M9e - Intermediate directories and retention

[x]

- Standalone preprocessing uses `preprocessed/work/<barcode-pair>/` for the
  Porechop FASTQ and explicitly routes `--temp_dir` to its `tmp/` child.
  Published FASTQs remain outside work. `--keep-intermediates` retains work;
  otherwise the stage removes work only after all samples succeed.
  Failures retain work.
- Both demux scripts now use `demux/work/`, retaining existing success-only
  cleanup and `--keep-intermediates` behavior. Prior `.work/` outputs are
  untouched.
- Driver `--keep-intermediates` forwards to demux/preprocessing and enables
  Nextflow `-resume` only when flagged. Demux/preprocessing completion and
  downstream read fingerprints exclude top-level `work/`.
  Nextflow uses explicit per-target `-work-dir`. Without the flag, successful
  targets remove only `work/` and `.nextflow/cache/` before completion recording.
  Failures/skips retain work; logs, results, other metadata and external caches
  remain. Cleanup refuses symlinked scratch or published links into scratch.
- Five preprocessing and 24 driver tests pass, covering routing, retention,
  failure preservation, link refusal, fingerprints and skips. Both demux scripts
  pass `bash -n`; CLI, Pylance signature and whitespace checks pass.
- Developer runs `results/m9e_keep` and `results/m9e_clean` pass in both modes:
  all input/output records match, work/cache retention or removal is correct,
  logs/results remain, and FASTQ records and ASV/assignment tables match M9d.
  Each unchanged repeat skips demux, preprocessing and all four targets;
  aggregation reruns and both tables/provenance log have newer timestamps.
  The final unimplemented-stage exit is expected.
- Developer approved closure; M9d repeat-run verification is complete.

## M10 - Read-count accounting

[ ]

**Goal:** a hierarchical read-flow report from the raw run through barcode
assignment, shared preprocessing, target inference and sequence/assignment QC.
Keep this separate from ASV and sample summaries, in a dedicated script.

### M10a - Preprocessing counting

[~]

- Preprocessing streams compressed four-line records and writes
  `preprocessed/read_counts.tsv` after all samples succeed, before cleanup.
  Columns: `sample_id`, `input_reads`, `porechop_reads`, `chopper_reads`,
  in manifest order. Missing inputs have blank counts; empty inputs have 0
  input reads and blank downstream counts. Measured empty outputs have 0.
  Counts do not depend on retention; failed runs publish no completed table.
- Six preprocessing tests pass: distinct measured stage counts, missing/empty
  inputs, measured zero output, retention-mode parity, incomplete records and
  tool failure preservation. Pylance and whitespace checks pass.
- **Review gate:** developer approves the count output and verification.

### M10b - QC labelling cleanup

[~]

- QC is consolidated into `summarise_asvs.py`; standalone QC scripts are
  retired. Required `--config-dir` supplies existing thresholds. The driver
  invokes aggregation once; all QC is computed before writing canonical tables.
- `asvs.tsv` holds sequence/taxonomy, length/coding evidence, `total_reads`
  and `contam_target_cutoff` once per `(target, asv_id)`.
  `assignments.tsv` holds its key/read count and separate lowercase
  `true`/`false` flags `contam_possible` and `locov`.
  No propagated evidence, combined labels or parallel QC TSVs.
- Existing scientific checks and strict thresholds are unchanged.
  Sequence passing uses length pass/not_checked and coding not_checked or
  one/multiple stop-free frames. Assignment passing requires both flags false;
  both-QC passing combines them by target-ASV. No stored combined pass flag.
- Contamination diagnostics join `summarise_asvs.log`; historical outputs
  remain untouched. M8 QC notes describe superseded standalone interfaces.
- All 50 offline tests pass, including migrated QC boundary/override tests,
  canonical schemas and driver regeneration. Pylance/signature and whitespace
  checks pass. Aggregation-only validation on `results/m9e_keep` preserves
  all sequence/taxonomy/count values (4 ASVs, 7 assignments) and verifies
  unique-ASV totals, using temporary outputs without changing that run.
- **Review gate:** review canonical outputs before closing M10b.

### M10c - Read counts aggregation

[ ]

Depends on M10a counts and the M10b QC contract. Implement one standalone
aggregation script and one feature test; no driver wiring or plotting.

**Tasks**
- Scope any additional read-flow fixture captures with the developer before
  replacing frozen data; M9 runs are available as validation evidence.
- Aggregate the raw-input count, `demux/demux_stats.txt`, shared preprocessing
  input/output counts, per-target `overall_summary.tsv`, Savont feature counts
  and the canonical ASV/assignment QC columns.
- Represent branches explicitly: run -> first barcode -> second barcode or
  unmatched -> manifest sample -> Porechop ABI -> Chopper -> target ->
  ampliseq preprocessing -> Savont -> QC-passing assignments. Preserve
  barcode-level counts without sample matches.
- Distinguish unmatched first barcodes, matched first barcodes with unmatched
  second barcodes, unexpected first barcodes, unexpected complete pairs and
  expected complete pairs. Evaluate second-barcode expectations within their
  first-barcode context; individually expected barcodes can form an unexpected
  pair. Use expected/unexpected, not correct/misclassified: manifest membership
  does not establish biological origin.
- Reconstruct first-barcode totals from their descendants. Compare the
  demux-reported total with an independent raw-input read count; do not assume
  conservation. Reconcile adjacent stages as well as counts across targets.
- Record Porechop ABI and Chopper counts once per sample at the shared
  preprocessing stage. Keep target-level counts, including
  `cutadapt_total_processed`, target-specific; warn on disagreements and never
  silently choose one target. Preserve missing values and `NA` separately from
  zero.
- Keep primer-selection, Savont and accepted-assignment counts target-specific.
  Targets receive the same sample reads: do not sum shared preprocessing counts
  across targets or assume target branches partition the reads.
- Include read counts before QC, after sequence QC and after both sequence
  and assignment QC, using the pass rules agreed in M10b. Sum each passing
  assignment's `read_count` once per subtotal. No per-label rejection rows
  or new acceptance thresholds.
- Write a hierarchical `read_flow.tsv` with `node_id`, `parent_id`, `stage`,
  `first_barcode`, `second_barcode`, `sample_id`, `target`, `label` and
  `read_count`; no percentages initially. Final QC stages follow M10b.
- Test expected/unexpected pairs, partial classification, missing values and
  reconciliation warnings using small fixtures approved before capture.

**Observed planning evidence:** the Illumina run in `results/m7a_test`
(2026-11-02) reports 931 demux reads: 169 first-barcode unmatched, 43
second-barcode unmatched, 1 unexpected complete pair and 718 expected-pair
reads. `barcode02_barcode64` is unexpected despite both barcodes occurring in
the manifest. `barcode08_barcode64` has 244 demux reads, while all target
summaries report 234 Chopper input reads for `BFHK_023_H8`; investigate rather
than normalize away this discrepancy. These observations are not thresholds.

**Review gate:** developer approves the hierarchy, output schemas, denominators,
missing-value handling and reconciliation/warning rules before implementation,
then reviews the fixture read-flow report.

## M11 - Summary plots

[ ]

**Goal:** `scripts/plot_summary.py` plots from the summary tables.

**Tasks:** agree the plots and the plotting dependency first, then implement.

**Review gate:** developer reviews the plots on real output.

## M12 - Full-scale validation

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

## M13 - First release

[ ]

**Tasks:** update the README from the real behaviour, document site setup,
add `CITATION.cff`, tag `v0.1.0`.

## M14 - Pipeline post-release optimisation

[ ]

**Goal:** investigate parameter sensitivity on small runs before full-scale
validation; distinguish changes in read recovery from changes in ASV inference.

### M14a - compare hybrid and illumina demux on test data

[ ]

- Compare hybrid and Illumina demultiplexing, Savont strand mode and minimum
  cluster size using controlled comparisons rather than changing them together.
  The Illumina fixture uses `--single-strand` for all targets and
  `--min-cluster-size 4` for ITS2/P3; these are not established production defaults.
- Investigate ITS2 sequence differences, per-sample supporting reads and
  taxonomy confidence. Hybrid ASV `4851853f76e4a41b` (466 nt, 59 reads from
  `ANLCx_002_A12`) and Illumina ASV `71e9be273051fbdb` (527 nt, 9 reads from
  `BFHK_023_A2` and 11 from `BFHK_023_H8`) arise from different samples; do not
  assume they represent the same template or attribute differences to parameters.
- Initial comparison (2026-11-02): unit-cost global Levenshtein distance is
  210 edits in forward orientation versus 264 against the reverse complement.
  One minimum-edit forward alignment has 115 substitutions, 17 deletions and
  78 insertions relative to the hybrid sequence. This is exploratory evidence,
  not a species-distance estimate or acceptance threshold.
- Investigate possible ITS2 reconstruction failure for `ANLCx_002_A12`:
  the old pipeline reported 76 usable ITS2 reads, the hybrid run assigned 59,
  but the Illumina run has 52 primer-passing reads and zero Savont assignments.
  These counts describe different stages; inspect supporting reads and inference
  logs before concluding that coverage or a particular parameter caused failure.
- Agree comparison metrics and production parameter choices with the developer;

### M14b - Sample-target interpretation

[ ]

- Rules for primary sequence selection based on upstream labels, 
  secondary-to-primary ratio and boundary,
  disqualifying labels and primary tie handling before implementation.
  Minimum-read/primary/secondary decisions and QC driver wiring remain pending.
- Start from every expected manifest sample-target pair. Report eligible-ASV
  counts, primary/secondary interpretation and explicit no-call outcomes.
- Compare eligible secondaries with the primary using supported taxonomy;
  distinguish disagreement from insufficient resolution and retain each
  comparison when several lineages occur.
- **Review gate:** approve outcome definitions, taxonomy rules and table schema,
  then implementation and fixture results.

### M14c - Sample level cross-target interpretation

[ ]

- Agree reference-taxonomy harmonisation and cross-target comparison rules
  after sample-target interpretation; do not infer conflicts from missing ranks.
- **Review gate:** approve rules and implementation scope before coding, then
  review integrated sample summaries before M10.

### M14d - Pipeline logging

[ ]

- Agree consistent logging across every pipeline stage: concise terminal
  progress and complete per-stage/sample tool diagnostics retained in logs.
- Preserve external-tool error messages and nonzero failures; do not hide
  failures when routing verbose output.
- **Review gate:** approve log layout and terminal behavior before coding,
  then verify successful runs and failure diagnostics.
