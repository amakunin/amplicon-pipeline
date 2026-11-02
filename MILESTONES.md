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
  fixture README; fixture-specific Savont settings are investigated in M11.
- **Review gate:** developer approves observed outcomes and Savont output
  layouts before M7 summary schemas are defined.

## M7 - ASV summaries

[~]

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

[~]

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

## M9 - Read-count accounting

[ ]

**Goal:** a hierarchical read-flow report from the raw run through barcode
assignment, sample preprocessing, target inference and M8 acceptance decisions.
Keep this separate from ASV and sample summaries, in a dedicated script.

**Tasks**
- Aggregate `demux/demux_stats.txt`, per-target `overall_summary.tsv`, Savont
  feature counts and downstream accepted assignments. Agree the script
  interface and output schemas before implementation.
- Represent branches explicitly: run -> first barcode -> second barcode or
  unmatched -> manifest sample -> preprocessing -> target -> Savont ->
  accepted assignments. Preserve barcode-level counts without sample matches.
- Distinguish unmatched first barcodes, matched first barcodes with unmatched
  second barcodes, unexpected first barcodes, unexpected complete pairs and
  expected complete pairs. Evaluate second-barcode expectations within their
  first-barcode context; individually expected barcodes can form an unexpected
  pair. Use expected/unexpected, not correct/misclassified: manifest membership
  does not establish biological origin.
- Reconstruct first-barcode totals from their descendants. Compare the
  demux-reported total with an independent raw-input read count; do not assume
  conservation. Reconcile adjacent stages as well as counts across targets.
- Record `chopper_input`, `chopper_output` and `cutadapt_total_processed` once
  per sample when all reporting targets agree. Warn on disagreements and
  retain target-specific observations; never silently choose one target.
  Preserve missing values and `NA` separately from zero.
- Keep primer-selection, Savont and accepted-assignment counts target-specific.
  Targets receive the same sample reads: do not sum shared preprocessing counts
  across targets or assume target branches partition the reads.
- Include read counts before and after M8 filtering of spurious sample-level
  ASV assignments. Do not introduce acceptance thresholds in this milestone.
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

## M10 - Summary plots

[ ]

**Goal:** `scripts/plot_summary.py` plots from the summary tables.

**Tasks:** agree the plots and the plotting dependency first, then implement.

**Review gate:** developer reviews the plots on real output.

## M11 - Pipeline parameter investigation

[ ]

**Goal:** investigate parameter sensitivity on small runs before full-scale
validation; distinguish changes in read recovery from changes in ASV inference.

**Tasks**
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
  preserve run provenance and do not tune solely to recover desired ASV counts.

**Review gate:** developer reviews parameter sensitivity and approves settings
for full-scale validation.

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
