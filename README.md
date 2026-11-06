# Amplicon pipeline

An Oxford Nanopore (ONT) amplicon sequencing pipeline. Given a raw run FASTQ
and a sample manifest, it will:

1. demultiplex reads by ONT and custom forward/reverse barcodes (Dorado),
2. preprocess each sample once with Porechop ABI and Chopper,
3. infer ASVs and assign taxonomy per target with
   [nf-core/ampliseq](https://nf-co.re/ampliseq) 3.0.0 (Savont),
4. summarise ASVs and per-sample outcomes into stable tables,
5. plot the results.

> **Status: in development.** Demultiplexing, preprocessing, per-target ampliseq
> ASV aggregation and sequence/assignment QC run. The pipeline then exits with an explicit error because
> sample summaries and plots are not yet implemented. See
> [MILESTONES.md](MILESTONES.md) for progress and the development plan.

## Requirements

- Python 3.12, PyYAML, [seqtk](https://github.com/lh3/seqtk), Porechop ABI
  0.5.0post1, Chopper 0.12.0b and [Nextflow](https://www.nextflow.io/) 26.04
  (provided by the project Conda environment; Nextflow is checked before the
  run starts)
- A container runtime (Docker, Singularity or Apptainer)
- [Dorado](https://github.com/nanoporetech/dorado) 0.8.2 (installed separately)
- A SINTAX reference database for each target
- Optionally an LSF cluster (`--profile lsf`)

## Installation

```bash
git clone https://github.com/<owner>/amplicon-pipeline.git
cd amplicon-pipeline
conda env create --file env.yaml
conda activate amplicon-pipeline
python run_pipeline.py --help
```

The single development and driver environment in [env.yaml](env.yaml) provides
Python 3.12, PyYAML, seqtk 1.4, Porechop ABI, Chopper and Nextflow with its Java
runtime. No separate test dependencies are needed. Update an existing
environment after dependency changes with
`conda env update --name amplicon-pipeline --file env.yaml`.

The project Conda environment does not supply Singularity/Apptainer (or Docker)
or Dorado 0.8.2; install these separately. Make the container runtime
available in the execution shell. Dorado is found through
`demux.dorado_exec` in your run configuration, or on `PATH` if unset; seqtk
always comes from the activated environment's `PATH`.

## Usage

```bash
python run_pipeline.py run.fastq.gz manifest.tsv config_dir output_dir --profile local
```

| Argument | Description |
|---|---|
| `run.fastq.gz` | Raw, undemultiplexed run FASTQ |
| `manifest.tsv` | Expected samples with their barcodes (schema not yet fixed) |
| `config_dir` | Run configuration directory (see below) |
| `output_dir` | Where results are written |
| `--profile` | `local` or `lsf` (default `local`) |

## Configuration

Copy [config/example/](config/example/) to a run-specific directory and edit:

- [pipeline.yaml](config/example/pipeline.yaml): shared `threads`, `demux`
  settings, optional shared `preprocess` tool arguments, native ampliseq
  parameters under `ampliseq.common` and `ampliseq.targets`, and QC settings
  under `summary.common`/`summary.targets`.
  Target values replace common values by parameter name, including entire
  option strings and lists; they are not concatenated or recursively merged.
  `ampliseq.targets` defines targets and execution order. The pipeline sets
  `input` and `outdir` last and writes flattened per-target `params.yaml`.
- [ampliseq.config](config/example/ampliseq.config): Nextflow configuration for
  your executor, container runtime and site resources.

Relative paths in the configuration resolve against the configuration
directory, so update the example barcode paths if you copy it elsewhere.
Separate target YAML files and the old top-level `targets` layout are no longer
supported. Summary settings are not passed to ampliseq.
Bundled barcode files are in [assets/barcodes/](assets/barcodes/); see their
[provenance notes](assets/barcodes/README.md).

## Read preprocessing

The driver runs preprocessing after demultiplexing and before preparing inputs
for all targets. It can also run standalone:

```bash
python scripts/preprocess_reads.py \
    --manifest manifest.tsv --config-dir config_dir \
    --demux-dir output_dir/demux --output-dir output_dir/preprocessed
```

Outputs retain the demultiplexed FASTQ filenames. Missing and empty expected
samples are listed in `excluded_samples.tsv`. The output directory must not
already exist; preprocessing does not overwrite prior outputs or resume a
partial run. Optional `preprocess.porechop_abi_options` and
`preprocess.chopper_options` in `pipeline.yaml` are passed as shell-style
argument strings; omitted or empty values use tool defaults. Both tools use
the shared `threads` value.
Standalone `--keep-intermediates` retains `preprocessed/work/`; otherwise it
is removed only after all samples succeed. Failed preprocessing retains work.
After all samples succeed, `preprocessed/read_counts.tsv` records
`sample_id`, `input_reads`, `porechop_reads` and `chopper_reads` in manifest
order, independently of retention mode. Missing inputs have blank counts;
empty inputs have input count 0 and blank downstream counts because tools
were not run. Measured empty tool outputs have count 0.

Input preparation requires `--preprocessed-dir`, not `--demux-dir`. It
copies upstream exclusion reasons and writes samplesheets pointing to the
remaining preprocessed FASTQs, leaving empty-read handling to ampliseq.
Both example configurations skip Porechop ABI and Chopper inside ampliseq.

## Outputs

```text
output_dir/
  .demux-completion.json    # demux input and output fingerprints
  demux/                    # barcode-named FASTQs and demux_stats.txt
  .preprocess-completion.json # preprocessing input and output fingerprints
  preprocessed/             # shared processed FASTQs and excluded_samples.tsv
  ampliseq/<target>/         # independent Nextflow launch directory
    samplesheet.tsv
    excluded_samples.tsv
    params.yaml
    resources.config
    .completion.json         # input/config fingerprints for safe target skip
    .nextflow/              # metadata; cache retained with --keep-intermediates
    .nextflow.log
    work/                   # retained with --keep-intermediates or on failure
    results/                # published ampliseq outputs
  summary/
    asvs.tsv                # unique target-ASVs and reported taxonomy
    assignments.tsv         # non-zero sample-target-ASV counts
    summarise_asvs.log       # successful aggregation paths and counts
  plots/                    # planned
```

The pooled CytB fixture run published these files under `results/`:

```text
savont/pooled_feature-table.tsv
savont/pooled_final_asvs.fasta
savont/pooled_stats.tsv
sintax/ASV_tax_sintax.user.tsv
sintax/ASV_tax_sintax.user.raw.tsv
summary_tables/ampliseq.taxonomy.SINTAX.user.tsv.gz
overall_summary.tsv
```

Counts, sequences and taxonomy share ASV IDs. These are observed pooled-mode,
custom-reference filenames, not yet a contract for other inference/database modes.
Future sample summaries must include every manifest sample, including those
without reads; input preparation currently records these in `excluded_samples.tsv`.

### ASV aggregation

The standalone aggregator reads one pooled Savont feature table and one
processed/raw SINTAX table pair per target:

```bash
python scripts/summarise_asvs.py \
    --target COI output_dir/ampliseq/COI/results \
    --target CytB output_dir/ampliseq/CytB/results \
    --output-dir output_dir/summary --config-dir config_dir
```

Repeat `--target TARGET RESULTS_DIR` for each target. The driver runs this
command after all configured ampliseq targets finish and includes it in
`--dry-run`. No external tools are needed to aggregate existing outputs.

- `asvs.tsv`: one row per target-ASV, with `target`, `asv_id`, `sequence`,
  `length`, the union of reported taxonomy-rank columns, `confidence`,
  `sintax_raw_taxonomy`, `sintax_strand` and `sintax_filtered_taxonomy`.
  Taxonomy strings and confidence values are unchanged; absent ranks are blank.
  Sequence QC adds `length_status`, `stop_free_frames`, `uncertain_frames`
  and `coding_status`. `total_reads` and `contam_target_cutoff` are stored
  once per target-ASV.
- `assignments.tsv`: one row per non-zero sample-target-ASV assignment, with
  `sample_id`, `target`, `asv_id`, integer `read_count`, `contam_possible`
  and `locov`. Both flags are explicit lowercase `true`/`false`.

Successful aggregation reports resolved target directories, input/output paths
and per-target unique ASVs, non-zero sample-ASV assignments and assigned reads
to stdout and `summarise_asvs.log`. The log is replaced on each successful run;
it is not an error log. Contamination diagnostics include target-ASV totals,
cutoffs and affected sample counts.

Rows follow supplied target and source-table order; rank columns follow first
appearance, including headers from empty targets. Empty targets contribute no
rows; both output files always have headers. Missing tables, duplicate IDs and
inconsistent IDs fail rather than silently dropping evidence. All inputs are
parsed and QC annotated before replacing the two canonical summary files.
No rows are filtered and taxonomy is not normalized; absent assignments are
not sample outcomes. Separate QC scripts/tables have been retired; historical
QC outputs are not removed or updated.

### Sequence QC

Sequence and assignment QC are part of aggregation, regenerating the two
canonical tables on every run.

`summary.common` and `summary.targets` merge by key, with target values winning.
`allowed_lengths` specifies discrete ASV lengths; `min_length`/`max_length`
specify inclusive bounds. If combined, all length checks must pass. Omitted
checks are `not_checked`. COI uses the approved discrete lengths; P3 uses
372-380 nt in the fixture to exercise a length failure, not as a production
recommendation. The fixture also checks CytB lengths 298, 301, 304, 307, 310, 313
and 316 nt with genetic code 5; ITS2 checks are deferred.
Length status is `too_short` or `too_long` outside the effective bounds,
`unexpected_length` for a gap in the allowed-length list, or `pass`.
Discrete lists supply their own minimum/maximum bounds; contradictory bounds
raise an error.

`genetic_code: 5` enables six-frame stop checking with the invertebrate
mitochondrial code (TAA/TAG stops; TGA is not a stop). Other codes are unsupported.
Shifts 0/1/2 are measured from each oriented sequence's start. Incomplete terminal
codons are ignored; ambiguous complete codons make a frame uncertain unless a
definite stop excludes it.

Sequence evidence appears only in `asvs.tsv`. Frame lists use
comma-separated strand/shift identifiers such as `+0,-2`. Coding status is
`not_checked`, `no_stop_free_frame`, `one_stop_free_frame`,
`multiple_stop_free_frames`, or `uncertain` if any candidate frame is uncertain.
All confirmed stop-free frames remain listed even when other frames are uncertain.
No ASVs are removed.

### Assignment QC

`total_reads` sums all assignments for a target-ASV, including sequence-QC
failures. `contam_target_cutoff` is `total_reads * max_contam_prop`.
`contam_possible` is true only when the total strictly exceeds
`min_contam_total_reads` and the assignment proportion is strictly below
`max_contam_prop`. Defaults are 10000 and 0.01. `locov` is independently true
when assignment reads are strictly below `max_locov` (default 100).
Target settings override common settings. Old contamination setting names
raise errors rather than being ignored.

Derived sequence passing requires length `pass`/`not_checked` and coding
`not_checked`/`one_stop_free_frame`/`multiple_stop_free_frames`. Assignment
passing requires both flags false; both-QC passing joins these conditions
on `(target, asv_id)`. No combined pass column or comma-separated QC labels
are stored.

### Restart handling

Demultiplexing is skipped only when its record matches the run, manifest,
configuration, barcode resources and all published files. An existing partial
or stale demux directory causes a clear error and is left untouched.
Preprocessing follows the same skip-or-refuse policy, with fingerprints of
demux files, the manifest, barcode mappings, preprocessing source/helpers,
effective options and threads, and all processed outputs. Summary or ampliseq
settings do not invalidate preprocessing. Standalone preprocessing does not
write a driver completion record.
Top-level demux/preprocessing `work/` is excluded from published-output
fingerprints and downstream input fingerprints.
For each target, ampliseq is skipped only when its input, configuration, selected profile
and published result files match; otherwise Nextflow runs again.
Pipeline `--keep-intermediates` retains stage work and Nextflow cache and
enables Nextflow `-resume`; without it Nextflow runs without `-resume`.
Each target uses its own `work/`. After success without the flag, only that
target's `work/` and `.nextflow/cache/` are removed, before completion recording.
Failures and skipped targets are not cleaned. Logs, results, other metadata
and external caches remain. Cleanup refuses symlinked scratch directories or
published links into scratch; use copy publication or retain intermediates.
Configuration fingerprints are stage-specific: demux uses its effective
settings/command and barcode resources; each ampliseq target uses its resolved
parameters, command, executor configuration and input/reference files. Unrelated
YAML changes do not invalidate upstream stages. Completion records from before
this fingerprint change are incompatible: existing demux outputs fail safely
without modification; use a fresh output directory for migration validation.
ASV aggregation always reruns, including when all ampliseq targets are skipped,
replacing its two summary tables from the current published results.

## Tests

```bash
python -m unittest discover -s tests -v
bash -n scripts/demux_hybrid.sh scripts/demux_illumina.sh
```

The tests check the scaffold interface, barcode consistency and the
[small mixed-target fixture](tests/fixtures/small/README.md).
They do not run Dorado or ampliseq.

## Contributing

Development conventions are in [AGENTS.md](AGENTS.md) and the plan in
[MILESTONES.md](MILESTONES.md).

## License

Released under the [MIT License](LICENSE).
