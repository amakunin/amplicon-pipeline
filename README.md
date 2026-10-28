# Amplicon pipeline

An Oxford Nanopore (ONT) amplicon sequencing pipeline. Given a raw run FASTQ
and a sample manifest, it will:

1. demultiplex reads by ONT and custom forward/reverse barcodes (Dorado),
2. infer ASVs and assign taxonomy per target with
   [nf-core/ampliseq](https://nf-co.re/ampliseq) 3.0.0 (Savont),
3. summarise ASVs and per-sample outcomes into stable tables,
4. plot the results.

> **Status: scaffold.** Only the command-line interface, configuration
> templates, barcode assets and a small test fixture exist. Running the
> pipeline exits with an explicit error and creates no outputs. See
> [MILESTONES.md](MILESTONES.md) for progress and the development plan.

## Requirements

Not all of these are used yet; they are what the finished pipeline will need.

- Python 3.12 (the scaffold uses the standard library only)
- [Nextflow](https://www.nextflow.io/) 26 or higher (required by ampliseq 3.0.0; checked before the run starts) and a container runtime (Docker,
  Singularity or Apptainer)
- [Dorado](https://github.com/nanoporetech/dorado) 0.8.2 and
  [seqtk](https://github.com/lh3/seqtk)
- A SINTAX reference database for each target (for COI, MIDORI2 or the
  ampliseq built-in `coidb`)
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
Python 3.12 and seqtk. Tests use the standard library, so no separate test
dependencies are needed. Update an existing environment after dependency changes
with `conda env update --name amplicon-pipeline --file env.yaml`.

Nextflow, its Java runtime, Singularity/Apptainer (or Docker), and Dorado 0.8.2
remain external tools supplied by your installation. Make the launcher and
container runtime available in the execution shell. Dorado is found through
`demux.dorado_exec` in your run configuration, or on `PATH` if unset; seqtk
always comes from the activated environment's `PATH`. The environment alone
does not make the scaffold an executable analysis pipeline.

### Using the environment with a coding agent

Selecting the environment with VS Code's **Python: Select Interpreter** command
configures Python editor tooling, but does not guarantee activation for agent
shell commands. Activating it in a separate terminal is also not sufficient.
After creating it, obtain its prefix and your Conda executable path:

```bash
conda env list
conda info --base
```

Give the agent the environment prefix and the absolute Conda executable path
(normally `bin/conda` under the Conda base directory). The agent can then run
each command explicitly, without relying on activation persisting between
shell calls:

```bash
/path/to/conda/bin/conda run --no-capture-output -p /path/to/env \
    python -m unittest discover -s tests -v
```

Run commands from the pipeline root. Keep machine-specific executable paths and
environment prefixes out of committed configuration and documentation.

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

- [pipeline.yaml](config/example/pipeline.yaml): threads (demultiplexing and ampliseq CPU limit), targets to run, demultiplexing
  method, optional Dorado path and barcode resources.
- [targets/COI.yaml](config/example/targets/COI.yaml): one file per target, using
  native ampliseq parameters (primers, reference database, Savont options).
  The pipeline sets `input` and `outdir` itself.
- [ampliseq.config](config/example/ampliseq.config): Nextflow configuration for
  your executor, container runtime and site resources.

Relative paths in the configuration resolve against the configuration
directory, so update the example barcode paths if you copy it elsewhere.
Bundled barcode files are in [assets/barcodes/](assets/barcodes/); see their
[provenance notes](assets/barcodes/README.md).

## Outputs

```text
output_dir/
  demux/                    # barcode-named FASTQs and demux_stats.txt
  ampliseq/<target>/         # independent Nextflow launch directory
    samplesheet.tsv
    excluded_samples.tsv
    params.yaml
    resources.config
    .nextflow/              # retained cache
    .nextflow.log
    work/                   # retained task files
    results/                # published ampliseq outputs
  summary/                  # planned ASV, count and sample/target QC tables
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
