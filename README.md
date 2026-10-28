# Amplicon pipeline

Scaffold for an ONT amplicon sequencing pipeline combining barcode
demultiplexing, nf-core/ampliseq, ASV summaries, sample annotations and plots.

**Status:** only the entrypoint interface and configuration templates exist.
Execution exits with an explicit error and creates no outputs. Stage scripts
are placeholders, not working analysis tools.

## Layout

- [run_pipeline.py](run_pipeline.py): fixed-order orchestration entrypoint.
- [scripts/demux_hybrid.sh](scripts/demux_hybrid.sh): planned adaptation of the
  existing Dorado hybrid demultiplexing logic.
- [scripts/prepare_inputs.py](scripts/prepare_inputs.py): barcode/sample mapping
  and target-specific ampliseq samplesheets.
- [scripts/summarise_asvs.py](scripts/summarise_asvs.py): Savont counts, sequences
  and taxonomy adapted to stable tables.
- [scripts/summarise_samples.py](scripts/summarise_samples.py): manifest-driven
  sample annotations and explicit QC outcomes.
- [scripts/plot_summary.py](scripts/plot_summary.py): plots from summary tables.
- [config/example/](config/example/): configuration to copy and customise.
- [assets/barcodes/](assets/barcodes/): bundled custom F/R arrangements,
  sequences and barcode mappings; see [asset provenance](assets/barcodes/README.md).
- [tests/fixtures/coi_small/](tests/fixtures/coi_small/): first COI fixture,
  expected demux outcomes and small-test LSF settings.
- [tests/helpers/](tests/helpers/): read-subsetting helpers, not production
  pipeline stages. See [fixture preparation](tests/fixtures/coi_small/README.md).
- [tests/test_scaffold.py](tests/test_scaffold.py): entrypoint and placeholder
  smoke tests using the Python standard library.

No module registry, configurable DAG or plugin framework is planned. Ampliseq
will be called directly by the driver, once per target, pinned to version 3.0.0.

## Intended interface

```bash
python run_pipeline.py --help
python run_pipeline.py run.fastq.gz manifest.tsv config_dir output_dir --profile local
```

The second command currently fails deliberately. Python 3.12 is the initial
development environment; no third-party Python dependencies are required by the
scaffold.

Before execution is implemented:

1. Copy [config/example/](config/example/) to a run-specific configuration directory.
2. Set barcode resource paths in [pipeline.yaml](config/example/pipeline.yaml).
3. Configure target primers and reference databases in
   [targets/COI.yaml](config/example/targets/COI.yaml).
4. Configure the executor, container runtime and site resources in
   [ampliseq.config](config/example/ampliseq.config).

`--profile local|lsf` selects the intended execution environment. Scheduler
submission, dependency checks, configuration loading, dry-run and resume support
are not implemented yet.

## Planned configuration contract

- Relative resource paths resolve against the user-supplied configuration directory.
  Update example barcode paths if copying the configuration elsewhere.
- `pipeline.yaml` lists targets explicitly and contains demultiplexing settings.
- Target YAML files use native ampliseq parameters. The driver owns `input`
  and `outdir` and will add them to generated, resolved parameter files.
- The driver will pass `ampliseq.config` to Nextflow with `-c`, and select the
  matching `local` or `lsf` profile. Choose a container profile explicitly in
  that file when configuring the installation.
- The COI template retains the prototype's permissive empty/filtering/trimming
  settings. Missing samples must be reported, not silently treated as success.
- Reporting thresholds remain unset until their scientific definitions and
  denominators are agreed.

## Planned data flow and outputs

```text
validate -> demultiplex -> prepare inputs -> ampliseq per target
         -> summarise ASVs -> summarise samples -> plot
```

```text
output_dir/
  demux/              # FASTQs, barcode mapping and read counts
  inputs/             # target samplesheets and resolved parameters
  ampliseq/<target>/  # published ampliseq outputs
  summary/            # ASV, count and sample/target QC tables
  plots/
  execution/<target>/ # stable Nextflow launch directory, cache, logs and work/
```

The manifest schema and target-routing policy are intentionally not fixed yet.
First verify whether each demultiplexed FASTQ contains one or several targets
and test primer-based selection before connecting it to ampliseq.

The ASV adapter must use Savont outputs: ampliseq 3.0.0's long-format counts
summary is DADA2-only. Verify actual published pooled-Savont outputs using a
small real run before fixing adapter filenames.

Sample summaries will start from all expected manifest sample/target pairs.
They must distinguish no reads, preprocessing loss, no ASVs, unresolved
taxonomy, QC rejection and processing failure. Preserve raw counts and
classifications alongside acceptance flags and reasons.

Restart handling will validate completion records and relevant inputs/config,
not merely test for output-file existence. Each ampliseq run will retain both
its Nextflow cache and work directory for `-resume`. Never automatically
remove previous results on a failed or changed run.

## Next implementation step

Adapt hybrid demultiplexing without hard-coded site paths or unconditional
environment-module loading, then run one target end to end on a small fixture.
Do not add summary rules or orchestration abstractions before verifying that
vertical slice.

## Scaffold tests

From this directory:

```bash
python -m unittest discover -s tests -v
bash -n scripts/demux_hybrid.sh
```

Tests cover the scaffold interface, barcode/mapping consistency, fixture
expectations and read-subsetting helpers. If the real fixture FASTQ is present,
they also validate its structure and selected read IDs. They do not execute
Dorado/ampliseq or assert biological accuracy.
