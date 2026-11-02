# Small mixed-target integration fixture

## Provenance

The developer supplied 931 raw reads from source run `00040` (an identifier
internal to the previous pipeline). Three samples were selected for low read
counts while recovering all four amplicons: COI, ITS2, CytB and P3.
Public redistribution of the updated reads, manifest and reference counts
has been approved by the developer.

- [run.fastq.gz](run.fastq.gz): recovered raw reads for the three selected
  samples, retaining the original reads for demultiplexing.
- [manifest.tsv](manifest.tsv): the three selected samples and two additional
  samples whose barcodes are absent from this dataset.
- [expected_demux.tsv](expected_demux.tsv): expected positive/zero sample
  outcomes and ONT/rear barcode indices, not measured demultiplexing counts.
- [old_pipeline_read_counts.tsv](old_pipeline_read_counts.tsv): reads usable
  in sample-level taxonomic assignment in the previous pipeline, roughly
  corresponding to this pipeline's M8.

The extraction ID list was deliberately not retained. The source subset
command and use of any motif filter were not recorded here; do not infer them.
Integrity checks use the delivered FASTQ and checksums, not a retained ID list.

## Expectations

| Sample | ONT index | Front index | Rear index | Expected demux reads |
|---|---|---|---|---|
| ANLCx_002_A12 | 10 | 02 | 89 | > 0 |
| BFHK_023_H8 | 16 | 08 | 64 | > 0 |
| BFHK_023_A2 | 16 | 08 | 09 | > 0 |
| ANLCx_001_A1 | 09 | 01 | 01 | 0 |
| ANLCx_001_B1 | 09 | 01 | 02 | 0 |

Both hybrid and Illumina demultiplexing are intended to use this dataset.
The captured multi-target run uses hybrid demultiplexing. Empty samples must
remain visible in downstream summaries even when they cannot be submitted to
ampliseq.

The previous pipeline reports 583 usable reads across the three samples.
These counts are sanity-check references, not exact acceptance thresholds:
they describe a later processing stage than the 931 recovered raw reads.
Observed target-level outputs are recorded below; taxonomy is not independently
verified.

## Preparation helpers

From the pipeline root, with source paths supplied:

```bash
bash tests/helpers/extract_read_ids.sh "$INPUT_FOR_ID_SELECTION" > selected_ids.txt
bash tests/helpers/retrieve_reads.sh selected_ids.txt "$RAW_RUN_FASTQ" \
    | gzip -c > new_run.fastq.gz
gzip -t new_run.fastq.gz
```

Use `set -euo pipefail` in shell scripts so failures cannot be hidden by gzip.
Use new output paths; do not regenerate or overwrite the fixture without
developer approval.

ID extraction emits all read IDs by default, sorted and deduplicated.
Optional `--motifs '<regex>'` filters sequence lines with a case-insensitive
AWK extended regular expression. Retrieval recovers complete four-line FASTQ
records in source order with original headers, and fails for empty or missing
IDs. Both helpers accept plain or gzipped FASTQ.

## Configuration and validation

The [configuration](config/pipeline.yaml) covers hybrid demultiplexing and
COI, CytB, ITS2 and P3, run sequentially in that order. Site executable paths,
executor settings and external reference databases must be configured for the execution environment.
The raw run FASTQ must never be passed directly to ampliseq.

The delivered FASTQ was verified to contain 931 complete records with unique
read IDs, nonempty sequences and matching sequence/quality lengths.
From this directory, verify input integrity with:

```bash
sha256sum -c checksums.sha256
```

From the pipeline root, run `python -m unittest discover -s tests -v`.
Tests check input checksums, FASTQ structure and record count, barcode
expectations, configuration/manifest loading and the dry-run plan without
creating outputs. They do not run demultiplexing or ampliseq.

## Captured multi-target outputs

The developer ran the following on 2026-11-02 with the default `local` profile,
then reran with corrected P3 taxonomic levels:

```bash
python run_pipeline.py tests/fixtures/small/run.fastq.gz \
    tests/fixtures/small/manifest.tsv tests/fixtures/small/config results/m6b_test
```

[ampliseq/](ampliseq/) contains four byte-preserved files per target, copied
from published outputs into `<target>/` as frozen inputs for downstream tests:

- `overall_summary.tsv`: per-sample read counts through preprocessing and ASV inference.
- `savont/pooled_feature-table.tsv`: ASV counts per sample.
- `sintax/ASV_tax_sintax.<database>.tsv`: consensus taxonomy and sequence per ASV ID.
- `sintax/ASV_tax_sintax.<database>.raw.tsv`: unparsed SINTAX calls, including
  confidences useful below the 0.8 cutoff.

`<database>` is `coidb` for COI and `user` for the other targets. No other run
outputs are captured. The published `results/` directory level is omitted;
large run outputs remain outside the fixture. These are regression examples, not biological acceptance
thresholds. Remaining cross-target taxonomy-reporting inconsistencies are
deferred to M8; P3 values are preserved without normalization.

| Target | ASVs | ASV ID (total reads) |
|---|---:|---|
| COI | 0 | none |
| CytB | 2 | `0aece8f365778d45` (308), `9e07fe1354e353c7` (29) |
| ITS2 | 1 | `4851853f76e4a41b` (59) |
| P3 | 2 | `0a4b74ea4eddec52` (18), `83ea72fd937a3d48` (29) |

COI has header-only feature/processed taxonomy tables and an empty raw taxonomy
file. Its overall summary reports `NA`, not zero, for `savont_output`.
The run's excluded-sample lists reported `ANLCx_001_A1` and `ANLCx_001_B1`
as `no_demux_file`; those lists are not captured.

Verify captured bytes from [ampliseq/](ampliseq/) with
`sha256sum -c checksums.sha256`. Unit tests also verify checksums, ASV counts,
matching unique IDs and nonempty sequences in nonempty targets.

### Historical M5c evidence

The local/Singularity CytB run on 2026-10-28 completed 26 tasks with no failed
or ignored tasks. Its counts agree with the captured multi-target CytB output:

| Sample | Chopper input | Chopper output | Primer-trimmed | Savont counts |
|---|---:|---:|---:|---:|
| ANLCx_002_A12 | 202 | 184 | 38 | 29 |
| BFHK_023_A2 | 398 | 356 | 187 | 177 |
| BFHK_023_H8 | 309 | 269 | 137 | 131 |

ASV `0aece8f365778d45` was assigned Homo sapiens (confidence 0.90);
`9e07fe1354e353c7` was assigned Bos indicus (0.97). Both had
`ampliseq_accept=TRUE`; assignments were not independently verified. IDs agreed
across counts, sequences and taxonomy. Two absent samples remained excluded.
Savont enabled low-polymorphism mode after an 88.51% no-SNPmer warning.

An earlier COI run completed with zero ASVs, an empty FASTA, header-only counts
and no published taxonomy; its export log contained an empty-FASTA error despite
workflow success. In contrast, the captured multi-target run publishes empty
taxonomy tables. Execution success must not be equated with ASV recovery.
No thresholds or analysis parameters were changed during the M5c inspection.
