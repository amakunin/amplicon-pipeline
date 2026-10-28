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

Both hybrid and Illumina demultiplexing are intended to use this dataset;
actual outcomes remain to be verified in M4. Empty samples must remain visible
in downstream summaries even when they cannot be submitted to ampliseq.

The previous pipeline reports 583 usable reads across the three samples.
These counts are sanity-check references, not exact acceptance thresholds:
they describe a later processing stage than the 931 recovered raw reads.
Exact target-level read counts, ASVs and taxonomy from this pipeline remain
unverified until analysis runs.

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

The existing [configuration](config/pipeline.yaml) currently covers hybrid
demultiplexing and [COI](config/targets/COI.yaml) only. Both-method
demultiplexing configuration is established in M4; ITS2, CytB and P3
configuration is added in M6. Site executable paths, executor settings and
external reference databases must be configured for the execution environment.
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
