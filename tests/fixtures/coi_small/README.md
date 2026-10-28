# Small COI integration fixture

## Inputs and expectations

Public redistribution of the fixture read data and read ID list has been
approved.

- [run.fastq.gz](run.fastq.gz): user-generated raw reads, restricted to COI by
  motif-based ID selection and read recovery. No demultiplexing has yet been
  performed by this scaffold.
- [read_ids.txt](read_ids.txt): the user-confirmed final COI selection, moved
  from `tests/ids.txt` without changing its contents.
- [manifest.tsv](manifest.tsv): five expected samples. Whitespace separators
  were normalised to tabs without changing the values.
- [expected_demux.tsv](expected_demux.tsv): user-confirmed expectations for
  demultiplexed read counts, not assertions about ASV inference success.

| Sample | ONT barcode | Rear barcode | Expected demux reads |
|---|---|---|---|
| ANLCx_001_A1 | 09 | 01 | > 0 |
| ANLCx_001_B1 | 09 | 02 | 0 |
| ANLCx_001_C1 | 09 | 03 | > 0 |
| ANLCx_001_D1 | 09 | 04 | 0 |
| ANLCx_001_E1 | 09 | 05 | > 0 |

Empty samples must remain visible in downstream sample summaries, even if they
cannot be submitted to ampliseq as nonempty FASTQs. Expected taxa and exact
read/ASV counts are not yet fixed. This selected COI dataset does not establish
target-routing performance on a mixed-target run.

## Fixture preparation

From the pipeline root, with the original input paths supplied:

```bash
# INPUT_FOR_ID_SELECTION may be a selected sample subset.
# RAW_RUN_FASTQ is the pre-demultiplexing source from which reads are recovered.
bash tests/helpers/extract_read_ids.sh "$INPUT_FOR_ID_SELECTION" > selected_ids.txt
bash tests/helpers/retrieve_reads.sh selected_ids.txt "$RAW_RUN_FASTQ" \
    | gzip -c > new_run.fastq.gz
gzip -t new_run.fastq.gz
```

In a shell script, use `set -euo pipefail` so helper errors cannot be hidden by
gzip. Use new output paths and validate before replacing an existing fixture.
Do not regenerate the current fixture just to test these commands.

The first helper emits all read IDs by default, sorted and deduplicated.
For optional sequence filtering, supply `--motifs` with a case-insensitive
AWK extended regular expression. The original COI selection can be requested
explicitly:

```bash
bash tests/helpers/extract_read_ids.sh \
    --motifs 'cataaagatattgg|gtgaccaaaaaatca' "$INPUT_FOR_ID_SELECTION" \
    > selected_ids.txt
```

This motif filter is a fixture-selection heuristic, not production primer
assignment or a complete test for COI specificity. For mixed-target selection,
omit the filter to retain all reads from the chosen sample subset.

The second helper recovers complete four-line FASTQ records by ID, preserving
source order and original headers. Both helpers accept plain or gzipped FASTQ.
The retrieval helper reports requested/matched/recovered/missing counts to
stderr and fails for an empty list or missing IDs. Failed runs may leave
partial redirected output; discard it after inspecting the error.

The owner confirmed use of these two helpers and the final ID list. The exact
original source FASTQ paths and any upstream sample-subsetting command have
not yet been recorded; do not infer them from the old pipeline's run config.

## Execution configuration

[config/pipeline.yaml](config/pipeline.yaml) uses the approved Dorado 0.8.2
binary, seqtk wrapper and bundled barcode assets. Relative resource paths are
relative to this fixture's configuration directory.

[config/targets/COI.yaml](config/targets/COI.yaml) retains pooled Savont
inference and the prototype primers/length settings. Before an analysis run,
set `sintax_ref` to an accessible local MIDORI2 COI SINTAX reference; the
reference is external and approximately 1.6 GB, so a tiny read input does not
imply tiny taxonomy memory requirements.

[config/ampliseq.config](config/ampliseq.config) defines profile `lsf` with
Singularity, queue `long`, queueSize 4, Savont 2 CPUs/8 GB, and SINTAX 10 GB.
These are unmeasured starting allocations, not guaranteed resource
requirements. Configure any site-specific bind mounts, queue settings and LSF
accounting options in a local copy. No job arrays, automatic retries or remote
weblog endpoint are enabled.

The pipeline driver and analysis stages are still placeholders. Once demux and
input preparation exist, launch ampliseq from a stable target execution directory
using this config and generated absolute input/output parameters. Do not launch
the raw run FASTQ directly into ampliseq.

## Validation

The initial fixture contains 1,073 records with unique IDs, all matching the
selected ID list. Input checksums are recorded in
[checksums.sha256](checksums.sha256); verify from this directory using
`sha256sum -c checksums.sha256`.

From the pipeline root:

```bash
python -m unittest discover -s tests -v
gzip -t tests/fixtures/coi_small/run.fastq.gz
```

The tests check expected barcode mapping, helper behaviour and, if present,
real FASTQ structure/ID recovery. Actual positive/zero demux read counts remain
unverified until Dorado runs.
