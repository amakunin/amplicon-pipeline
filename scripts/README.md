# Demultiplexing scripts

`demux_hybrid.sh` (ONT kit barcode, then Sanger R) and `demux_illumina.sh`
(Sanger F, then Sanger R) share one interface and output contract. Both are
runnable standalone.

```
demux_*.sh --input FASTQ --output-dir DIR --barcode-dir DIR \
           [--dorado PATH] [--threads N] [--keep-intermediates]
```

- `dorado` (`--dorado`, default from `PATH`) and `seqtk` (from `PATH`, e.g. a
  conda environment) must be available.
- `--output-dir` must not exist. Nothing is deleted automatically.
- Exit 0 on success. Bad input or missing tools give exit 1 and a one-line
  message; a failed Dorado run stops the script.

## Output

| File | Content |
|---|---|
| `ont_barcodeNN_barcodeMM.fastq.gz` | hybrid: ONT barcode NN + Sanger R barcode MM |
| `barcodeNN_barcodeMM.fastq.gz` | Illumina: Sanger F barcode NN + Sanger R barcode MM |
| `ont_barcodeNN_unclassified.fastq.gz`, `barcodeNN_unclassified.fastq.gz` | first-level barcode found, no R barcode |
| `unclassified.fastq.gz` | no first-level barcode |
| `demux_stats.txt` | `name<TAB>read count`, sorted by name |

Barcode pairs without reads have no file and no stats line. `.work/` holds
intermediates and is kept only with `--keep-intermediates` or after a failure.

## Method-specific

Barcode files used: hybrid uses the Dorado kit `SQK-NBD114-24` and
`Sanger_custom_barcode_R.*`; Illumina uses `Sanger_custom_barcode_{F,R}.*`.
Mapping manifest barcodes to these file names belongs to the driver (M4d).
