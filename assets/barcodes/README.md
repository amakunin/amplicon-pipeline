# Barcode assets

These files were imported from an existing project prototype with approval for
use in this repository. Public redistribution has been approved for this
repository.

- `Sanger_custom_barcode_F.{fasta,tmol,tsv}`: 16 forward 8 bp Illumina-scheme
  barcodes, Dorado arrangement and sequence/index/plate mapping.
- `Sanger_custom_barcode_R.{fasta,tmol,tsv}`: 96 rear 8 bp Illumina-scheme
  barcodes, arrangement and sequence/index/well mapping.
- [ONT_nb24_barcode.tsv](ONT_nb24_barcode.tsv): standard ONT barcode01-24
  sequence/index/plate mapping used with `SQK-NBD114-24`.

The hybrid method uses Dorado's built-in ONT kit at the front and the custom R
arrangement for its second demultiplexing round. It does **not** use the custom
F arrangement. The Illumina-style method uses the custom F and R resources;
that method is not implemented in the new scaffold.

The manifest's `ont_barcode` column, not `front_barcode`, identifies the hybrid
forward barcode. Plate labels describe the source mapping and must not override
the sample identity or well recorded in the manifest.

Sequences, mappings and scoring settings are preserved from the source.
Text files have LF line endings and a final newline. Future barcode changes should be
reviewed with mapping/sequence consistency tests.
