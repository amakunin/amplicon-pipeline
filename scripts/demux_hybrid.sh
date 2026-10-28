#!/usr/bin/env bash
# Hybrid demultiplexing: ONT kit barcodes (Dorado) followed by Sanger R
# barcodes at either read end. Requires dorado and seqtk.
set -euo pipefail
shopt -s nullglob

usage() {
    echo "Usage: $0 --input FASTQ --output-dir DIR --barcode-dir DIR [--dorado PATH] [--threads N]" >&2
    exit 1
}

input="" outdir="" bcdir="" dorado="dorado" threads=16
while [ $# -gt 0 ]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output-dir) outdir="$2"; shift 2 ;;
        --barcode-dir) bcdir="$2"; shift 2 ;;
        --dorado) dorado="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *) usage ;;
    esac
done
[ -n "$input" ] && [ -n "$outdir" ] && [ -n "$bcdir" ] || usage

for tool in "$dorado" seqtk; do
    command -v "$tool" > /dev/null || { echo "Required tool not found: $tool" >&2; exit 1; }
done
[ -f "$input" ] || { echo "Input FASTQ not found: $input" >&2; exit 1; }
for ext in tmol fasta; do
    [ -f "$bcdir/Sanger_custom_barcode_R.$ext" ] \
        || { echo "Barcode file not found: $bcdir/Sanger_custom_barcode_R.$ext" >&2; exit 1; }
done
[ ! -e "$outdir" ] || { echo "Output directory exists, terminating: $outdir" >&2; exit 1; }

work="$outdir/.work"
mkdir -p "$work"

# Sanger R demux of one FASTQ into a directory
demux_r() {
    "$dorado" demux \
        --emit-fastq \
        --barcode-arrangement "$bcdir/Sanger_custom_barcode_R.tmol" \
        --kit-name SciOps \
        --barcode-sequences "$bcdir/Sanger_custom_barcode_R.fasta" \
        --threads "$threads" \
        --output-dir "$2" \
        "$1"
}

###############################################################################
# Stage 1: ONT kit demultiplex
###############################################################################

# --no-trim mimics Illumina-only demultiplexing
"$dorado" demux \
    --emit-fastq \
    --kit-name SQK-NBD114-24 \
    --no-trim \
    --threads "$threads" \
    --output-dir "$work/1_ont" \
    "$input"

###############################################################################
# Stage 2: Sanger R barcode at read start (rr), then at read end (rf)
###############################################################################

for file in "$work"/1_ont/*barcode*.fastq; do
    [[ "$file" != *unclassified* ]] || continue
    barcode=${file##*_}
    barcode=${barcode%.fastq}
    rrdir="$work/2_rr/$barcode"

    demux_r "$file" "$rrdir"

    # reads without a start barcode: reverse complement and look again
    for unclass in "$rrdir"/*unclassified.fastq; do
        seqtk seq -r "$unclass" > "$rrdir/unclassified_rc.fastq"
        demux_r "$rrdir/unclassified_rc.fastq" "$work/3_rf/$barcode"
        rm "$unclass" "$rrdir/unclassified_rc.fastq"
    done
done

###############################################################################
# Stage 3: merge rr and rf into <ont>_<R>.fastq, count, compress
###############################################################################

for file in "$work"/2_rr/*/*.fastq "$work"/3_rf/*/*.fastq; do
    ont_barcode=$(basename "$(dirname "$file")")
    r_barcode=${file##*_}
    cat "$file" >> "$outdir/${ont_barcode}_${r_barcode}"
done

: > "$outdir/unclassified.fastq"
for file in "$work"/1_ont/*unclassified.fastq; do
    cat "$file" >> "$outdir/unclassified.fastq"
done

: > "$outdir/demux_stats.txt"
for file in "$outdir"/*.fastq; do
    printf '%s\t%s\n' "$(basename "$file" .fastq)" "$(( $(wc -l < "$file") / 4 ))" \
        >> "$outdir/demux_stats.txt"
    gzip "$file"
done

rm -r "$work"
