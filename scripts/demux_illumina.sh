#!/usr/bin/env bash
# Illumina-style demultiplexing: Sanger F barcodes then Sanger R barcodes
# (Dorado). Requires dorado and seqtk.
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
for end in F R; do
    for ext in tmol fasta; do
        [ -f "$bcdir/Sanger_custom_barcode_$end.$ext" ] \
            || { echo "Barcode file not found: $bcdir/Sanger_custom_barcode_$end.$ext" >&2; exit 1; }
    done
done
[ ! -e "$outdir" ] || { echo "Output directory exists, terminating: $outdir" >&2; exit 1; }

work="$outdir/.work"
mkdir -p "$work"

# Sanger demux of one FASTQ into a directory; $3 is F or R
demux_sanger() {
    "$dorado" demux \
        --emit-fastq \
        --barcode-arrangement "$bcdir/Sanger_custom_barcode_$3.tmol" \
        --kit-name SciOps \
        --barcode-sequences "$bcdir/Sanger_custom_barcode_$3.fasta" \
        --threads "$threads" \
        --output-dir "$2" \
        "$1"
}

###############################################################################
# Stage 1: F barcodes (i5, fewer barcodes so faster than R first)
###############################################################################

# ff: F barcode found at read end
demux_sanger "$input" "$work/1_ff" F

# fr: reverse complement the rest, so that a barcode that was at read start
# is now at the end
for unclass in "$work"/1_ff/*unclassified.fastq; do
    seqtk seq -r "$unclass" > "$work/unclassified_rc.fastq"
done
demux_sanger "$work/unclassified_rc.fastq" "$work/1_fr" F

###############################################################################
# Stage 2: R barcodes (i7, 96) for each F barcode separately
###############################################################################

for fbarcode in $(ls "$work"/1_ff/*barcode*.fastq "$work"/1_fr/*barcode*.fastq \
                    | sed 's/.*_\(barcode[0-9]*\)\.fastq$/\1/' | sort -u); do
    mkdir -p "$work/2_r/$fbarcode"
    cat "$work"/1_f[fr]/*_"$fbarcode".fastq > "$work/${fbarcode}_concat.fastq"
    demux_sanger "$work/${fbarcode}_concat.fastq" "$work/2_r/$fbarcode" R
done

###############################################################################
# Stage 3: merge to <F>_<R>.fastq, count, compress
###############################################################################

# several files per barcode pair (read groups) are concatenated
for file in "$work"/2_r/*/*.fastq; do
    f_barcode=$(basename "$(dirname "$file")")
    r_barcode=${file##*_}
    cat "$file" >> "$outdir/${f_barcode}_${r_barcode}"
done

: > "$outdir/unclassified.fastq"
for file in "$work"/1_fr/*unclassified.fastq; do
    cat "$file" >> "$outdir/unclassified.fastq"
done

: > "$outdir/demux_stats.txt"
for file in "$outdir"/*.fastq; do
    printf '%s\t%s\n' "$(basename "$file" .fastq)" "$(( $(wc -l < "$file") / 4 ))" \
        >> "$outdir/demux_stats.txt"
    gzip "$file"
done

rm -r "$work"
