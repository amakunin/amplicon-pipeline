#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'Usage: %s <fastq[.gz]>\n' "$0" >&2
    exit 2
fi

fastq="$1"
[[ -r "$fastq" ]] || { printf 'Cannot read FASTQ: %s\n' "$fastq" >&2; exit 1; }
motifs='cataaagatattgg|gtgaccaaaaaatca'

zcat -f "$fastq" |
    awk -v motifs="$motifs" '
        FNR % 4 == 1 { header = $0 }
        FNR % 4 == 2 {
            if (tolower($0) ~ motifs) {
                sub(/^@/, "", header)
                sub(/[[:space:]].*$/, "", header)
                print header
            }
        }
    ' |
    sort -u
