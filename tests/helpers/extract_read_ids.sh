#!/usr/bin/env bash
set -euo pipefail

motifs=''
if [[ $# -eq 3 && "$1" == "--motifs" && -n "$2" ]]; then
    motifs="$2"
    fastq="$3"
elif [[ $# -eq 1 && "$1" != --* ]]; then
    fastq="$1"
else
    printf 'Usage: %s [--motifs <regex>] <fastq[.gz]>\n' "$0" >&2
    exit 2
fi

[[ -r "$fastq" ]] || { printf 'Cannot read FASTQ: %s\n' "$fastq" >&2; exit 1; }

zcat -f "$fastq" |
    awk -v motifs="$motifs" '
        BEGIN { motifs = tolower(motifs) }
        FNR % 4 == 1 { header = $0 }
        FNR % 4 == 2 {
            if (motifs == "" || tolower($0) ~ motifs) {
                sub(/^@/, "", header)
                sub(/[[:space:]].*$/, "", header)
                print header
            }
        }
    ' |
    sort -u
