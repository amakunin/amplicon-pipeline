#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    printf 'Usage: %s <ids.txt> <fastq[.gz]>\n' "$0" >&2
    exit 2
fi

ids_file="$1"
fastq="$2"
[[ -r "$ids_file" ]] || { printf 'Cannot read ID list: %s\n' "$ids_file" >&2; exit 1; }
[[ -s "$ids_file" ]] || { printf 'ID list is empty: %s\n' "$ids_file" >&2; exit 1; }
[[ -r "$fastq" ]] || { printf 'Cannot read FASTQ: %s\n' "$fastq" >&2; exit 1; }

zcat -f "$fastq" |
    awk '
        FILENAME != "-" {
            if (NF == 0) next
            ids[$1] = 1
            next
        }

        FNR % 4 == 1 { header = $0 }
        FNR % 4 == 2 { seq = $0 }
        FNR % 4 == 3 { plus = $0 }
        FNR % 4 == 0 {
            qual = $0
            id = header
            sub(/^@/, "", id)
            sub(/[[:space:]].*$/, "", id)

            if (id in ids) {
                print header
                print seq
                print plus
                print qual
                found[id] = 1
                matches++
            }
        }

        END {
            for (id in ids)
                requested++

            for (id in found)
                recovered++

            printf "Requested IDs: %d\nMatching records: %d\nUnique IDs recovered: %d\nMissing IDs: %d\n",
                   requested, matches, recovered, requested - recovered > "/dev/stderr"
            if (requested == 0 || requested != recovered) exit 1
        }
    ' "$ids_file" -
