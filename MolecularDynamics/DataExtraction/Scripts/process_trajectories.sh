#!/usr/bin/env bash
set -euo pipefail

IN_DIR="../Raw_Data"
OUT_DIR="../Processed_Data"

mkdir -p "$OUT_DIR"

for xtc in "$IN_DIR"/*.xtc; do
    [[ -e "$xtc" ]] || continue
    base=$(basename "$xtc" .xtc)

    tpr="$IN_DIR/${base}.tpr"
    out="$OUT_DIR/${base}_Proc.xtc"

    if [[ ! -f "$out" ]]; then
        if [[ -f "$tpr" ]]; then
            echo "Processing $base ..."
            echo -e "DNA\nDNA" | gmx trjconv \
                -s "$tpr" \
                -f "$xtc" \
                -o "$out" \
                -pbc nojump -center
        else
            echo "Skipping $base (no matching .tpr found)"
        fi
    else
        echo "Skipping $base (already processed)"
    fi
done