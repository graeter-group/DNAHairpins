#!/bin/bash
set -euo pipefail
module load gromacs/2025.2

DATA_DIRECTORY_TPRS="../RawData/"
DATA_DIRECTORY_PROCESSED="../ProcessedData/"
INDEX_DIRECTORY="../IndexFiles/"
OUTPUT_DIRECTORY="../../Analysis/ExtractedData/"

mkdir -p "$OUTPUT_DIRECTORY"

for tprfile in "$DATA_DIRECTORY_TPRS"/*.tpr; do
    [ -e "$tprfile" ] || continue

    filename=$(basename "$tprfile" .tpr)

    if [[ "$filename" == HP060* ]]; then
        echo "Skipping $filename — ignored by rule."
        continue
    fi

    # expected inputs
    trajfile="${DATA_DIRECTORY_PROCESSED}/${filename}_Proc.xtc"
    ndxfile="${INDEX_DIRECTORY}/${filename}_BackboneBonds.ndx"
    outfile="${OUTPUT_DIRECTORY}/${filename}_BackboneBonds.txt"

    # skip if already processed
    if [ -f "$outfile" ]; then
        echo "Skipping $filename — output already exists."
        continue
    fi

    echo "Processing: $filename"

    num_groups=$(grep -c '^\[' "$ndxfile")
    arr2=($(seq 0 $((num_groups-1))))

    gmx distance \
        -f "$trajfile" \
        -s "$tprfile" \
        -n "$ndxfile" \
        -select "${arr2[@]}" \
        -b 20000 \
        > "$outfile"
done
