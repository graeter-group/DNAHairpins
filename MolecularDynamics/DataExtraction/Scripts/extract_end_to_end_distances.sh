#!/bin/bash
set -euo pipefail
module load gromacs/2025.2

DATA_DIRECTORY_TPRS="../RawData/"
DATA_DIRECTORY_PROCESSED="../ProcessedData/"
INDEX_DIRECTORY="../IndexFiles/"
OUTPUT_DIRECTORY="../../3_Data_Analysis/ExtractedData/"

mkdir -p "$OUTPUT_DIRECTORY"

for tprfile in "$DATA_DIRECTORY_TPRS"/*.tpr; do
    [ -e "$tprfile" ] || continue

    filename=$(basename "$tprfile" .tpr)

    echo "Processing: $filename"

    trajfile="${DATA_DIRECTORY_PROCESSED}/${filename}_Proc.xtc"
    ndxfile="${INDEX_DIRECTORY}/${filename}_EtE.ndx"
    outfile="${OUTPUT_DIRECTORY}/${filename}_EndToEndDistances.xvg"

    gmx distance \
        -f "$trajfile" \
        -s "$tprfile" \
        -n "$ndxfile" \
        -oall "$outfile" \
        -select 0 \
        -nopbc
done
