#!/bin/bash
#SBATCH --job-name=bender
#SBATCH --output=bender.out
#SBATCH --error=bender.err
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --partition=gpu

conda activate cellbender

BASE_DIR="/home/jwang/file/Dafei/a20/cellbender"

for FOLDER in "$BASE_DIR"/*; do
    if [ -d "$FOLDER" ]; then
        INPUT_FILE="$FOLDER/raw_feature_bc_matrix.h5"
        OUTPUT_FILE="$FOLDER/$(basename "$FOLDER")_cellbender_output.h5"
        COMPRESSED_FILE="$FOLDER/$(basename "$FOLDER")_cellbender_output_compressed.h5"

        if [ -f "$INPUT_FILE" ]; then
            echo "Processing: $(basename "$FOLDER")"

            cellbender remove-background --cuda \
                --input "$INPUT_FILE" \
                --output "$OUTPUT_FILE"


            if [ -f "$OUTPUT_FILE" ]; then

                ptrepack --complevel 5 "$OUTPUT_FILE":/matrix "$COMPRESSED_FILE":/matrix
                echo "Compression complete: $COMPRESSED_FILE"
            else
                echo "Error: $OUTPUT_FILE was not created. Skipping compression."
            fi

            echo "Complete: $OUTPUT_FILE"
        else
            echo "Skipping $(basename "$FOLDER"), missing input file: $INPUT_FILE"
        fi
    fi
done