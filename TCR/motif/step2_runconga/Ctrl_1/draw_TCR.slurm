#!/bin/sh
#SBATCH --job-name=draw_p20
#SBATCH --output=draw_p20.out
#SBATCH --error=draw_p20.err
#SBATCH --time=24:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=12



script_path="/mount/gfs/jwang/A20/conga/conga/scripts/make_tcr_logos.py"

output_prefix="/mount/gfs/jwang/A20/conga/ctrl/ctrl1/"

python /mount/gfs/jwang/A20/conga/conga/scripts/setup_10x_for_conga.py --filtered_contig_annotations_csvfile ./filtered_contig_annotations.csv --organism mouse

