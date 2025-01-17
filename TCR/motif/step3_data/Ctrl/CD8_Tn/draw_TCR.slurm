#!/bin/sh
#SBATCH --job-name=draw_p20
#SBATCH --output=draw_p20.out
#SBATCH --error=draw_p20.err
#SBATCH --time=24:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=12


python /home/jwang/file/conga/p20/conga/scripts/make_tcr_logos.py --tcrs_tsvfile Ctrl_CD8-Tn_TCR_conga.tsv --outfile_prefix ./Ctrl_CD8-Tn_TCR_conga_test --organism mouse
