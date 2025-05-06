#!/bin/sh
#SBATCH --job-name=cyto
#SBATCH --output=cyto.out
#SBATCH --error=cyto.err
#SBATCH --time=24:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=12

###library(DropletUtils)
###imm<-readRDS("./")
###write10xCounts(x = imm@assays$RNA@counts, path = "./raw", version="3")
###write.table(imm@meta.data, file.path("./raw", 'metadata'), sep='\t', quote=F)


/mount/gfs/jwang/A20/cytosig/CytoSig/CytoSig/CytoSig_run.py \
-o ./output  -s 1  -c 0 \
-i /mount/gfs/jwang/A20/cytosig/data/all/barcodes.tsv.gz, \
/mount/gfs/jwang/A20/cytosig/data/all/features.tsv.gz, \
/mount/gfs/jwang/A20/cytosig/data/all/matrix.mtx.gz
