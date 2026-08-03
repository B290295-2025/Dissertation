#!/bin/bash
#$ -N EW2_paga
#$ -cwd
#$ -pe sharedmem 8
#$ -l h_vmem=8G
#$ -l h_rt=04:00:00
#$ -o paga_job.out
#$ -e paga_job.err

set -euo pipefail

echo "start："
date

echo "running node："
hostname

echo "current working directory："
pwd

PYTHON=/exports/eddie/scratch/$USER/conda_envs/scanpy_paga/bin/python

echo "Python:"
$PYTHON -c "import sys; print(sys.executable)"

echo "test scanpy:"
$PYTHON -c "import scanpy as sc; print(sc.__version__)"

echo "check input files:"
ls -lh paga_input/EW2_metadata.csv
ls -lh paga_input/EW2_pca.csv
ls -lh paga_input/EW2_umap.csv

echo "run PAGA:"
$PYTHON paga_final.py \
  --metadata paga_input/EW2_metadata.csv \
  --pca paga_input/EW2_pca.csv \
  --umap paga_input/EW2_umap.csv \
  --celltype-col celltype \
  --output-prefix EW2_paga \
  --outdir EW2_paga_final \
  --n-neighbors 15 \
  --n-pcs 30 \
  --paga-threshold 0.03 \
  --point-size 8

echo "finish："
date

