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
ls -lh EW2_metadata_expression_only.csv
ls -lh EW2_pca_expression_only.csv
ls -lh EW2_umap_expression_only.csv

echo "run PAGA:"
$PYTHON paga_final_script.py \
  --metadata EW2_metadata_expression_only.csv \
  --pca EW2_pca_expression_only.csv \
  --umap EW2_umap_expression_only.csv \
  --outdir /exports/eddie/scratch/s2845297/PAGA/paga_expression_only_dpt \
  --output-prefix EW2_paga_expression_only \
  --cell-id-col cell_id \
  --celltype-col celltype \
  --group-col celltype \
  --n-neighbors 15 \
  --n-pcs 100 \
  --paga-threshold 0.03 \
  --root-group "NAIVE 2" \
  --dpt-root-strategy centroid \
  --random-state 0
echo "finish："
date

