#!/bin/bash
#$ -N EW1_resume
#$ -cwd
#$ -pe sharedmem 6
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -o /exports/eddie/scratch/snumber/logs/
#$ -e /exports/eddie/scratch/snumber/logs/

. /etc/profile.d/modules.sh
module load igmm/apps/cellranger/7.0.0

cd /exports/eddie/scratch/s2845297

cellranger multi \
  --id=EW1_output \
  --csv=EW1_multi.csv \
  --localcores=6 \
  --localmem=45
