#!/bin/bash
#$ -N cr_multi
#$ -cwd

# us 6 cores
#$ -pe sharedmem 6

# th maximum running time
#$ -l h_rt=48:00:00

# each core using 8G, overall using 48G
#$ -l h_vmem=8G

# Job array:run 3saple in parallel
#$ -t 1-3

# top running: run 3 sample at one time
#$ -tc 3

# logs
#$ -o logs/
#$ -e logs/

# email notification
#$ -M s2845297@ed.ac.uk
#$ -m be

####################################
# environent initialization
####################################

. /etc/profile.d/modules.sh

module load igmm/apps/cellranger/7.0.0

####################################
# working directory setup
####################################

WORKDIR=/exports/eddie/scratch/s2845297

cd $WORKDIR

####################################
# sample list
####################################

SAMPLES=("EW1" "EW2" "EW3")

SAMPLE=${SAMPLES[$SGE_TASK_ID-1]}

echo "====================================="
echo "Starting Cell Ranger Multi"
echo "Sample: ${SAMPLE}"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "====================================="

####################################
# Cell Ranger Multi
####################################

cellranger multi \
    --id=${SAMPLE}_output \
    --csv=${WORKDIR}/${SAMPLE}_multi.csv \
    --localcores=6 \
    --localmem=42

####################################
# finish notification
####################################

echo "====================================="
echo "${SAMPLE} completed"
echo "End time: $(date)"
echo "====================================="
