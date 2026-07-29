#!/bin/bash
#$ -N cr_multi
#$ -cwd

# 使用 6 核
#$ -pe sharedmem 6

# 最长运行时间
#$ -l h_rt=48:00:00

# 每核 8G，总共 48G
#$ -l h_vmem=8G

# Job array：三个样本并行
#$ -t 1-3

# 最多同时跑 3 个
#$ -tc 3

# 日志
#$ -o logs/
#$ -e logs/

# 邮件通知
#$ -M s2845297@ed.ac.uk
#$ -m be

####################################
# 初始化环境
####################################

. /etc/profile.d/modules.sh

module load igmm/apps/cellranger/7.0.0

####################################
# 路径设置
####################################

WORKDIR=/exports/eddie/scratch/s2845297

cd $WORKDIR

####################################
# 样本列表
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
# 完成
####################################

echo "====================================="
echo "${SAMPLE} completed"
echo "End time: $(date)"
echo "====================================="
