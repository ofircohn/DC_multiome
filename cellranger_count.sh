#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P sankaranlab
#$ -l h_vmem=40g
#$ -pe smp 12
#$ -binding linear:12
#$ -l h_rt=30:00:00
#$ -e cellrangercount.err
#$ -o cellrangercount.out
#$ -R y

# Load required software
source /broad/software/scripts/useuse
reuse -q .bcl2fastq2-v2.20.0
reuse Python-2.7
reuse UGER
export PATH=/broad/sankaranlab/cohn/apps/cellranger-arc-2.0.2:$PATH

cellranger-arc count --id=DCmultiome \
                       --reference=/broad/sankaranlab/cohn/apps/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/broad/sankaranlab/cohn/apps/new/cellranger-arc-count.csv \
                       --localcores=16 \
                       --localmem=54
