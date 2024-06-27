#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=runHiCpro

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=168:00:00

### Requirements
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
#SBATCH --output=hicpro-%j.out

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

start0=`date +%s`

#Path to HiC pro singularity img
IMG="/shared/software/singularity/images/hicpro-3.1.0.sif"
#Path to raw files
RAW="/shared/projects/primate_hic/captureHiC/rhes/delHERVK/raw"
#Folder path for outputs
OUTPATH="/shared/projects/primate_hic/captureHiC/rhes/delHERVK/hicpro_results"
#Path to config file
CONFIG="/shared/projects/primate_hic/captureHiC/rhes/delHERVK/config-hicpro.txt"


singularity exec ${IMG} HiC-Pro -i ${RAW} -o ${OUTPATH} -c ${CONFIG}


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`
runtime=$((end-start0))
minute=60

