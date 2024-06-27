#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=SliceFastq

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
#SBATCH --output=SplitFastq-%j.out

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

module load seqkit/2.1.0

READS1=$1
READS2=$2
SLICE=$3

mkdir sliced_${SLICE}

seqkit split2 -1 ${READS1} -2 ${READS2} -s ${SLICE} -O sliced_${SLICE} -f
##-f overwrittes the output directoty

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`
runtime=$((end-start0))
minute=60
echo "---- Total runtime $runtime s ; $((runtime/minute)) min ----"