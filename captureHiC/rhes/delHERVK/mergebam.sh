#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=mergeBams

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
#SBATCH --output=mergebam-%j.out

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

module load samtools/1.15.1

#Path to HiC pro singularity img
IMG="/shared/software/singularity/images/hicpro-3.1.0.sif"
#Path to raw files
RAW="/shared/projects/primate_hic/captureHiC/rhes/delHERVK/raw"
#Folder path for outputs
OUTPATH="/shared/projects/primate_hic/captureHiC/rhes/delHERVK/hicpro_results"
#Path to config file
CONFIG="/shared/projects/primate_hic/captureHiC/rhes/delHERVK/config-hicpro.txt"


singularity exec ${IMG} HiC-Pro -i ${RAW} -o ${OUTPATH} -c ${CONFIG}

samtools merge -@ 16 -n -f hicpro_results/bowtie_results/bwt2/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R1_rheMac10_masked.bwt2merged.bam \
hicpro_results/bowtie_results/bwt2_global/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R1_rheMac10_masked.bwt2glob.bam \
hicpro_results/bowtie_results/bwt2_local/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R1_rheMac10_masked.bwt2glob.unmap_bwt2loc.bam

samtools merge -@ 16 -n -f hicpro_results/bowtie_results/bwt2/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R2_rheMac10_masked.bwt2merged.bam \
hicpro_results/bowtie_results/bwt2_global/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R2_rheMac10_masked.bwt2glob.bam \
hicpro_results/bowtie_results/bwt2_local/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R2_rheMac10_masked.bwt2glob.unmap_bwt2loc.bam

samtools sort -@ 16 -m 62M -n -T hicpro_results/tmp/rhES_E9.9.delHERVK.Hom.R2_rheMac10_masked \
-o hicpro_results/bowtie_results/bwt2/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R2_rheMac10_masked.bwt2merged.sorted.bam \
hicpro_results/bowtie_results/bwt2/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R2_rheMac10_masked.bwt2merged.bam
samtools sort -@ 16 -m 62M -n -T hicpro_results/tmp/rhES_E9.9.delHERVK.Hom.R1_rheMac10_masked \
-o hicpro_results/bowtie_results/bwt2/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R1_rheMac10_masked.bwt2merged.sorted.bam \
hicpro_results/bowtie_results/bwt2/rhES_E9.9/rhES_E9.9.delHERVK.Hom.R1_rheMac10_masked.bwt2merged.bam


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`
runtime=$((end-start0))
minute=60

