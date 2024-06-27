#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=Bowtie2indexForAllelicHiCMapping

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
#SBATCH --output=Bowtie2index-%j.out

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

# modules loading
module load bowtie2/2.4.4
module load bedtools/2.30.0
module load samtools/1.13


REF="/shared/projects/primate_hic/genomes/homo_sapiens/hg38/hg38.fa"
MASKEDREF="/shared/projects/primate_hic/captureHiC/h9/genome/hg38.masked.fa"
VCF="/shared/projects/primate_hic/captureHiC/h9/snp_h9/h9.wgs.snp.biallelic.PASS.Het.chrX.Phased.Xa.Xi.vcf"
BOWTIE2="/shared/projects/primate_hic/captureHiC/h9/bowtie2/hg38_masked"

mkdir -p genome bowtie2

#Mask the reference fasta
bedtools maskfasta -fi ${REF} -bed ${VCF} -fo ${MASKEDREF}

#Index the masked fasta
samtools faidx ${MASKEDREF}

#Create index folder
#mkdir $2

# index any genome with bowtie2
bowtie2-build --threads 20 ${MASKEDREF} ${BOWTIE2}


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`
