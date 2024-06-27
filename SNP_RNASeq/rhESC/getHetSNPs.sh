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
module load bcftools/1.9

REF="/shared/projects/primate_hic/genomes/homo_sapiens/hg38/hg38.fa"
MASKEDREF="/shared/projects/primate_hic/allelic_hic/hg38.na12878.masked.fa"
VCF="/shared/projects/primate_hic/allelic_hic/snp_na12878/NA12878.vcf.gz"
OUTVCF="/shared/projects/primate_hic/allelic_hic/snp_na12878/NA12878.SNPs.Het.PASS.vcf"
BOWTIE2="/shared/projects/primate_hic/allelic_hic/bowtie2/hg38_na12878"

mkdir bowtie2

#Get Het SNPs
bcftools view -v snps -g het -f PASS -O v -o ${OUTVCF} ${VCF}

#Mask the reference fasta
bedtools maskfasta -fi ${REF} -bed ${OUTVCF} -fo ${MASKEDREF}

#Index the masked fasta
samtools faidx ${MASKEDREF}

#Create index folder
#mkdir $2

# index any genome with bowtie2
bowtie2-build --threads 20 ${MASKEDREF} ${BOWTIE2}


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`

