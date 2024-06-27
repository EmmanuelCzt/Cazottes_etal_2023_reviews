#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=TADcalling

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
#SBATCH --output=TADcalling-%j.out

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

module load hicexplorer/3.7.2

MCOOL=$1 # Path to the multi-resolution cooler file (.mcool)
OUTDIR=$2

BED2WIG=("/shared/projects/primate_hic/captureHiC/rhes/WT/scripts/bedGraphToBigWig")
CHROM=("/shared/projects/primate_hic/genomes/rhemac10/rheMac10chrom.noScaff.sorted.sizes")


# Check arguments
if [ -z ${MCOOL} ]; then
  echo "Error: Please provide the path to the multi-resolution cooler file (.mcool) as the first argument."
  exit 1
fi

if [ -z ${OUTDIR} ]; then
  OUTDIR="."
fi

mkdir -p ${OUTDIR}

#Extract cool absolute paths
COOL=$(cooler ls ${MCOOL})


# Loop through each resolution
for i in ${COOL}; 
do
	# Extract resolutions
	RESOLUTION=$(basename $i)

  mkdir ${OUTDIR}/${RESOLUTION}

  STEP=$((${RESOLUTION}*2))
  MIN=$((${RESOLUTION}*3))
  MAX=$((((${RESOLUTION}*10)+${STEP})))

  PREFIX=$(echo "rhES_WT"\_min${MIN}\_max${MAX}\_step${STEP}\_bin${RESOLUTION})

  hicFindTADs -m $i \
  --outPrefix ${OUTDIR}/${RESOLUTION}/${PREFIX} \
  --minDepth ${MIN} --maxDepth ${MAX} --step ${STEP} --correctForMultipleTesting fdr

  ${BED2WIG} ${OUTDIR}/${RESOLUTION}/${PREFIX}\_score.bedgraph ${CHROM} \
  ${OUTDIR}/${RESOLUTION}/${PREFIX}\_score.bw
  #Generate contact map using cooler show
done


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`

