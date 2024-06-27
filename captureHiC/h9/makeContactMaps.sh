#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=showHiC

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
#SBATCH --output=showHicMap-%j.out

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

module load cooltools/0.5.4

MCOOL=$1 # Path to the multi-resolution cooler file (.mcool)
OUTDIR=$2
HMCOORD="chrX:73077515-75326610"

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
  #RESOLUTION=$(grep -Eo '::/resolutions/[0-9]+' $i | cut -d'/' -f3)
  # Generate output filename
  output_file="${OUTDIR}/$(basename ${MCOOL} .mcool).${RESOLUTION}.pdf"

  #Generate contact map using cooler show

  cooler show -b --dpi 320 -o ${output_file} $i ${HMCOORD}

  echo "Generated contact map for resolution ${RESOLUTION}: ${output_file}"
done

echo "All contact maps generated!"



echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
end=`date +%s`
