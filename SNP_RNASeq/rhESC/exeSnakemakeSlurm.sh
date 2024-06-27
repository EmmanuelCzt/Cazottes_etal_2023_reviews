#!/bin/bash

module load samtools/1.15.1
module load snakemake/5.7.4
module load picard/2.22.0
module load bcftools/1.16
module load gatk4/4.2.0.0
module load bedtools/2.30.0
# module load sra-tools/2.10.3
# module load star/2.7.2b
# module load kallisto/0.46.2
# module load htseq/0.13.5

snakemake --unlock
snakemake -rp -j 200 --resources load=100 --cluster "sbatch -p {cluster.partition} --cpus-per-task {cluster.cpu} --mem {cluster.ram}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.yaml
