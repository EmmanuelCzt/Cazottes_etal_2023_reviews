#!/bin/bash

PAIRS2COOL=("/home/emmanuel/software/HiCPro/HiC-Pro/bin/utils/allValidPairs2cooler.sh")
CHROM=("/home/emmanuel/Documents/annotations/chrom.sizes/rheMac10chrom.sizes" "/home/emmanuel/Documents/annotations/chrom.sizes/panTro6chrom.sizes" "/home/emmanuel/Documents/annotations/chrom.sizes/hg38.chrom.sizes")
RES=("250" "500" "750")
OUTDIR=("cooler/merge")
INDIR=("hicproQ30_rhES/hic_results/merge" "hicproQ30_chiPSC/hic_results/merge" "hicproQ30_h9/hic_results/merge")
INFILE=("rhES" "chiPSC" "H9") # 

i=0

mkdir ${OUTDIR}


while [ $i -lt `echo "${#INFILE[@]}"` ]
do
	j=0
	mkdir ${OUTDIR}/${INFILE[$i]}
	while [ $j -lt `echo "${#RES[@]}"` ]
	do
    	mkdir ${OUTDIR}/${INFILE[$i]}/${RES[$j]}
    	${PAIRS2COOL} -i ${INDIR[$i]}/${INFILE[$i]}_ontarget.XIConly.q23.allValidPairs \
    	-c ${CHROM[$i]} -r ${RES[$j]} -p 16 -o ${OUTDIR}/${INFILE[$i]}/${RES[$j]}
    	j=$(($j+1))
	done
	i=$(($i+1))
done
