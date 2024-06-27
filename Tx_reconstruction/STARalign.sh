#!/bin/bash
# If paired-end RNA-SEQ fastq1 list : coma separated file names. Fastq2 list : space separated file names
# conda activate STARalign

SAMPLE=("SRR_Acc_List.txt")
INDEX=("/home/emmanuel/Documents/Genomes/calJac4/calJac4_STAR/")
ASSEMBLY=("calJac4")
COMP=("zcat")
OUTDIR=("map")
#SJ_tab=$5 # two pass mapping for SNP calling

mkdir ${OUTDIR}

for i in `cat ${SAMPLE}`
do
	STAR --runMode alignReads \
	--genomeLoad LoadAndKeep \
	--outSAMstrandField intronMotif \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonical \
	--readFilesCommand ${COMP} \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir ${INDEX} \
	--readFilesIn fastq/$i\_1.fastq.gz fastq/$i\_2.fastq.gz \
	--runThreadN 12 \
	--limitBAMsortRAM 15000000000 \
	--outFileNamePrefix ${OUTDIR}/`basename $i`\_${ASSEMBLY} 
done


STAR --genomeDir ${INDEX} --genomeLoad Remove

#PS RNA-seq data with a lot of unmapped reads so I will keep them to see what's up
#	--outReadsUnmapped Fastx \
#	--sjdbFileChrStartEnd ${SJ_tab}




