#!/bin/bash

# ATAC-seq preprocessing and mapping
# Dataset:PRJNA725028

GSE=PRJNA725028
SRR_LIST=$(cat reads/${GSE}/SRR_list_${GSE})
READS_DIRECTORY=reads/${GSE}
SRA_BACKUP_DIRECTORY=backup/sra/${GSE}
THREADS=15

# 1. Download .sra from NCBI

for i in $SRR_LIST; do
	prefetch -p $i -O READS_DIRECTORY;
	fasterq-dump -v --threads 15 READS_DIRECTORY/$i/${i}.sra -O READS_DIRECTORY/$i/; #extract .fastq files from .sra and save them in the same folder
	gzip -v READS_DIRECTORY/$i/${i}.fastq;
	mkdir -p SRA_BACKUP_DIRECTORY/$i;
	mv -v READS_DIRECTORY/$i/${i}.sra SRA_BACKUP_DIRECTORY/$i;
        #do echo $i;
done

# 2. Mapping with Bowtie2

for i in $(ls READS_DIRECTORY/SRR*/*_*.fastq.gz); do
	gzip -dv $i; 
done

for i in $(ls READS_DIRECTORY/SRR*/*fastq | perl -ne 's/_\d\.fastq//; print' | sort | uniq); do
	F=$i"_1.fastq"
	R=$i"_2.fastq"
	index="ref/Bwt2_indexed_Human_chm13_v2"
	echo "Mappping $i"

## Mapping
	bowtie2-align-s -p 55 -x $index"/chm13v2.0" -1 $F -2 $R -S mapped/$i.sam 

## Converting SAM to BAM
	if [ -f mapped/$i".sam" ]
	then
	samtools view -bS -@ 55 "mapped/"$i".sam" > "mapped/"$i".bam"
	else
        echo -e "\nSomething went wrong during the mapping of $i" || exit
	fi
done

# 3. Processing BAM file

for i in $(ls mapped | grep .sam | perl -ne 's/\.sam//; print'); do
	echo -e "\n\nSorting $i\n\n"
	samtools sort -@ 55 mapped/$i".bam" -o "mapped/"$i"_sorted.bam" 
# | samtools index -b --threads 15 - "mapped/"$i".bai" 
    # Getting indices for the following step as in
    # https://github.com/CebolaLab/ATAC-seq
	samtools index -b -@ 55 mapped/$i".bam" -o "mapped/"$i".bai"
	if [ -f mapped/$i"_sorted.bam" ]
	then
		rm mapped/$i".sam"
		rm mapped/$i".bam"
	else
		echo -e "\nSomething went wrong during the processing of $i" || exit
	fi
done
