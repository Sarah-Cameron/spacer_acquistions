#!/bin/bash

seqkit grep -f $1_matches.BLAST.Result.txtAll_Tags_New.csvread_ids.txt $1_trimmed_rename.fasta -o $1_acquisition_reads.fasta

bwa-mem2 index $2

bwa-mem2 mem $2 $1_matches.BLAST.Result.txtAll_Tags_New.csvUnique_Tags_New.fasta > $1.sam

samtools sort $1.sam -n -o $1.bam

samtools sort $1.bam -o sorted_$1.bam

samtools index sorted_$1.bam

samtools depth -a sorted_$1.bam -o $1_read_depth.coverage 

bwa-mem2 mem $2 $1_matches.BLAST.Result.txtAll_Tags_New.csvTotal_Tags_New.fasta > $1_all.sam

samtools sort $1_all.sam -n -o $1_all.bam

samtools sort $1_all.bam -o sorted_$1_all.bam

samtools index sorted_$1_all.bam

samtools depth -a sorted_$1_all.bam -o $1_all_read_depth.coverage
