#!/bin/bash

sampleID="Vand_Pt38"
nSample="SRR4289743"
tSample="SRR4289744"

muTect="/fh/fast/sun_w/bin/muTect-1/muTect-1.1.7.jar"
bamDir="/fh/fast/sun_w/research/Immuno/data/Hugo_2016/bams_processed/"
outDir="/fh/fast/sun_w/research/Immuno/data/Hugo_2016/muTect/"
reference="/fh/fast/sun_w/research/data/human/hg38/Homo_sapiens_assembly38.fasta"
gatkBundle="/fh/fast/sun_w/research/data/GATK_bundle/hg38/"
cosmicDir="/fh/fast/sun_w/research/data/human/COSMIC/GRCH38_v83/"

ml java/jdk1.7.0_51

java -Xmx7g -Djava.io.tmpdir=/fh/fast/sun_w/tmp4java/ \
 -jar /fh/fast/sun_w/bin/muTect-1/muTect-1.1.7.jar \
 --disable_auto_index_creation_and_locking_when_reading_rods \
 --analysis_type MuTect \
 --reference_sequence ${reference} \
 --cosmic ${cosmicDir}Cosmic_sorted_by_GATK_hg38.vcf \
 --dbsnp ${gatkBundle}dbsnp_146.hg38.vcf.gz \
 --input_file:normal ${bamDir}${nSample}_sorted_dedup_realigned_recaled.bam \
 --input_file:tumor ${bamDir}${tSample}_sorted_dedup_realigned_recaled.bam \
 --out ${outDir}${sampleID}_hg38_muTect_call_stats.txt\
 --vcf ${outDir}${sampleID}_hg38_muTect.vcf
 

