#!/bin/bash

sampleName="SRR3083873"
picard="/fh/fast/sun_w/bin/picard.jar"
gatk="/fh/fast/sun_w/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
gatkBundle="/fh/fast/sun_w/research/data/GATK_bundle/hg38"
projDir="/fh/fast/sun_w/research/Immuno/data/Hugo_2016"
javaTmpDir="/fh/fast/sun_w/tmp4java/"
reference="/fh/fast/sun_w/research/data/human/hg38/Homo_sapiens_assembly38.fasta"

ml java

# Add read groups, sort by coordinate, save as a bam file and index the bam file
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${picard} AddOrReplaceReadGroups \
 INPUT=${projDir}/sams/${sampleName}.sam \
 OUTPUT=${projDir}/bams/${sampleName}_sorted_rg.bam \
 SORT_ORDER=coordinate \
 CREATE_INDEX=true \
 ID=${sampleName}.LANE001 SM=${sampleName} LB=${sampleName} \
 PL=ILLUMINA PU=${sampleName} &&
#
# Mark duplicated reads
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${picard} MarkDuplicates \
 I=${projDir}/bams/${sampleName}_sorted_rg.bam \
 O=${projDir}/bams/${sampleName}_sorted_dedup.bam \
 M=${projDir}/bams/${sampleName}_sorted_dedup_metric.txt \
 ASSUME_SORT_ORDER=coordinate \
 CREATE_INDEX=true &&
#
# RealignerTargetCreator
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${gatk} \
 -T RealignerTargetCreator \
 -R ${reference} \
 -I ${projDir}/bams/${sampleName}_sorted_dedup.bam \
 -o ${projDir}/bams/${sampleName}_sorted_dedup.bam.intervals \
 -known ${gatkBundle}/Homo_sapiens_assembly38.known_indels.vcf.gz \
 -known ${gatkBundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz &&
#
# Realign
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${gatk} \
 -T IndelRealigner \
 -R ${reference} \
 -I ${projDir}/bams/${sampleName}_sorted_dedup.bam \
 -targetIntervals ${projDir}/bams/${sampleName}_sorted_dedup.bam.intervals \
 -o ${projDir}/bams/${sampleName}_sorted_dedup_realigned.bam \
 -known ${gatkBundle}/Homo_sapiens_assembly38.known_indels.vcf.gz \
 -known ${gatkBundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
 --filter_bases_not_stored &&
#
# Build index
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${picard} BuildBamIndex \
 INPUT=${projDir}/bams/${sampleName}_sorted_dedup_realigned.bam &&
#
# Base quality recalibration: create table
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${gatk} \
 -T BaseRecalibrator \
 -R ${reference} \
 -knownSites ${gatkBundle}/dbsnp_146.hg38.vcf.gz \
 -knownSites ${gatkBundle}/Homo_sapiens_assembly38.known_indels.vcf.gz \
 -knownSites ${gatkBundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
 -I ${projDir}/bams/${sampleName}_sorted_dedup_realigned.bam \
 -o ${projDir}/bams/${sampleName}_sorted_dedup_realigned_recal_table.txt &&
#
# Base quality recalibration
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${gatk} \
 -T PrintReads \
 -R ${reference} \
 -I ${projDir}/bams/${sampleName}_sorted_dedup_realigned.bam \
 -BQSR ${projDir}/bams/${sampleName}_sorted_dedup_realigned_recal_table.txt \
 -o ${projDir}/bams/${sampleName}_sorted_dedup_realigned_recaled.bam &&
#
# Build index
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${picard} BuildBamIndex \
 INPUT=${projDir}/bams/${sampleName}_sorted_dedup_realigned_recaled.bam
 

