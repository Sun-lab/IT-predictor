#!/bin/bash

# first run the commands in file handle_HLA_GATK_hg38.sh to prepare fasta file 
# without HLA and the interval of the genomic regions without HLA

sampleName="SRR4289744"
picard="/fh/fast/sun_w/bin/picard.jar"
samtools="/fh/fast/sun_w/bin/samtools/bin/samtools"
projDir="/fh/fast/sun_w/research/Immuno/data/Hugo_2016/bams_processed/"
javaTmpDir="/fh/fast/sun_w/tmp4java/"
reference="/fh/fast/sun_w/research/data/human/hg38/Homo_sapiens_assembly38.fasta"
interval="/fh/fast/sun_w/research/data/human/hg38/Homo_sapiens_assembly38_without_HLA.interval_list"

ml java

# filter out reads that are not located on HLA regions
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${picard} FilterSamReads \
 I=${projDir}${sampleName}_sorted_dedup_realigned_recaled.bam \
 O=${projDir}${sampleName}_without_HLA.bam \
 INTERVAL_LIST=${interval} \
 FILTER=includePairedIntervals \
 WRITE_READS_FILES=false &&
#
#
${samtools} view -H \
 ${projDir}${sampleName}_without_HLA.bam \
 | grep -v "SN:HLA" > ${projDir}${sampleName}_header.sam &&
#
#
java -Xmx7g -Djava.io.tmpdir=${javaTmpDir} \
 -jar ${picard} ReplaceSamHeader \
 I=${projDir}${sampleName}_without_HLA.bam \
 HEADER=${projDir}${sampleName}_header.sam \
 O=${projDir}${sampleName}_without_HLA_fixed_header.bam &&
#
#
${samtools} index ${projDir}${sampleName}_without_HLA_fixed_header.bam &&
#
#
rm ${projDir}${sampleName}_header.sam &&
#
#
rm ${projDir}${sampleName}_without_HLA.bam

