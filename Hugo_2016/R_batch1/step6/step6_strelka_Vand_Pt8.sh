#!/bin/bash
#

sampleID="Vand_Pt8"
nSample="SRR4289714"
tSample="SRR4289715"

strelka="/fh/fast/sun_w/bin/strelka/bin/configureStrelkaWorkflow.pl"
bamDir="/fh/fast/sun_w/research/Immuno/data/Hugo_2016/bams_processed/"
outDir="/fh/fast/sun_w/research/Immuno/data/Hugo_2016/strelka/"
reference="/fh/fast/sun_w/research/data/human/hg38/Homo_sapiens_assembly38_no_HLA.fasta"
configFile="/fh/fast/sun_w/bin/strelka/bin/_strelka_config_bwa_default.ini"

rm -rf ${outDir}${sampleID} &&
#
# 
${strelka} \
--normal=${bamDir}${nSample}_without_HLA_fixed_header.bam \
--tumor=${bamDir}${tSample}_without_HLA_fixed_header.bam \
--ref=${reference} \
--config=${configFile} \
--output-dir=${outDir}${sampleID} &&
# 
#
make -C ${outDir}${sampleID} 

