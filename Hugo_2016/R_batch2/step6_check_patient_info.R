
#--------------------------------------------------------------------
# Step 6: check patient information
#
# _supp/mmc2.xlsx is the supplementary table 2 of the paper
#
# there are also two clinial or sample information from 
# https://github.com/hugon/bms038_analysis/blob/master/data/: 
# bms038_clinical_data.csv and SampleTableCorrected.9.19.16.csv 
#--------------------------------------------------------------------

library(magrittr)
library(stringi)
library(Biostrings)
library(openxlsx)

#--------------------------------------------------------------------
# 0. patient with neoantigen information
#--------------------------------------------------------------------

sample_mb = read.table(file = "hugo_patient_mb_info.txt", 
                               header = TRUE, sep = " ")
dim(sample_mb)
sample_mb[1:2,]

sample_with_neoAg = read.table(file = "hugo_nMut_with_peptides.txt", 
                              header = TRUE, sep = " ")
dim(sample_with_neoAg)
sample_with_neoAg[1:2,]

names(sample_with_neoAg)[1] = "sample"
sample_mb = merge(sample_mb, sample_with_neoAg, by="sample", all.x=TRUE)
dim(sample_mb)
sample_mb[1:2,]

summary(sample_mb$nMuts_with_peptides)
sample_mb$nMuts_with_peptides[which(is.na(sample_mb$nMuts_with_peptides))] = 0

summary(sample_mb$n_mutations_neoAg - sample_mb$nMuts_with_peptides)

sample_mb$lowerID = tolower(sample_mb$sample)
dim(sample_mb)
sample_mb[1:2,]

#--------------------------------------------------------------------
# 2. upload clinical data 
#--------------------------------------------------------------------
hugo.patient.info = read.delim("patient_info_with_mutations.txt", sep="\t")
hugo.patient.info[1:2, 1:26]
dim(hugo.patient.info)

#--------------------------------------------------------------------
# 3. check response information from sample file
#--------------------------------------------------------------------
table(hugo.patient.info$irRECIST)

#remove Vand_Pt27_2
hugo.patient.info=hugo.patient.info[which(hugo.patient.info$matchID !="Vand_Pt27_2"),]
dim(hugo.patient.info)
table(hugo.patient.info$irRECIST)

sessionInfo()
q(save="no")

