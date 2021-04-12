
#--------------------------------------------------------------------
# Step 6: check patient information
#
# _supp/mmc2.xlsx is the supplementary table 2 of the paper
#
# there are also two clinial or sample information from 
# https://github.com/riazn/bms038_analysis/blob/master/data/: 
# bms038_clinical_data.csv and SampleTableCorrected.9.19.16.csv 
#--------------------------------------------------------------------

library(magrittr)
library(stringi)
library(Biostrings)
library(openxlsx)

#--------------------------------------------------------------------
# 0. patient with neoantigen information
#--------------------------------------------------------------------

sample_mb = read.table(file = "riaz_patient_mb_info.txt", 
                               header = TRUE, sep = " ")
dim(sample_mb)
sample_mb[1:2,]

sample_with_neoAg = read.table(file = "riaz_nMut_with_peptides.txt", 
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
sample_mb$PreOn = sub(".*_", "", sample_mb$sample)
sample_mb$Patient = sub("_.*", "", sample_mb$sample)
dim(sample_mb)
sample_mb[1:2,]

table(sample_mb$PreOn)

#--------------------------------------------------------------------
# 1. compare clinical information from paper supp and GitHub
#--------------------------------------------------------------------
# import patient info 

supp.clinic = read.xlsx("_supp/mmc2.xlsx", startRow=2)
dim(supp.clinic)
supp.clinic[1:2,]

riaz.clinic = read.csv("_github/bms038_clinical_data.csv")
dim(riaz.clinic)
riaz.clinic[1:2,]
table(riaz.clinic$SampleType)

table(supp.clinic$Patient == riaz.clinic$PatientID, useNA="ifany")
table(supp.clinic$Response, riaz.clinic$BOR,        useNA="ifany")
table(supp.clinic$Response, riaz.clinic$myBOR,      useNA="ifany")
table(supp.clinic$Response, riaz.clinic$myBOR2,     useNA="ifany")
table(supp.clinic$Response, riaz.clinic$myBOR3,     useNA="ifany")
table(supp.clinic$Response, riaz.clinic$IBOR,       useNA="ifany")

# so BOR is the same as response, 
# myBOR* are different ways to collapse patient groups
# it is not clear what is IBOR

table(sample_mb$Patient %in% supp.clinic$Patient)

#--------------------------------------------------------------------
# 2. compare clinical and sample information from GitHub
#--------------------------------------------------------------------

riaz.sample = read.csv("_github/SampleTableCorrected.9.19.16.csv")
dim(riaz.sample)
riaz.sample[1:2,]
table(riaz.sample$X == riaz.sample$Sample)
riaz.sample = riaz.sample[,-1]
dim(riaz.sample)
riaz.sample[1:2,]

table(table(riaz.sample$PatientID))
table(table(riaz.clinic$PatientID))

table(unique(tolower(riaz.sample$PatientID)) %in% tolower(riaz.clinic$PatientID))
setdiff(riaz.sample$PatientID, riaz.clinic$PatientID)

table(tolower(riaz.clinic$PatientID) %in% tolower(riaz.sample$PatientID))

riaz.clinic$Sample = tolower(riaz.clinic$Sample)
riaz.sample$Sample = tolower(riaz.sample$Sample)

table(riaz.clinic$Sample %in% riaz.sample$Sample)
setdiff(riaz.clinic$Sample, riaz.sample$Sample)

table(riaz.sample$Sample %in% riaz.clinic$Sample)
setdiff(riaz.sample$Sample, riaz.clinic$Sample)

#--------------------------------------------------------------------
# 3. check response information from sample file
#--------------------------------------------------------------------

table(riaz.sample$BOR, riaz.sample$Response)
response.sample = unique(riaz.sample[,c("PatientID", "BOR")])
dim(response.sample)
response.sample[1:2,]
table(table(response.sample$PatientID))

response2compare = merge(riaz.clinic, response.sample, by="PatientID")
dim(response2compare)
response2compare[1:2,]

table(response2compare$BOR.x, response2compare$BOR.y, useNA="ifany")

sessionInfo()
q(save="no")

