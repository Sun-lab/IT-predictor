
library(data.table)
library(stringr)

# ------------------------------------------------------------------------
# first check the number of cells per sample for Sade_Feldman_2018
# ------------------------------------------------------------------------

dir = "../scRNAseq/Sade_Feldman_2018/"
fnm = paste0(dir, "GSE120575_patient_ID_single_cells_rm_headers.txt")

info.SF = fread(fnm)
dim(info.SF)
info.SF[1:2,]

names(info.SF) = gsub("characteristics: ", "", names(info.SF))
names(info.SF)[5] = "sampleID"

table(info.SF$`molecule`)
table(info.SF$`description`)
table(info.SF$`processed data file`)
table(info.SF$`raw file`)

info.SF = info.SF[,1:7]

pIDs = str_split_fixed(info.SF$sampleID, "_", 2)
dim(pIDs)
pIDs[1:2,]

info.SF$patientID = pIDs[,2]
info.SF$type = pIDs[,1]

dim(info.SF)
info.SF[1:2,]

length(unique(info.SF$title))
length(unique(info.SF$patientID))

info.SF.ind = unique(info.SF[,.(patientID, sampleID, type, response, therapy)])
dim(info.SF.ind)
info.SF.ind[1:2,]

table(info.SF.ind$response, info.SF.ind$therapy)
table(info.SF.ind$response, info.SF.ind$type)

sort(table(info.SF$sampleID))
table(tapply(info.SF.ind$type, info.SF.ind$patientID, paste, collapse=":"))

# ------------------------------------------------------------------------
# Next check the data Jerby_Arnon_2018
# ------------------------------------------------------------------------

dir = "../scRNAseq/Jerby_Arnon_2018/"
fnm = paste0(dir, "GSE115978_cell.annotations.csv.gz")

info.JA = fread(fnm)
dim(info.JA)
info.JA[1:2,]

length(unique(info.JA$cells))
length(unique(info.JA$samples))

# there are 31 individuals, but 33 samples

info.JA.ind = unique(info.JA[,.(samples,Cohort,treatment.group)])
dim(info.JA.ind)
info.JA.ind = info.JA.ind[order(info.JA.ind$samples),]
info.JA.ind


sort(table(info.JA$samples))
table(info.JA.ind$Cohort)
table(info.JA.ind$treatment.group)
table(info.JA.ind[,.(Cohort,treatment.group)])

sessionInfo()
q(save="no")


