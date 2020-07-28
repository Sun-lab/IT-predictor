
library(data.table)

# ------------------------------------------------------------------------
# first check the number of cells per sample for Sade_Feldman_2018
# ------------------------------------------------------------------------

dir = "../scRNAseq/Sade_Feldman_2018/"
fnm = paste0(dir, "GSE120575_patient_ID_single_cells_rm_headers.txt")

info.SF = fread(fnm)
dim(info.SF)
info.SF[1:2,]

length(unique(info.SF$title))

names(info.SF) = gsub("characteristics: ", "", names(dat))
names(info.SF)[5] = "patientID"

sort(table(info.SF$patientID))

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

sort(table(info.JA$samples))

table(info.JA$Cohort)
info.JA.ind = unique(info.JA[,.(samples,Cohort,treatment.group)])
dim(info.JA.ind)
table(info.JA.ind$Cohort)
table(info.JA.ind$treatment.group)
table(info.JA.ind[,.(Cohort,treatment.group)])

sessionInfo()
q(save="no")


