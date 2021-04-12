
#----------------------------------------------------------------
# Step 3. Match predicted HLA per subject/SRR for Riaz Data
#----------------------------------------------------------------

#----------------------------------------------------------------
# Step 1. Import data 
#----------------------------------------------------------------

library(data.table)
library(magrittr)
library(Biostrings)
library(dplyr)

# patient info with SRR info (should match with mutdata3_sm)
riaz.patient.mut = fread("../data/srr_sample2.txt")
dim(riaz.patient.mut)
riaz.patient.mut[1:2,]

# load OptiType results
riaz_pred_hla = fread("../data/riaz_optitype_result_withSRR.txt")
dim(riaz_pred_hla)
riaz_pred_hla[1:2,]
names(riaz_pred_hla)[1] = "SRR_ref"

# load HLA-HD results
load("../data/riaz_hlahd_result.RData") 
ls()
length(riaz_hlahd_out)
riaz_hlahd_out[1]

table(riaz_pred_hla$SRR_ref %in% names(riaz_hlahd_out))
x = setdiff(names(riaz_hlahd_out), riaz_pred_hla$SRR_ref)
x

# fill the optiType HLA information of this sample by the HLA-HD results

r1 = riaz_hlahd_out[[x]]
r1
HLA1 = as.character(t(r1[1:3,]))
HLA1

HLA1 = gsub("HLA-", "", HLA1, fixed=TRUE)
HLA1 = substr(HLA1, 1, 7)

HLA1 = as.data.frame(matrix(HLA1, nrow=1))
HLA1

HLA1 = cbind(x, HLA1, NA, NA)
names(HLA1) = names(riaz_pred_hla)
HLA1

riaz_pred_hla = rbind(riaz_pred_hla, HLA1)
dim(riaz_pred_hla)
riaz_pred_hla[c(1,nrow(riaz_pred_hla)),]

setequal(riaz_pred_hla$SRR_ref, names(riaz_hlahd_out))

#------------------------------------------------------------------------
# Step 2. Manage predictions and create 'riaz_pateintID_with_HLA.txt'
#------------------------------------------------------------------------

length(unique(riaz.patient.mut$SRR_ref))
length(unique(riaz_pred_hla$SRR_ref))

# hla1
riaz_pt_with_hla = merge(riaz.patient.mut, riaz_pred_hla[,1:7], by="SRR_ref")
dim(riaz_pt_with_hla)
riaz_pt_with_hla[1:2,]

riaz_pt_with_hla = as.data.frame(riaz_pt_with_hla)
# remove * in the HLA 
riaz_pt_with_hla = riaz_pt_with_hla[,c(2, 1, 3:ncol(riaz_pt_with_hla))]
riaz_pt_with_hla$A1 = gsub("[\\*, :]", "", riaz_pt_with_hla$A1)
riaz_pt_with_hla$A2 = gsub("[\\*, :]", "", riaz_pt_with_hla$A2)
riaz_pt_with_hla$B1 = gsub("[\\*, :]", "", riaz_pt_with_hla$B1)
riaz_pt_with_hla$B2 = gsub("[\\*, :]", "", riaz_pt_with_hla$B2)
riaz_pt_with_hla$C1 = gsub("[\\*, :]", "", riaz_pt_with_hla$C1)
riaz_pt_with_hla$C2 = gsub("[\\*, :]", "", riaz_pt_with_hla$C2)

# hla2
head(riaz_hlahd_out[[1]])
riaz_pred_hlaii= as.data.frame(matrix(NA, nrow=nrow(riaz.patient.mut), ncol=11))

colnames(riaz_pred_hlaii)= c("SRR_ref", "DRB11", "DRB12", "DQA11", "DQB11", 
                             "DQA12", "DQB12", "DPA11", "DPB11", "DPA12", "DPB12")

for(i in 1:nrow(riaz.patient.mut)){
  srr = toString(riaz.patient.mut[i,"SRR_ref"])
  riaz_pred_hlaii$SRR_ref[i] = srr
  riaz_pred_hlaii$DRB11[i] = riaz_hlahd_out[[srr]]["DRB1",1]
  riaz_pred_hlaii$DRB12[i] = riaz_hlahd_out[[srr]]["DRB1",2]
  riaz_pred_hlaii$DQA11[i] = riaz_hlahd_out[[srr]]["DQA1",1]
  riaz_pred_hlaii$DQA12[i] = riaz_hlahd_out[[srr]]["DQA1",2]
  riaz_pred_hlaii$DQB11[i] = riaz_hlahd_out[[srr]]["DQB1",1]
  riaz_pred_hlaii$DQB12[i] = riaz_hlahd_out[[srr]]["DQB1",2]
  riaz_pred_hlaii$DPA11[i] = riaz_hlahd_out[[srr]]["DPA1",1]
  riaz_pred_hlaii$DPA12[i] = riaz_hlahd_out[[srr]]["DPA1",2]
  riaz_pred_hlaii$DPB11[i] = riaz_hlahd_out[[srr]]["DPB1",1]
  riaz_pred_hlaii$DPB12[i] = riaz_hlahd_out[[srr]]["DPB1",2]
}

# remove repeat observations 
dim(riaz_pred_hlaii)
riaz_pred_hlaii = distinct(riaz_pred_hlaii)
dim(riaz_pred_hlaii)

# check that it is right
riaz_pred_hlaii[25,]
riaz_hlahd_out[["SRR5134763"]]

# -merge hlai and hlaii by SRR
riaz_pt_with_hla = merge(riaz_pt_with_hla, riaz_pred_hlaii, by="SRR_ref")
riaz_pt_with_hla = riaz_pt_with_hla[, c(2,1, 3:ncol(riaz_pt_with_hla))]
dim(riaz_pt_with_hla)
riaz_pt_with_hla[1:2,]

# format hlaii for netmhc-ii prediction
riaz_pt_with_hla$DRB11 = gsub("[:]", "", substr(riaz_pt_with_hla$DRB11, 5, 14))
riaz_pt_with_hla$DRB11 = gsub("[\\*]", "_", riaz_pt_with_hla$DRB11)
riaz_pt_with_hla$DRB12 = gsub("[:]", "", substr(riaz_pt_with_hla$DRB12, 5, 14))
riaz_pt_with_hla$DRB12 = gsub("[\\*]", "_", riaz_pt_with_hla$DRB12)

riaz_pt_with_hla$DQA11 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DQA11, 5, 14))
riaz_pt_with_hla$DQA12 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DQA12, 5, 14))

riaz_pt_with_hla$DQB11 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DQB11, 5, 14))
riaz_pt_with_hla$DQB12 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DQB12, 5, 14))

riaz_pt_with_hla$DPA11 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DPA11, 5, 14))
riaz_pt_with_hla$DPA12 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DPA12, 5, 14))

riaz_pt_with_hla$DPB11 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DPB11, 5, 14))
riaz_pt_with_hla$DPB12 = gsub("[:, \\*]", "", substr(riaz_pt_with_hla$DPB12, 5, 14))

#output table 
write.table(riaz_pt_with_hla, file = "riaz_patientID_with_HLA.txt", 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

sessionInfo()
q(save="no")

