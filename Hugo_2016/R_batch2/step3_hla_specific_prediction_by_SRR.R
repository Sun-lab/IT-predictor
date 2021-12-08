
#----------------------------------------------------------------
# Step 3. Match predicted HLA per subject/SRR for hugo Data
#----------------------------------------------------------------

#----------------------------------------------------------------
# Step 1. Import data 
#----------------------------------------------------------------

library(magrittr)
library(Biostrings)
library(dplyr)

# patient info with SRR info (should match with hugo_mutdata3_sm)
hugo.patient.mut = read.table("patient_info_with_mutations.txt", header=TRUE, sep="\t")
dim(hugo.patient.mut)
names(hugo.patient.mut)
hugo.patient.mut[1:2,1:26]

hugo_matchID_SRR = hugo.patient.mut[!duplicated(hugo.patient.mut$matchID),c(1:3, 19:20)]
length(unique(hugo_matchID_SRR$SRA.normal.WES)) + length(unique(hugo_matchID_SRR$SRA.tumor.WES))

# load OptiType results
hugo_pred_hla = read.delim(file="hla_pred/hugo_optitype_result.txt", sep= "\t")
dim(hugo_pred_hla)
hugo_pred_hla[1:2,]
names(hugo_pred_hla)[1] = "SRR"
# hugo_pred_hla$SRR_ref = row.names(hugo_pred_hla)

# load HLA-HD results
load("hla_pred/hugo_hlahd_result.RData") 
ls()
hugo_hlahd_out = hlahd_out
length(hugo_hlahd_out)
hugo_hlahd_out[1]

table(hugo_pred_hla$SRR %in% names(hugo_hlahd_out))
x = setdiff(names(hugo_hlahd_out), hugo_pred_hla$SRR)
x

# do not need to fill the optiType HLA information 
# with HLA-HD results since no missing!

#------------------------------------------------------------------------
# Step 2. Manage predictions and create 'hugo_pateintID_with_HLA.txt'
#------------------------------------------------------------------------

#note: these are different bc Vand_Pt27 and Vand_Pt27_2 have same SRR_ref
length(unique(hugo.patient.mut$SRA.normal.WES))
length(unique(hugo.patient.mut$matchID))

#note: this includes ref and mut SRR 
length(unique(hugo_pred_hla$SRR))

#------hla1
names(hugo_pred_hla)[1]= "SRA.normal.WES"
hugo_pt_with_hla = merge(hugo.patient.mut[1:26], hugo_pred_hla[,1:7], by="SRA.normal.WES", 
                         all.x = TRUE, all.y=FALSE)
dim(hugo_pt_with_hla)
hugo_pt_with_hla[1:2,]

# remove * in the HLA 
hugo_pt_with_hla = hugo_pt_with_hla[,c(2, 1, 3:ncol(hugo_pt_with_hla))]
hugo_pt_with_hla$A1 = gsub("[\\*, :]", "", hugo_pt_with_hla$A1)
hugo_pt_with_hla$A2 = gsub("[\\*, :]", "", hugo_pt_with_hla$A2)
hugo_pt_with_hla$B1 = gsub("[\\*, :]", "", hugo_pt_with_hla$B1)
hugo_pt_with_hla$B2 = gsub("[\\*, :]", "", hugo_pt_with_hla$B2)
hugo_pt_with_hla$C1 = gsub("[\\*, :]", "", hugo_pt_with_hla$C1)
hugo_pt_with_hla$C2 = gsub("[\\*, :]", "", hugo_pt_with_hla$C2)

hugo_pt_with_hla[1:2,]

#------hla2
hugo_hlahd_out[[1]]
hugo_pred_hlaii= as.data.frame(matrix(NA, nrow=nrow(hugo.patient.mut), ncol=11))

colnames(hugo_pred_hlaii)= c("SRA.normal.WES", "DRB11", "DRB12", "DQA11", "DQB11", 
                        "DQA12", "DQB12", "DPA11", "DPB11", "DPA12", "DPB12")

for(i in 1:nrow(hugo.patient.mut)){
  srr = toString(hugo.patient.mut[i,"SRA.normal.WES"])
  hugo_pred_hlaii$SRA.normal.WES[i] = srr
  hugo_pred_hlaii$DRB11[i] = hugo_hlahd_out[[srr]]["DRB1",1]
  hugo_pred_hlaii$DRB12[i] = hugo_hlahd_out[[srr]]["DRB1",2]
  hugo_pred_hlaii$DQA11[i] = hugo_hlahd_out[[srr]]["DQA1",1]
  hugo_pred_hlaii$DQA12[i] = hugo_hlahd_out[[srr]]["DQA1",2]
  hugo_pred_hlaii$DQB11[i] = hugo_hlahd_out[[srr]]["DQB1",1]
  hugo_pred_hlaii$DQB12[i] = hugo_hlahd_out[[srr]]["DQB1",2]
  hugo_pred_hlaii$DPA11[i] = hugo_hlahd_out[[srr]]["DPA1",1]
  hugo_pred_hlaii$DPA12[i] = hugo_hlahd_out[[srr]]["DPA1",2]
  hugo_pred_hlaii$DPB11[i] = hugo_hlahd_out[[srr]]["DPB1",1]
  hugo_pred_hlaii$DPB12[i] = hugo_hlahd_out[[srr]]["DPB1",2]
}

#remove repeat observations 
dim(hugo_pred_hlaii)
hugo_pred_hlaii = unique(hugo_pred_hlaii)
dim(hugo_pred_hlaii)

# check that it is right
hugo_pred_hlaii[25,]
hugo_hlahd_out[["SRR4289724"]]

#----------merge hlai and hlaii by SRR
hugo_pt_with_hla = merge(hugo_pt_with_hla, hugo_pred_hlaii, by="SRA.normal.WES")
hugo_pt_with_hla = hugo_pt_with_hla[, c(2,1, 3:ncol(hugo_pt_with_hla))]
dim(hugo_pt_with_hla)
hugo_pt_with_hla[1:2,]

# format hlaii for netmhc-ii prediction
hugo_pt_with_hla$DRB11 = gsub("[:]", "", substr(hugo_pt_with_hla$DRB11, 5, 14))
hugo_pt_with_hla$DRB11 = gsub("[\\*]", "_", hugo_pt_with_hla$DRB11)
hugo_pt_with_hla$DRB12 = gsub("[:]", "", substr(hugo_pt_with_hla$DRB12, 5, 14))
hugo_pt_with_hla$DRB12 = gsub("[\\*]", "_", hugo_pt_with_hla$DRB12)

hugo_pt_with_hla$DQA11 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DQA11, 5, 14))
hugo_pt_with_hla$DQA12 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DQA12, 5, 14))

hugo_pt_with_hla$DQB11 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DQB11, 5, 14))
hugo_pt_with_hla$DQB12 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DQB12, 5, 14))

hugo_pt_with_hla$DPA11 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DPA11, 5, 14))
hugo_pt_with_hla$DPA12 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DPA12, 5, 14))

hugo_pt_with_hla$DPB11 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DPB11, 5, 14))
hugo_pt_with_hla$DPB12 = gsub("[:, \\*]", "", substr(hugo_pt_with_hla$DPB12, 5, 14))

hugo_pt_with_hla[1:2,]
#output table 
write.table(hugo_pt_with_hla, file = "hugo_patientID_with_HLA.txt", 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

sessionInfo()
q(save="no")

