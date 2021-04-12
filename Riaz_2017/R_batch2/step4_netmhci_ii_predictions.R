
#------------------------------------------------------
# Preperation of NetMHC I and II data for prediction and manage 
#  predictions for RIAZ data 
# 1. Import data 
# 2. Output somatic mutations (sm) for NetMHC prediction
# 3. Output sbatch code (emlinating HLA not valid for netmhc)
# Note: Using NetMHCI-4.0 and NetMHCIIpan-3.2
#------------------------------------------------------

library(magrittr)
library(Biostrings)

#--------------------------------------------------------------------
# Step 1: Import Data 
#--------------------------------------------------------------------

#----- somatic mutation data 

load("riaz_mutdata3_sm.RData")
ls()

dim(riaz_mutdata3_sm)
riaz_mutdata3_sm[1:2,]

#----- patient with HLA prediction

riaz_pt_with_hla = read.delim(file="riaz_patientID_with_HLA.txt", 
                              sep= "\t", as.is=TRUE)
dim(riaz_pt_with_hla)
riaz_pt_with_hla[1:2,]

#----- valid hla alleles for netmhci and netmhcii

# HLA-i
valid_hlai = read.delim(file= "../NetMHCI_valid_HLA.txt", 
                        sep= "\t", header=FALSE, as.is=TRUE, 
                        col.names =c("terminal", "standard", "hla_type"))
dim(valid_hlai)
valid_hlai[1:2,]

# modify format of terminal 
valid_hlai$terminal1 = substr(valid_hlai$terminal, 5,9)

# HLA-ii
valid_hlaii = read.delim(file= "../NetMHCIIpan_valid_HLA.txt", 
                         sep= "\t", as.is=TRUE, header=FALSE)
dim(valid_hlaii)
valid_hlaii[1:2,]

# modify format
for(i in 1:nrow(valid_hlaii)){
  if(substr(valid_hlaii$V1[i],1,3)=="DRB"){
    valid_hlaii$terminal[i] = gsub("[:]", "", valid_hlaii$V1[i])
    valid_hlaii$terminal[i] = gsub("[\\*]", "_", valid_hlaii$terminal[i])
  } else{
    valid_hlaii$terminal[i] = gsub("[:, \\*]", "", valid_hlaii$V1[i])
  }
}

dim(valid_hlaii)
valid_hlaii[1:2,]

#--------------------------------------------------------------------
# Step 2: Output somatic mutations for predictions (in terminal)
#--------------------------------------------------------------------

riaz_sm = riaz_mutdata3_sm

sum(riaz_sm$som_mut9 == "ENST not found")
sum(riaz_sm$som_mut9 == "Sequence unavailable")

head(riaz_sm[riaz_sm$som_mut9 == "ENST not found",])
head(riaz_sm[riaz_sm$som_mut9 == "Sequence unavailable",])

# remove observations with missing sequence 

riaz_sm1 = riaz_sm[which(riaz_sm$som_mut9!="ENST not found" &
                          riaz_sm$som_mut9!= "Sequence unavailable" ),]
dim(riaz_sm1)
riaz_sm1[1:2,]

#-------calcualte number of peptides per patient

tb1 = table(riaz_sm1$id)
sort(tb1)

riaz_obs_per_pt1 = data.frame(tb1, stringsAsFactors=FALSE)
names(riaz_obs_per_pt1) = c("matchID", "nMuts_with_peptides")
dim(riaz_obs_per_pt1)
riaz_obs_per_pt1[1:2,]

riaz_obs_per_pt1$matchID = as.character(riaz_obs_per_pt1$matchID)

write.table(riaz_obs_per_pt1, file = "riaz_nMut_with_peptides.txt", 
            quote=FALSE, col.names= TRUE, row.names = FALSE)

# comparing with previous Laura results 
nrow(riaz_sm[riaz_sm$som_mut9 == "ENST not found" & riaz_sm$id=="Pt102_on",])
nrow(riaz_sm[riaz_sm$som_mut9 == "Sequence unavailable" & riaz_sm$id=="Pt102_on",])
riaz_obs_per_pt1[1,]
a= riaz_sm[riaz_sm$som_mut9 == "ENST not found",]
b = riaz_sm[riaz_sm$som_mut9 == "Sequence unavailable",]
unique(unique(a$id), unique(b$id))
unique(b$id)

#------- function to write sm file for netmhci and netmhcii-pan

# if start is numeric, the past function will add some space
riaz_sm1$start = as.character(riaz_sm1$start)
key.cols       = c("seqnames", "start", "REF", "ALT")
riaz_sm1$key   = apply(riaz_sm1[,key.cols], 1, paste, collapse=":")
riaz_sm1$start = as.numeric(riaz_sm1$start)
riaz_sm1$key   = gsub("^chr", "", riaz_sm1$key)
dim(riaz_sm1)
riaz_sm1[1:2,]

# prepare fasta files
for( i in 1:nrow(riaz_obs_per_pt1)){
  
  matchid_string = toString(riaz_obs_per_pt1$matchID[i])
  
  pt = riaz_sm1[which(riaz_sm1$id == matchid_string), ]
  
  # for HLA-I need to be at least length 9 
  pt9 = pt[nchar(pt$som_mut9)>=9,]
  
  fnm9 = paste("for_netmhci/", matchid_string, "_sm9.txt", sep="")
  cat("", file=fnm9)
  for(j in 1:nrow(pt9)){
    cat(paste(">", pt9$key[j], "\n", pt9$som_mut9[j],"\n", sep=""), 
        file=fnm9, append=TRUE)
  }

  fnm9r = paste("for_netmhci/", matchid_string, "_sm9r.txt", sep="")
  cat("", file=fnm9r)
  for(j in 1:nrow(pt9)){
    cat(paste(">", pt9$key[j], "\n", pt9$som_mut9r[j],"\n", sep=""), 
        file=fnm9r, append=TRUE)
  }

  # for HLA-II need to be at least length 15
  pt15 = pt[nchar(pt$som_mut15)>=15,]
  
  fnm15 = paste("for_netmhciipan/", matchid_string, "_sm15.txt", sep="")
  cat("", file=fnm15)
  for(j in 1:nrow(pt15)){
    cat(paste(">", pt15$key[j], "\n", pt15$som_mut15[j],"\n", sep=""),
        file=fnm15, append=TRUE)
  }

  fnm15r = paste("for_netmhciipan/", matchid_string, "_sm15r.txt", sep="")
  cat("", file=fnm15r)
  for(j in 1:nrow(pt15)){
    cat(paste(">", pt15$key[j], "\n", pt15$som_mut15r[j],"\n", sep=""), 
        file=fnm15r, append=TRUE)
  }
}

#--------------------------------------------------------------------
# Step 3: Output sbatch code for terminal
#--------------------------------------------------------------------

#------------HLA-I
# check how many HLA alleles can be predicted 

HLA1  = as.matrix(riaz_pt_with_hla[,7:12])
dim(HLA1)
HLA1[1:2,]

avail = t(apply(HLA1, 1, function(v){ v %in% valid_hlai$terminal1 }))
dim(avail)
avail[1:2,]

sum(c(avail))/(6*nrow(avail))
table(rowSums(avail))

# create shell script to submit jobs 

fnm = "for_netmhci/_NetMHCI_pred_validHLA.sh"
cat("", file=fnm)

for(i in 1:nrow(riaz_obs_per_pt1)){
  id = riaz_obs_per_pt1[i,1]
  t  = riaz_pt_with_hla[which(riaz_pt_with_hla$SampleName == id),7:12] 
  t1 = t[, t %in% valid_hlai$terminal1]
  t2 = paste("HLA-", t1, sep="")
  t3 = paste(t2, collapse=",")
  
  cat(paste("../runNetMHC_riaz.sh", id, t3, "\n", sep=" "), 
      file=fnm, append=TRUE)
  
  # cat(paste("sbatch ../runNetMHC_riaz_sbatch.sh", id, t3, "\n", sep=" "), 
  #     file=fnm, append=TRUE)
}

#------------HLA-II
# check how many HLA alleles can be predicted 

HLA2  = as.matrix(riaz_pt_with_hla[,13:ncol(riaz_pt_with_hla)])
dim(HLA2)
HLA2[1:2,]

avail = t(apply(HLA2, 1, function(v){ v %in% valid_hlaii$terminal }))
dim(avail)
avail[1:2,]

sum(c(avail))/(10*nrow(avail))
table(rowSums(avail))

#../netMHCIIpan -f example.fsa -a DRB1_0101 > example.fsa.myout

# create shell script to submit jobs 
# one job per 20 samples

jobid = 1

for(i in 1:nrow(riaz_obs_per_pt1)){
  
  if(i %% 20 == 1){
    fnm = sprintf("for_netmhciipan/_netMHCIIpan_pred_validHLA_%d.sh", jobid)
    cat("", file=fnm)
    jobid = jobid + 1
  }
  
  id = riaz_obs_per_pt1[i,1]
  wi = which(riaz_pt_with_hla$SampleName == id)
  t  = riaz_pt_with_hla[wi,13:ncol(riaz_pt_with_hla)]
  t1 = t[, t %in% valid_hlaii$terminal]
  
  t2 = paste(t1$DRB11, t1$DRB12, sep=",")
  t3 = paste("HLA", t1$DQA11, t1$DQB11, sep="-")
  t4 = paste("HLA", t1$DQA12, t1$DQB12, sep="-")
  t5 = paste("HLA", t1$DQA11, t1$DQB12, sep="-")
  t6 = paste("HLA", t1$DQA12, t1$DQB11, sep="-")
  t7 = paste("HLA", t1$DPA11, t1$DPB11, sep="-")
  t8 = paste("HLA", t1$DPA12, t1$DPB12, sep="-") 
  t9 = paste("HLA", t1$DPA11, t1$DPB12, sep="-")
  t10 = paste("HLA", t1$DPA12, t1$DPB11, sep="-")
  t11 = c(t2, t3, t4, t5, t6, t7, t8, t9, t10)
  w2rm = union(grep("--", t11, fixed=TRUE), grep("-$", t11))
  if(length(w2rm) > 0){ t11 = t11[-w2rm] }
  t11 = paste(t11, collapse=",")
 
  cat(paste("../runNetMHCIIpan_riaz.sh", id , t11, "\n", sep=" "),
      file=fnm, append=TRUE)
  
  # cat(paste("sbatch ../runNetMHCIIpan_riaz_sbatch.sh", id , t11, "\n", sep=" "),
  #     file=fnm, append=TRUE)
  
}

sessionInfo()
q(save="no")



