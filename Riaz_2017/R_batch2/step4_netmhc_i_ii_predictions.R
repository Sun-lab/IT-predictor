
#--------------------------------------------------------------------
# Preperation of NetMHC I and II data for prediction and manage 
#  predictions for RIAZ data 
# 1. Import data 
# 2. Output somatic mutations (sm) for NetMHC prediction
# 3. Output sbatch code (emlinating HLA not valid for netmhc)
# Note: Using NetMHCpan4.1 and NetMHCIIpan-4.0
#--------------------------------------------------------------------

library(magrittr)
library(Biostrings)
library(data.table)

#--------------------------------------------------------------------
# Step 1: Import Data 
#--------------------------------------------------------------------

# somatic mutation data 

load("../data/riaz_mutdata_sm.RData")
ls()

dim(riaz_mutdata_sm)
riaz_mutdata_sm[1:2,]

#--------------------------------------------------------------------
# Step 2: Output somatic mutations for predictions
#--------------------------------------------------------------------

sum(riaz_mutdata_sm$som_mut17 == "ENST not found")
sum(riaz_mutdata_sm$som_mut17 == "Sequence unavailable")

head(riaz_mutdata_sm[riaz_mutdata_sm$som_mut17 == "ENST not found",])
head(riaz_mutdata_sm[riaz_mutdata_sm$som_mut17 == "Sequence unavailable",])

# remove observations with missing sequence 
w2rm = c("ENST not found", "Sequence unavailable")
riaz_sm1 = riaz_mutdata_sm[which(! riaz_mutdata_sm$som_mut17 %in% w2rm),]
dim(riaz_sm1)
riaz_sm1[1:2,]

# calcualte the number of peptides per patient

tb1 = table(riaz_sm1$id)
sort(tb1)

riaz_obs_per_pt1 = data.frame(tb1, stringsAsFactors=FALSE)
names(riaz_obs_per_pt1) = c("matchID", "nMuts_with_peptides")
dim(riaz_obs_per_pt1)
riaz_obs_per_pt1[1:2,]

riaz_obs_per_pt1$matchID = as.character(riaz_obs_per_pt1$matchID)

write.table(riaz_obs_per_pt1, file = "../data/riaz_nMut_with_peptides.txt", 
            quote=FALSE, col.names= TRUE, row.names = FALSE)

# write out fastq format input for netmhci-pan and netmhcii-pan

# if 'start' is numeric, the past function will add some space
riaz_sm1$start = as.character(riaz_sm1$start)
key.cols       = c("seqnames", "start", "REF", "ALT")
riaz_sm1$key   = apply(riaz_sm1[,key.cols], 1, paste, collapse=":")
riaz_sm1$start = as.numeric(riaz_sm1$start)
riaz_sm1$key   = gsub("^chr", "", riaz_sm1$key)
dim(riaz_sm1)
riaz_sm1[1:2,]

# prepare fasta files
for( i in 1:nrow(riaz_obs_per_pt1)){
  
  match_id = riaz_obs_per_pt1$matchID[i]
  pt       = riaz_sm1[which(riaz_sm1$id == match_id), ]
  
  # for HLA-I need to be at least length 9 
  pt_i = pt[which(nchar(pt$som_mut17)>=9),]
  
  fnm1 = paste("../data/netmhci/", match_id, "_mut.txt", sep="")
  cat("", file=fnm1)
  
  for(j in 1:nrow(pt_i)){
    cat(paste(">", pt_i$key[j], "\n", pt_i$som_mut17[j], "\n", sep=""), 
        file=fnm1, append=TRUE)
  }

  fnm1r = paste("../data/netmhci/", match_id, "_ref.txt", sep="")
  cat("", file=fnm1r)
  
  for(j in 1:nrow(pt_i)){
    cat(paste(">", pt_i$key[j], "\n", pt_i$som_mut17r[j],"\n", sep=""), 
        file=fnm1r, append=TRUE)
  }

  # for HLA-II need to be at least length 21
  pt_i = pt[which(nchar(pt$som_mut35)>=21),]
  
  fnm1 = paste("../data/netmhcii/", match_id, "_mut.txt", sep="")
  cat("", file=fnm1)
  
  for(j in 1:nrow(pt_i)){
    cat(paste(">", pt_i$key[j], "\n", pt_i$som_mut35[j],"\n", sep=""),
        file=fnm1, append=TRUE)
  }

  fnm1r = paste("../data/netmhcii/", match_id, "_ref.txt", sep="")
  cat("", file=fnm1r)
  
  for(j in 1:nrow(pt_i)){
    cat(paste(">", pt_i$key[j], "\n", pt_i$som_mut35r[j],"\n", sep=""), 
        file=fnm1r, append=TRUE)
  }
}

#--------------------------------------------------------------------
# Step 3: Output sbatch code for terminal
#--------------------------------------------------------------------

# patient with HLA prediction

riaz_pt_with_hla = fread("../data/riaz_patientID_with_HLA.txt")
dim(riaz_pt_with_hla)
riaz_pt_with_hla[1:2,]
riaz_pt_with_hla = as.data.frame(riaz_pt_with_hla)
  
# valid hla alleles for netmhci and netmhcii

# HLA-i
valid_hlai = scan("../data/netMHC_pan_4.1_allele_names.txt", 
                  what = character())
length(valid_hlai)
valid_hlai[1:2]

# HLA-ii
valid_hlaii = fread("../data/netMHC_II_pan_4.0_allele_names.list", 
                    fill=TRUE)
dim(valid_hlaii)
valid_hlaii[c(1:2, 101:102),]

valid_hlaii$DR = gsub("[:]", "", valid_hlaii$DR)
valid_hlaii$DR = gsub("[\\*]", "_", valid_hlaii$DR)

for(i in 2:5){
  valid_hlaii[[i]] = gsub("[:, \\*]", "", valid_hlaii[[i]])
}
dim(valid_hlaii)
valid_hlaii[c(1:2, 101:102),]

valid_hlaii = unlist(valid_hlaii)
length(valid_hlaii)

# HLA-I
# check how many HLA alleles can be predicted 

table(riaz_obs_per_pt1$matchID  %in% riaz_pt_with_hla$SampleName)
mat1 = match(riaz_obs_per_pt1$matchID, riaz_pt_with_hla$SampleName)
riaz_pt_with_hla = riaz_pt_with_hla[mat1,]
dim(riaz_pt_with_hla)

HLA1  = as.matrix(riaz_pt_with_hla[,7:12])
dim(HLA1)
HLA1[1:2,]

avail = t(apply(HLA1, 1, function(v){ v %in% valid_hlai}))
dim(avail)
avail[1:2,]

sum(c(avail))/(6*nrow(avail))
table(rowSums(avail))

nHLA = t(apply(HLA1, 1, function(v){ length(unique(v)) }))
table(nHLA)

nrow(unique(HLA1))

# create shell script to submit jobs 

fnm = "step4_netMHC_pan_submit_jobs.sh"
cat("", file=fnm)

for(i in 1:nrow(riaz_obs_per_pt1)){
  id = riaz_obs_per_pt1[i,1]
  t0 = riaz_pt_with_hla[which(riaz_pt_with_hla$SampleName == id),7:12]
  t0 = unique(unlist(t0))
  t1 = sort(t0[t0 %in% valid_hlai])
  t2 = paste(t1, collapse=",")
  
  cat(paste("./step4_run_NetMHCpan.sh", id, t2, "\n", sep=" "), 
      file=fnm, append=TRUE)
}

# HLA-II
# check how many HLA alleles can be predicted 

HLA2  = as.matrix(riaz_pt_with_hla[,13:ncol(riaz_pt_with_hla)])
dim(HLA2)
HLA2[1:2,]

avail = t(apply(HLA2, 1, function(v){ v %in% valid_hlaii }))
dim(avail)
avail[1:2,]

sum(c(avail))/(10*nrow(avail))
table(rowSums(avail))

#../netMHCIIpan -f example.fsa -a DRB1_0101 > example.fsa.myout

# create shell script to submit jobs 
# one job per 15 samples

jobid = 1

for(i in 1:nrow(riaz_obs_per_pt1)){
  
  if(i %% 15 == 1){
    fnm = sprintf("step4_netMHC_II_pan_submit_job_%d.sh", jobid)
    cat("", file=fnm)
    jobid = jobid + 1
  }
  
  id = riaz_obs_per_pt1[i,1]
  wi = which(riaz_pt_with_hla$SampleName == id)
  t0 = riaz_pt_with_hla[wi,13:ncol(riaz_pt_with_hla)]
  
  t1 = list()
  t1[["DRB1"]] = sort(unique(c(t0$DRB11, t0$DRB12)))
  t1[["DPA1"]] = sort(unique(c(t0$DPA11, t0$DPA12)))
  t1[["DPB1"]] = sort(unique(c(t0$DPB11, t0$DPB12)))
  t1[["DQA1"]] = sort(unique(c(t0$DQA11, t0$DQA12)))
  t1[["DQB1"]] = sort(unique(c(t0$DQB11, t0$DQB12)))
  
  for(i in 1:length(t1)){
    t1[[i]] = intersect(t1[[i]], valid_hlaii)
  }
  
  tR = paste(t1[["DRB1"]], collapse=",")
  tP = tQ = NULL
  
  for(alpha in t1[["DPA1"]]){
    for(beta in t1[["DPB1"]]){
      tP = c(tP, paste("HLA", alpha, beta, sep="-"))
    }
  }
  
  for(alpha in t1[["DQA1"]]){
    for(beta in t1[["DQB1"]]){
      tQ = c(tQ, paste("HLA", alpha, beta, sep="-"))
    }
  }
  
  hlaII = paste(c(tR, tP, tQ), collapse=",")
 
  cat(paste("./step4_run_NetMHCIIpan.sh", id , hlaII, "\n", sep=" "),
      file=fnm, append=TRUE)
  
}

sessionInfo()
q(save="no")



