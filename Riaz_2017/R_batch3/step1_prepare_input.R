
#--------------------------------------------------------------------
# Preperation Input data for Peppermint 
# 1. Import data 
# 2. Output somatic mutations (sm) for NetMHC prediction
# 3. Output sbatch code (emlinating HLA that are not valid 
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

# write out input for PEPPRMINT

# if 'start' is numeric, the past function will add some space
riaz_sm1$start = as.character(riaz_sm1$start)
key.cols       = c("seqnames", "start", "REF", "ALT")
riaz_sm1$key   = apply(riaz_sm1[,key.cols], 1, paste, collapse=":")
riaz_sm1$start = as.numeric(riaz_sm1$start)
riaz_sm1$key   = gsub("^chr", "", riaz_sm1$key)
dim(riaz_sm1)

riaz_sm1$mut17_len = nchar(riaz_sm1$som_mut17)
table(riaz_sm1$mut17_len)
riaz_sm1 = riaz_sm1[which(riaz_sm1$mut17_len >= 9),]
dim(riaz_sm1)
riaz_sm1[1:2,]

# prepare input files for PEPPRMINT

idx2kp = NULL
seq_mut = seq_ref = NULL

encode9 <- function(x){
  paste0(substr(x, 1, 4), "XXX", substr(x, 5, 5), "XXX", substr(x, 6, 9))
}

for( i in 1:nrow(riaz_sm1)){
  if(i %% 2000 == 0) {cat(i, date(), "\n")}
  len_i = riaz_sm1$mut17_len[i]
  
  for(j in 1:(len_i - 8)){
    seq_ij_mut = encode9(substr(riaz_sm1$som_mut17[i], j, j+8))
    seq_ij_ref = encode9(substr(riaz_sm1$som_mut17r[i], j, j+8))
    
    seq_mut = c(seq_mut, seq_ij_mut)
    seq_ref = c(seq_ref, seq_ij_ref)
    idx2kp  = c(idx2kp, i)
  }
}

pdat_mut = data.frame(peptide = seq_mut)
pdat_mut = cbind(pdat_mut, riaz_sm1[idx2kp, c("key", "id")])

pdat_ref = data.frame(peptide = seq_ref)
pdat_ref = cbind(pdat_ref, riaz_sm1[idx2kp, c("key", "id")])

names(pdat_mut)[3] = names(pdat_ref)[3] = "sample"

dim(pdat_mut)
pdat_mut[1:2,]

dim(pdat_ref)
pdat_ref[1:2,]

u_sample = unique(pdat_mut$sample)

fwrite(pdat_mut, "../data/riaz_peptide_mut.txt", sep="\t")
fwrite(pdat_ref, "../data/riaz_peptide_ref.txt", sep="\t")

#--------------------------------------------------------------------
# Step 3: Output sbatch code for terminal
#--------------------------------------------------------------------

# patient with HLA prediction

riaz_pt_with_hla = fread("../data/riaz_patientID_with_HLA.txt")
dim(riaz_pt_with_hla)
riaz_pt_with_hla[1:2,]
riaz_pt_with_hla = as.data.frame(riaz_pt_with_hla)

table(u_sample %in% riaz_pt_with_hla$SampleName)

# valid hla alleles

# HLA-i
valid_hlai = scan("../data/netMHC_pan_4.1_allele_names.txt", 
                  what = character())
length(valid_hlai)
valid_hlai[1:2]

# HLA-I
# check how many HLA alleles can be predicted 

HLA1  = as.matrix(riaz_pt_with_hla[,7:12])
dim(HLA1)
HLA1[1:2,]

avail = t(apply(HLA1, 1, function(v){ v %in% valid_hlai}))
dim(avail)
avail[1:2,]

sum(c(avail))/(6*nrow(avail))
table(rowSums(avail))

# create hla allele list file

hlas = rep(NA, length(u_sample))

for(i in 1:length(u_sample)){
  wi = which(riaz_pt_with_hla$SampleName == u_sample[i])
  if(length(wi) != 1){ stop('unexpected length') }
  t0 = riaz_pt_with_hla[wi,7:12]
  t0 = unique(unlist(t0))
  t1 = sort(t0[t0 %in% valid_hlai])
  t2 = paste(t1, collapse=",")
  hlas[i] = t2
}

hla_df = data.frame(sample=u_sample, hlas=hlas)
dim(hla_df)
hla_df[1:2,]

fnm = "../data/riaz_HLA_I_allele_list.txt"
fwrite(hla_df, file=fnm, sep=" ", col.names = FALSE)


sessionInfo()
q(save="no")



