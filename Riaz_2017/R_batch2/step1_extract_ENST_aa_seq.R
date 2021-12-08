
#--------------------------------------------------------------------
# Extract amino acid (aa) sequence using ENST 
# Steps: 
  #1. Import data
  #2. Extract aa reference sequence using biomaRt
  #3. Save riaz_mutdata_refaa.RData
#--------------------------------------------------------------------

library(biomaRt)
library(Biostrings)
library(GenomicRanges)
library(stringr)
library(plyr)

#--------------------------------------------------------------------
# Step 1: Import Data and calculate start position of aa
#--------------------------------------------------------------------

load("../data/riaz_mutdata.RData")
ls()

dim(riaz_mutdata)
riaz_mutdata[1:2,c(1:10,97:104)]

#--------------------------------------------------------------------
# Step 2. Get aa sequence based on ENST
#--------------------------------------------------------------------

# since the mutations were annotated at 2019 using ANNOVAR
# we should use an archived version of ensembl
# here we check trnascript ENST00000405322 to check
# its annotated wild type aa is "S"
# and also transcript ENST00000304418
# its annotated mutation is at position 321

urls  = listEnsemblArchives()$url[7:20]
urls
wt_aa1 = rep("", length(urls))
len_aa = rep(NA, length(urls))

for(i in 1:length(urls)){
  cat(i, date(), "\n")
  u1 = urls[i]
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", 
                    host=u1)
  seq_i = getSequence(id=c("ENST00000405322", "ENST00000304418"), 
                      type="ensembl_transcript_id",
                      seqType="peptide", mart=ensembl)
  len_aa[i] = nchar(seq_i[1,1])
  wt_aa1[i] = substr(seq_i[2,1], 238, 238)
}

cbind(wt_aa1, len_aa, urls)

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", 
                  host="https://may2017.archive.ensembl.org")

uniq.EnsembleID = unique(riaz_mutdata$EnsembleID)
length(riaz_mutdata$EnsembleID)
length(uniq.EnsembleID)

# it seems if we try to get all the sequence at once, 
# there may be error for "Operation timed out"

starts = seq(1, length(uniq.EnsembleID), by=1000)
ends = c(seq(1000, length(uniq.EnsembleID), by=1000), length(uniq.EnsembleID))

cbind(starts, ends)

seq.prot0 = NULL

for(i in 1:length(starts)){
  cat(i, date(), "\n")
  seq_i = getSequence(id=uniq.EnsembleID[starts[i]:ends[i]], 
                      type="ensembl_transcript_id", 
                      seqType="peptide", mart=ensembl)
  seq.prot0 = rbind(seq.prot0, seq_i)
}

dim(seq.prot0)
seq.prot0[1:2,]

tb1 = table(seq.prot0$peptide)
table(tb1)

tb2 = tb1[tb1 > 1]
names(tb2) = substr(names(tb2), 1, 20)
tb2

g1 = grep("X", seq.prot0$peptide)
seq.prot0$peptide[g1]
seq.prot0$peptide[g1] = "Sequence unavailable"

table(seq.prot0$peptide=="Sequence unavailable")
setequal(uniq.EnsembleID, seq.prot0$ensembl_transcript_id)

#--------------------------------------------------------------------
# Step 3. Save dataset
#--------------------------------------------------------------------

riaz_mutdata_refaa = riaz_mutdata[,c(1, 9:13, 28:29, 98:104)] 
riaz_mutdata_refaa$prot.seq = rep("", nrow(riaz_mutdata_refaa))

mat1 = match(riaz_mutdata_refaa$EnsembleID, seq.prot0$ensembl_transcript_id)
table(is.na(mat1))

w2update = which(!is.na(mat1))
table(riaz_mutdata_refaa$EnsembleID[w2update] == 
        seq.prot0$ensembl_transcript_id[mat1[w2update]])
riaz_mutdata_refaa$prot.seq[w2update] = seq.prot0$peptide[mat1[w2update]]

dim(riaz_mutdata_refaa)
riaz_mutdata_refaa[1:2,]
table(riaz_mutdata_refaa$prot.seq == "")

save(riaz_mutdata_refaa, file = "../data/riaz_mutdata_refaa.RData")

sessionInfo()
q(save="no")

