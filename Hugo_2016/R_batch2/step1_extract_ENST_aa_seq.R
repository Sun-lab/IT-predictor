#--------------------------------------------------------------------
# Extract amino acid (aa) sequence using ENST 
# Steps: 
  #1. Import data
  #2. Extract aa reference sequence using biomaRt
  #3. Save mutdata_refaa.RData
#--------------------------------------------------------------------

library("biomaRt")
library("Biostrings")
library("GenomicRanges")
library("stringr")
library("plyr")

#--------------------------------------------------------------------
# Step 1: Import Data and calculate start position of aa
#--------------------------------------------------------------------

load("hugo_mutdata3.RData")
ls()

dim(hugo_mutdata3)
hugo_mutdata3[1:2,c(1:10,97:104)]

#--------------------------------------------------------------------
# Step 2. Get aa sequence based on ENST
#--------------------------------------------------------------------

# first try to get all the sequences in one batch

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

uniq.EnsembleID = unique(hugo_mutdata3$EnsembleID)
length(hugo_mutdata3$EnsembleID)
length(uniq.EnsembleID)

date()
seq.prot0 = getSequence(id=uniq.EnsembleID, 
                            type="ensembl_transcript_id", 
                            seqType="peptide", 
                            mart=ensembl)
date()

dim(seq.prot0)
seq.prot0[1:2,]
table(table(seq.prot0$peptide))

table(seq.prot0$peptide=="Sequence unavailable")

table(seq.prot0$ensembl_transcript_id %in% uniq.EnsembleID)
missing.ens = setdiff(uniq.EnsembleID, seq.prot0$ensembl_transcript_id)
length(missing.ens)

# next try to get the sequences one by one for thos missing ones

seq.prot1 = list()

for(i in 1:length(missing.ens)){
  seq.prot1[[i]] = getSequence(id=missing.ens[i],
                              type="ensembl_transcript_id",
                              seqType="peptide", 
                              mart=ensembl)
}

table(sapply(seq.prot1, nrow))

#--------------------------------------------------------------------
# Step 3. Save dataset
#--------------------------------------------------------------------

hugo_mutdata3_refaa = hugo_mutdata3[,c(1, 9:13, 28:29, 98:104)] 
hugo_mutdata3_refaa$prot.seq = rep("", nrow(hugo_mutdata3_refaa))

mat1 = match(hugo_mutdata3_refaa$EnsembleID, seq.prot0$ensembl_transcript_id)
table(is.na(mat1))
w2update = which(!is.na(mat1))
table(hugo_mutdata3_refaa$EnsembleID[w2update] == seq.prot0$ensembl_transcript_id[mat1[w2update]])
hugo_mutdata3_refaa$prot.seq[w2update] = seq.prot0$peptide[mat1[w2update]]

dim(hugo_mutdata3_refaa)
hugo_mutdata3_refaa[1:2,]
table(hugo_mutdata3_refaa$prot.seq == "")

save(hugo_mutdata3_refaa, file = "hugo_mutdata3_refaa.RData")

sessionInfo()
q(save="no")

