
#--------------------------------------------------------------------
# Create somatic mutation chains 
# 1. Import data
# 2. extrac peptides that cover the mutated amino acid and is 17 or 
#    35 aa long. 17 aa long sequence is for HLA-I (9 aa peptide)
#    and 35 aa long sequence is for HLA-II (15 aa peptides + 
#    3aa context on either side)
# 3. save riaz_mutdata_sm.RData
#--------------------------------------------------------------------

library(stringr)

#--------------------------------------------------------------------
# Step 1: Import Data 
#--------------------------------------------------------------------

# input dataset with the reference seq

load("../data/riaz_mutdata_refaa.RData")
dim(riaz_mutdata_refaa)
riaz_mutdata_refaa[1,]

# check stop codons, and remove them

fun1 <- function(x){substr(x, nchar(x), nchar(x))}
last_aa = sapply(riaz_mutdata_refaa$prot.seq, fun1)
table(last_aa)

prot.seq.nostop = gsub("\\*.*", "", riaz_mutdata_refaa$prot.seq)
riaz_mutdata_refaa$prot.seq[1]
prot.seq.nostop[1]

n0 = nchar(riaz_mutdata_refaa$prot.seq)
n1 = nchar(prot.seq.nostop)
table(n0 - n1)

if(any(n0 - n1) > 1){ stop("stop condon in the middle of sequence\n") }

# it looks like those transcripts without stop codon are immunoglobulins
riaz_mutdata_refaa$prot.seq = prot.seq.nostop
riaz_mutdata_refaa[which(n0 - n1==0)[1:2],]

#--------------------------------------------------------------------
# Step 2: Replace mutated AA and extract peptides for prediction
#--------------------------------------------------------------------

som_mut17  = som_mut17r = rep("", nrow(riaz_mutdata_refaa))
som_mut35  = som_mut35r = rep("", nrow(riaz_mutdata_refaa))
mut17_start = mut17_end = rep(NA, nrow(riaz_mutdata_refaa))
mut35_start = mut35_end = rep(NA, nrow(riaz_mutdata_refaa))

for(i in 1:nrow(riaz_mutdata_refaa)){
  b = riaz_mutdata_refaa[i,]
  
  if(b$prot.seq == ""){ 
    som_mut17[i] = som_mut17r[i] = "ENST not found"
    som_mut35[i] = som_mut35r[i] = "ENST not found"
  } else if (b$prot.seq=="Sequence unavailable"){
    som_mut17[i] = som_mut17r[i] = "Sequence unavailable"
    som_mut35[i] = som_mut35r[i] = "Sequence unavailable"
  } else{

    if(substr(b$prot.seq, b$pos, b$pos) != b$W){
      stop("mismatch of amino acid for wild type allelel")
    }
    
    # replace aa at position 'pos' with mutated aa
    b$prot.seq.ref = b$prot.seq
    substr(b$prot.seq, as.numeric(b$pos), as.numeric(b$pos)) = as.character(b$M)
    
    # extract aa of interest per somatic mutation
    mut17_start[i] = max(as.numeric(b$pos)-8, 1)
    mut17_end[i]   = min(as.numeric(b$pos)+8, nchar(b$prot.seq))
    som_mut17[i]   = substr(b$prot.seq, mut17_start[i], mut17_end[i])
    som_mut17r[i]  = substr(b$prot.seq.ref, mut17_start[i], mut17_end[i])
    
    mut35_start[i] = max(as.numeric(b$pos)-17, 1)
    mut35_end[i]   = min(as.numeric(b$pos)+17, nchar(b$prot.seq))
    som_mut35[i]   = substr(b$prot.seq, mut35_start[i], mut35_end[i])
    som_mut35r[i]  = substr(b$prot.seq.ref, mut35_start[i], mut35_end[i])
  }
}

riaz_mutdata_sm = data.frame(riaz_mutdata_refaa[,1:15], som_mut17, som_mut17r, 
                         som_mut35, som_mut35r, mut17_start, mut17_end, 
                         mut35_start, mut35_end, stringsAsFactors=FALSE)
dim(riaz_mutdata_sm)
riaz_mutdata_sm[1:2,]

table(riaz_mutdata_sm$som_mut17 == "Sequence unavailable")
table(riaz_mutdata_sm$som_mut17 == "ENST not found")

w1 = which(riaz_mutdata_sm$som_mut17 == "Sequence unavailable")
riaz_mutdata_sm[w1,]

riaz_mutdata_sm = riaz_mutdata_sm[-w1,]
dim(riaz_mutdata_sm)
riaz_mutdata_sm[1:2,]

save(riaz_mutdata_sm, file = "../data/riaz_mutdata_sm.RData")

sessionInfo()
q(save="no")

