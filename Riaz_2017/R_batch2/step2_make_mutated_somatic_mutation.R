
#--------------------------------------------------------------------
# Create somatic mutation chains 
# 1. Import data
# 2. Replace mutated AA and Extract somatic mutations for 9 and 15   
#    amino-acid length, includes replacing * with X in stop codon
# 3. save riaz_mutdata_sm.RData
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Step 1: Import Data 
#--------------------------------------------------------------------
#input dataset with the reference seq

load("../data/riaz_mutdata_refaa.RData")
dim(riaz_mutdata_refaa)
riaz_mutdata_refaa[1,]

#--------------------------------------------------------------------
# Step 2: Replace mutated AA and extract sm chain for prediction
# includes replacing * with X for stop codon
#--------------------------------------------------------------------

som_mut9  = som_mut9r  = rep("", nrow(riaz_mutdata_refaa))
som_mut15 = som_mut15r = rep("", nrow(riaz_mutdata_refaa))

for(i in 1:nrow(riaz_mutdata_refaa)){
  b = riaz_mutdata_refaa[i,]
  
  if(b$prot.seq == ""){ 
    som_mut9[i]  = som_mut9r[i]  = "ENST not found"
    som_mut15[i] = som_mut15r[i] = "ENST not found"
  } else if (b$prot.seq=="Sequence unavailable"){
    som_mut9[i]  = som_mut9r[i]  = "Sequence unavailable"
    som_mut15[i] = som_mut15r[i] = "Sequence unavailable"
  } else{

    # replace aa at position 'pos' with mutated aa
    b$prot.seq.ref = b$prot.seq
    if(substr(b$prot.seq, b$pos, b$pos) != b$W){
      stop("mismatch of amino acid for wild type allelel")
    }
    
    substr(b$prot.seq, b$pos, b$pos) = b$M
    
    # extract aa of interest per somatic mutation
    # for 9 sliding window 
    som_mut9[i]  = substr(b$prot.seq, max(as.numeric(b$pos)-8, 1), 
                      min(as.numeric(b$pos)+8, nchar(b$prot.seq)))
    
    som_mut9r[i] = substr(b$prot.seq.ref, max(as.numeric(b$pos)-8, 1), 
                      min(as.numeric(b$pos)+8, nchar(b$prot.seq.ref)))
    
    # for 15 sliding window - want 14 aa upstream/downstream
    som_mut15[i]  = substr(b$prot.seq, max(as.numeric(b$pos)-14, 1), 
                       min(as.numeric(b$pos)+14, nchar(b$prot.seq)))
    
    som_mut15r[i] = substr(b$prot.seq.ref, max(as.numeric(b$pos)-14, 1), 
                       min(as.numeric(b$pos)+14, nchar(b$prot.seq.ref)))
  }
}

riaz_mutdata_sm = data.frame(riaz_mutdata_refaa[,1:15], som_mut9, som_mut9r, 
                         som_mut15, som_mut15r, stringsAsFactors=FALSE)
dim(riaz_mutdata_sm)
riaz_mutdata_sm[1:2,]

save(riaz_mutdata_sm, file = "riaz_mutdata_sm.RData")

sessionInfo()
q(save="no")

