
#------------------------------------------------------
# Create somatic mutation chains 
# 1. Import data
# 2. Replace mutated AA and Extract somatic mutations for 9 and 15 peptide  
#    length, includes replacing * with X in stop codon
# 3. save hugo_mutdata3_sm.RData
#------------------------------------------------------

#--------------------------------------------------------------------
# Step 1: Import Data 
#--------------------------------------------------------------------
#input dataset with the reference seq

load("hugo_mutdata3_refaa.RData")
dim(hugo_mutdata3_refaa)
hugo_mutdata3_refaa[1,]

#--------------------------------------------------------------------
# Step 2: Replace mutated AA and extract sm chain for prediction
# includes replacing * with X for stop codon
#--------------------------------------------------------------------

som_mut9 = som_mut9r = som_mut15 = som_mut15r = rep("", nrow(hugo_mutdata3_refaa))

for(i in 1:nrow(hugo_mutdata3_refaa)){
  b = hugo_mutdata3_refaa[i,]
  
  if(b$prot.seq == ""){ 
    som_mut9[i]  = som_mut9r[i]  = "ENST not found"
    som_mut15[i] = som_mut15r[i] = "ENST not found"
  } else if (b$prot.seq=="Sequence unavailable"){
    som_mut9[i]  = som_mut9r[i]  = "Sequence unavailable"
    som_mut15[i] = som_mut15r[i] = "Sequence unavailable"
  } else{
    b$prot.seq = as.character(b$prot.seq)
    # replace * with X 
    b$prot.seq = gsub("\\*", "X", b$prot.seq)
    
    # replace aa at position 'pos' with mutated aa 'M'
    b$prot.seq.ref = b$prot.seq
    substr(b$prot.seq, as.numeric(b$pos), as.numeric(b$pos)) <- as.character(b$M)
    
    # extract aa of interest per somatic mutation
    # for 9 sliding window 
    som_mut9[i]  = substr(b$prot.seq, max(as.numeric(b$pos)-8, 1), 
                      min(as.numeric(b$pos)+8, nchar(b$prot.seq)))
    
    som_mut9r[i] = substr(b$prot.seq.ref, max(as.numeric(b$pos)-8, 1), 
                      min(as.numeric(b$pos)+8, nchar(b$prot.seq.ref)))
    
    # for 15 sliding window - want 14 aa upstream of pos and 14 aa downstream of pos
    som_mut15[i]  = substr(b$prot.seq, max(as.numeric(b$pos)-14, 1), 
                       min(as.numeric(b$pos)+14, nchar(b$prot.seq)))
    
    som_mut15r[i] = substr(b$prot.seq.ref, max(as.numeric(b$pos)-14, 1), 
                       min(as.numeric(b$pos)+14, nchar(b$prot.seq.ref)))
  }
}

hugo_mutdata3_sm = data.frame(hugo_mutdata3_refaa[,1:15], som_mut9, som_mut9r, 
                         som_mut15, som_mut15r, stringsAsFactors=FALSE)
dim(hugo_mutdata3_sm)
hugo_mutdata3_sm[1:2,]

save(hugo_mutdata3_sm, file = "hugo_mutdata3_sm.RData")

sessionInfo()
q(save="no")

