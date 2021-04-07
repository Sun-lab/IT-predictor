
library(VariantAnnotation)

# ------------------------------------------------------------------------
# read in sample information
# ------------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

info1 = read.table("SraRunTable.txt", sep="\t", as.is=TRUE, header=TRUE)
dim(info1)
info1[1:2,]

info2 = read.table("SraRunTable_Vanderbilt.txt", sep="\t", as.is=TRUE, 
                   header=TRUE)
dim(info2)
info2[1:2,]

ww1 = which(info1$tissue == "melanoma biopsy")
subjects1 = info1$Sample_Name[ww1]
subjects1 = gsub("-baseline", "", subjects1)
subjects1

ww2 = which(info2$tissue != "PBMC")
subjects = info2$Library_Name[ww2]
subjects = gsub("-PD1Cell2016-WES", "", subjects)
subjects = gsub("-tumor", "", subjects)
subjects = gsub("-", "_", subjects)
subjects = paste0("Vand_", subjects)
subjects

subjects = c(subjects1, subjects)
length(subjects)

# ----------------------------------------------------------------------
# load data of strelka, check the number of called mutations
# ----------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016/strelka_results")

snv.vcfs = paste0(subjects, "_passed.somatic.snvs.vcf")
length(snv.vcfs)
snv.vcfs[1:2]

indel.vcfs = paste0(subjects, "_passed.somatic.indels.vcf")
length(indel.vcfs)
indel.vcfs[1:2]

snvs = indels = list()

for(i in 1:length(subjects)){
  id1  = subjects[i]

  if(! file.exists(snv.vcfs[i])){
    stop(sprintf("file %s does not exist\n", snv.vcfs[i]))
  }
  
  if(! file.exists(indel.vcfs[i])){
    stop(sprintf("file %s does not exist\n", indel.vcfs[i]))
  }
  
  snv1  = readVcf(snv.vcfs[i], "hg38")
  snvs[[id1]] = snv1
  
  indel1  = readVcf(indel.vcfs[i], "hg38")
  indels[[id1]] = indel1
}

n.snvs = sapply(snvs, nrow)
summary(n.snvs)

n.indels = sapply(indels, nrow)
summary(n.indels)

# ----------------------------------------------------------------------
# read in mutation call data
# ----------------------------------------------------------------------

for(kk in 1:length(subjects)){
  
  id1 = subjects[kk]
  cat("\n", kk, id1, "\n")
  
  # ----------------------------------------------------------------------
  # first work on SNVs
  # ----------------------------------------------------------------------
  
  rBT = rowRanges(snvs[[id1]])
  iBT = info(snvs[[id1]])
  gBT = geno(snvs[[id1]])
  
  # --------------------------------------------------------------------
  # obtain read counts for reference or alternative allele
  # --------------------------------------------------------------------

  rBT = as.data.frame(rBT)
  ref = as.character(rBT$REF)
  alt = sapply(rBT$ALT, function(v1){paste(as.character(v1), collapse=":")})

  tier   = iBT$TQSS_NT
  alCt   = list()
  
  alCt[["A"]] = gBT$AU
  alCt[["C"]] = gBT$CU
  alCt[["G"]] = gBT$GU
  alCt[["T"]] = gBT$TU
  
  dat2kp = as.data.frame(iBT)
  dat2kp$seqnames = rBT$seqnames
  dat2kp$start    = rBT$start
  dat2kp$end      = rBT$end
  
  dat2kp$REF = ref
  dat2kp$ALT = alt
  
  dat2kp$nRefTumor  = dat2kp$nAltTumor  = rep(NA, nrow(dat2kp))
  dat2kp$nRefNormal = dat2kp$nAltNormal = rep(NA, nrow(dat2kp))

  for(ii in 1:nrow(dat2kp)){
    dat2kp$nRefNormal[ii] = alCt[[ref[ii]]][ii, 1, tier[ii]]
    dat2kp$nRefTumor[ii]  = alCt[[ref[ii]]][ii, 2, tier[ii]]
    
    if(! alt[ii] %in% c("A", "C", "G", "T")){
      next
    }

    dat2kp$nAltNormal[ii] = alCt[[alt[ii]]][ii, 1, tier[ii]]
    dat2kp$nAltTumor[ii]  = alCt[[alt[ii]]][ii, 2, tier[ii]]
  }
  
  DP1 = gBT$DP - gBT$FDP
  dim(DP1)
  DP1[1:2,]
  
  dat2kp$rdNormalFiltered = DP1[,1]
  dat2kp$rdTumorFiltered  = DP1[,2]

  dat2kp.snv = dat2kp
  
  # ----------------------------------------------------------------------
  # next work on indels
  # ----------------------------------------------------------------------
  
  rBT = rowRanges(indels[[id1]])
  iBT = info(indels[[id1]])
  gBT = geno(indels[[id1]])
  
  # --------------------------------------------------------------------
  # obtain read counts for reference or alternative allele
  # --------------------------------------------------------------------
  
  rBT = as.data.frame(rBT)
  
  if(nrow(rBT) == 0){
    dat2kp.indel = NULL
  }else{
    ref = as.character(rBT$REF)
    alt = sapply(rBT$ALT, function(v1){paste(as.character(v1), collapse=":")})
    
    tier = iBT$TQSI
    
    dat2kp          = as.data.frame(iBT)
    dat2kp$seqnames = rBT$seqnames
    dat2kp$start    = rBT$start
    dat2kp$end      = rBT$end
    
    dat2kp$REF = ref
    dat2kp$ALT = alt
    
    dat2kp$nRefTumor  = dat2kp$nAltTumor  = rep(NA, nrow(dat2kp))
    dat2kp$nRefNormal = dat2kp$nAltNormal = rep(NA, nrow(dat2kp))
    
    for(ii in 1:nrow(dat2kp)){
      
      dat2kp$nRefNormal[ii] = gBT$TAR[ii, 1, tier[ii]]
      dat2kp$nRefTumor[ii]  = gBT$TAR[ii, 2, tier[ii]]
      
      dat2kp$nAltNormal[ii] = gBT$TIR[ii, 1, tier[ii]]
      dat2kp$nAltTumor[ii]  = gBT$TIR[ii, 2, tier[ii]]
      
    }
    
    DP1 = gBT$DP - gBT$FDP
    dim(DP1)
    
    dat2kp$rdNormalFiltered = DP1[,1]
    dat2kp$rdTumorFiltered  = DP1[,2]
    
    dat2kp.indel = dat2kp
  }

  # --------------------------------------------------------------------
  # combine SNV and indel information
  # --------------------------------------------------------------------
  
  if(is.null(dat2kp.indel)){
    dat2kp = dat2kp.snv
    dat2kp$type = rep("SNV", times=nrow(dat2kp))
    
    print(dim(dat2kp.snv))
    print(dat2kp.snv[1,,FALSE])
    
  }else{
    print("dim(dat2kp.snv)")
    print(dim(dat2kp.snv))
    print("dim(dat2kp.indel)")
    print(dim(dat2kp.indel))
    
    print("dat2kp.snv[1,,FALSE]")
    print(dat2kp.snv[1,,FALSE])
    
    print("dat2kp.indel[1,,FALSE]")
    print(dat2kp.indel[1,,FALSE])

    table(dat2kp.snv$SOMATIC)
    table(dat2kp.indel$SOMATIC)
    table(dat2kp.indel$SVTYPE, useNA="ifany")
    table(dat2kp.indel$OVERLAP, useNA="ifany")
    
    col2rm = c("RU", "RC", "IC", "IHP", "SVTYPE", "OVERLAP")
    dat2kp.indel = dat2kp.indel[,-which(names(dat2kp.indel) %in% col2rm)]
    
    dim(dat2kp.snv)
    dim(dat2kp.indel)
    dat2kp.snv[1,,FALSE]
    dat2kp.indel[1,,FALSE]
    
    print(cbind(names(dat2kp.snv), names(dat2kp.indel)))
    
    names(dat2kp.indel) = names(dat2kp.snv)
    dat2kp = rbind(dat2kp.snv, dat2kp.indel)
    
    ntype = c(nrow(dat2kp.snv), nrow(dat2kp.indel))
    dat2kp$type = rep(c("SNV", "indel"), times=ntype)
  }

  dat2kp$vafTumor  = dat2kp$nAltTumor/(dat2kp$nAltTumor + dat2kp$nRefTumor)
  dat2kp$vafNormal = dat2kp$nAltNormal/(dat2kp$nAltNormal + dat2kp$nRefNormal)

  ff1 = sprintf("../data_strelka/strelka_%s.txt", id1)

  write.table(dat2kp, ff1, append = FALSE, quote = FALSE,
            sep = "\t", row.names = TRUE, col.names = TRUE)
}


q(save="no")

