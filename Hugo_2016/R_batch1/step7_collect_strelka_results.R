
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

# ------------------------------------------------------------------------
# collect strelka_results
# ------------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

nn = length(subjects)

for(sub1 in subjects){
  dir1 = paste0("strelka/", sub1, "/results")
  
  f1 = sprintf("%s/passed.somatic.snvs.vcf", dir1)
  f2 = sprintf("%s/passed.somatic.indels.vcf", dir1)

  if(! file.exists(f1)){
    cat("cannot find file", f1, "\n")
  }
  
  if(! file.exists(f2)){
    cat("cannot find file", f2, "\n")
  }
  
  cmd1 = sprintf("cp %s strelka_results/%s_passed.somatic.snvs.vcf", f1, sub1)
  cmd2 = sprintf("cp %s strelka_results/%s_passed.somatic.indels.vcf", f2, sub1)
  
  system(cmd1)
  system(cmd2)
}

q(save="no")

