
# ------------------------------------------------------------------------
# read in sample information
# ------------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

srr1 = scan(file="SRR_Acc_List.txt", what=character(0))
srr2 = scan(file="SRR_Acc_List_Vanderbilt.txt", what=character(0))

srr1[1:5]
srr2[1:5]

srr = c(srr1, srr2)

info1 = read.table("SraRunTable.txt", sep="\t", as.is=TRUE, header=TRUE)
dim(info1)
info1[1:2,]

info2 = read.table("SraRunTable_Vanderbilt.txt", sep="\t", as.is=TRUE, 
                   header=TRUE)
dim(info2)
info2[1:2,]

table(info1$tissue)
t1 = table(info1$isolate)
t1[1:5]
table(t1)

info1[which(info1$isolate == names(t1)[t1==3]),]

table(tapply(info1$tissue, info1$isolate, paste, collapse=":"))

table(info2$tissue)
table(info2$disease)

t2 = table(info2$isolate)
t2[1:5]
table(t2)

info2[which(info2$isolate == names(t2)[t2==3]),]

dim(info1)
dim(info2)
length(unique(c(info1$Run, info2$Run)))
length(srr)
setequal(c(info1$Run, info2$Run), srr)

# ------------------------------------------------------------------------
# prepare strelka codes for batch 1
# ------------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/R_batch2")

codes = scan("step6_strelka_ex.sh", what=character(), sep="\n", 
             blank.lines.skip=FALSE)
length(codes)
codes[1:8]

dim(info1)
info1[1:2,]

table(info1$tissue)
ww1 = which(info1$tissue == "melanoma biopsy")
info1$Sample_Name[ww1]
table(info1$Library_Name[ww1] == info1$Sample_Name[ww1])

tN = "peripheral blood mononuclear cells"

subjects = info1$Sample_Name[ww1]
subjects = gsub("-baseline", "", subjects)
shells = paste0("step6/step6_strelka_", subjects, ".sh")

for(k in 1:length(subjects)){
  rowId   = ww1[k]
  sub1    = subjects[k]
  id1     = info1$isolate[rowId]
  tSample = info1$Run[rowId]
  nSample = info1$Run[which(info1$tissue == tN & info1$isolate == id1)]
  
  if(length(nSample) != 1 || length(tSample) != 1){
    stop("unexpected number of sample names\n")
  }
  
  codes[4] = sprintf("sampleID=\"%s\"", sub1)
  codes[5] = sprintf("nSample=\"%s\"", nSample)
  codes[6] = sprintf("tSample=\"%s\"", tSample)
  
  cat(codes, file=shells[k], sep="\n")
  cat("\n", file=shells[k], append=TRUE)
}

logFile = paste0("strelka_", subjects, ".txt")
cmds = "sbatch -n 1 -c 2 --mem=16192 --error="
cmds = paste0(cmds, logFile, " --output=", logFile, " ", shells)

# ------------------------------------------------------------------------
# prepare strelka codes for batch 2
# ------------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/R_batch2")

codes = scan("step6_strelka_ex.sh", what=character(), sep="\n", 
             blank.lines.skip=FALSE)
length(codes)
codes[1:8]

dim(info2)
info2[1:2,]

table(info2$tissue)
ww2 = which(info2$tissue != "PBMC")
info2[ww2, c("Sample_Name", "Library_Name", "tissue")]

tN = "PBMC"

subjects = info2$Library_Name[ww2]
subjects = gsub("-PD1Cell2016-WES", "", subjects)
subjects = gsub("-tumor", "", subjects)
subjects = gsub("-", "_", subjects)
subjects = paste0("Vand_", subjects)

shells = paste0("step6/step6_strelka_", subjects, ".sh")

for(k in 1:length(subjects)){
  rowId   = ww2[k]
  sub1    = subjects[k]
  id1     = info2$isolate[rowId]
  tSample = info2$Run[rowId]
  nSample = info2$Run[which(info2$tissue == tN & info2$isolate == id1)]
  
  if(length(nSample) != 1 || length(tSample) != 1){
    stop("unexpected number of sample names\n")
  }
  
  codes[4] = sprintf("sampleID=\"%s\"", sub1)
  codes[5] = sprintf("nSample=\"%s\"", nSample)
  codes[6] = sprintf("tSample=\"%s\"", tSample)
  
  cat(codes, file=shells[k], sep="\n")
  cat("\n", file=shells[k], append=TRUE)
}

logFile = paste0("strelka_", subjects, ".txt")
cmds2 = "sbatch -n 1 -c 2 --mem=16192 --error="
cmds2 = paste0(cmds2, logFile, " --output=", logFile, " ", shells)

cmds = c(cmds, cmds2)

cat(cmds, file="step6_strelka_batch.sh", sep="\n")

q(save="no")

