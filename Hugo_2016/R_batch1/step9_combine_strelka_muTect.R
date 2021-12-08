
library(data.table)

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

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016/data_strelka")

setequal(paste0("strelka_", subjects, ".txt"), list.files(pattern="strelka_"))

strelka = list()

for(sub1 in subjects){
  fnm1 = paste0("strelka_", sub1, ".txt")
  v1   = read.table(fnm1, header=TRUE, sep="\t", as.is=TRUE)
  strelka[[sub1]] = v1
}

# ----------------------------------------------------------------------
# load data of muTect, check the number of called mutations
# here we look at all the mutations instead of those filtered ones
# ----------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016/muTect")

setequal(paste0(subjects, "_hg38_muTect_call_stats.txt"), 
         list.files(pattern="_call_stats.txt"))

muTect = list()
nInt   = rep(0, length(subjects))

for(i in 1:length(subjects)){
  if(i %% 5 == 0){
    cat(i, date(), "\n")
  }
  
  sub1 = subjects[i]
  fnm1 = paste0(sub1, "_hg38_muTect_call_stats.txt")
  v1   = fread(fnm1)
  
  w2kp = which(v1$judgement=="KEEP")
  v1   = v1[w2kp,]
  
  nms  = paste0(v1$contig, ":", v1$position, "_", v1$ref_allele, "/", v1$alt_allele)
  rownames(v1) = nms
  
  muTect[[sub1]] = v1
  
  nInt[i] = length(intersect(rownames(strelka[[sub1]]), rownames(v1)))
}

nMt0 = sapply(strelka, function(x){length(which(x$type=="SNV"))})
summary(nMt0)

nMt1 = sapply(muTect, nrow)
summary(nMt1)

summary(nInt/nMt0)
summary(nInt/nMt1)

# ----------------------------------------------------------------------
# check the count
# ----------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016/")

summary(log10(nMt0 + 1))
summary(log10(nMt1 + 1))

pdf("figures/n_mutations_strealka_vs_muTect.pdf",
      width=12, height=8)

par(mfrow=c(2,3), mar=c(5,4,1,1), bty="n", cex=0.9)

hist(log10(nMt0), xlab="log10(# of mutations by strelka)", main="")
hist(log10(nMt1), xlab="log10(# of mutations by muTect", main="")
hist(log10(nInt), xlab="log10(# of mutations by both methods)", main="")

plot(log10(nMt0 + 1), log10(nMt1 + 1), xlim=c(0,5), ylim=c(0,5), 
     xlab="log10(# of mutations by strelka + 1)",
     ylab="log10(# of mutations by muTect + 1)")
abline(0, 1, col="darkred", lwd=2)

hist(nInt/nMt0, xlab="Prop. strelka mutations called by MuTect", 
     ylab="Freq", main="", breaks=seq(0,1,by=0.1))

hist(nInt/nMt1, xlab="Prop. MuTect mutations called by strelka", 
     ylab="Freq", main="", breaks=seq(0,1,by=0.1))

dev.off()

# ----------------------------------------------------------------------
# take intersection
# ----------------------------------------------------------------------

snvs = NULL
indels = NULL

col1 = rgb(0.8, 0.2, 0.2, 0.5)

filter1 <- function(snvs, rd.n, rd.t, nalt, lod, qss.nt){
  w2kp = (snvs$n_q20_count >= rd.n & snvs$rdNormalFiltered >= rd.n)
  w2kp = w2kp & (snvs$n_ref_count + snvs$n_alt_count >= rd.n)
  w2kp = w2kp & (snvs$nAltNormal + snvs$nRefNormal >= rd.n)
  
  w2kp = (w2kp & snvs$t_q20_count >= rd.t & snvs$rdTumorFiltered >= rd.t)
  w2kp = w2kp & (snvs$t_ref_count + snvs$t_alt_count >= rd.t)
  w2kp = w2kp & (snvs$nAltTumor  + snvs$nRefTumor >= rd.t)
  
  w2kp = (w2kp & snvs$t_alt_count >= nalt & snvs$nAltTumor >= nalt)
  w2kp = (w2kp & snvs$t_lod_fstar >= lod & snvs$QSS_NT >= qss.nt)
  
  w2kp
}

w2kp1 <- function(snvs){
  filter1(snvs, rd.n=20, rd.t=20, nalt=5, lod=10, qss.nt=20)
}


for(id1 in subjects){
  
  v1  = strelka[[id1]]
  v2  = muTect[[id1]]
  
  windel = which(v1$type == "indel")
  if(length(windel) > 0){
    vindel = v1[windel,,drop=FALSE]
    vindel = cbind(rep(id1, nrow(vindel)), vindel)
    names(vindel)[1] = "id"
    indels = rbind(indels, vindel)
  }
  
  if(nrow(v1) == 0 || nrow(v2) == 0){ next }
  
  ww1 = which(rownames(v1) %in% rownames(v2))
  if(length(ww1) == 0){ next }
  
  v1  = v1[ww1,]
  v2  = v2[match(rownames(v1), rownames(v2)),]
  
  v3  = cbind(rep(id1, nrow(v1)), v1, v2)
  names(v3)[1] = "id"
  
  w2kp = w2kp1(v3)

  ff1 = sprintf("figures/_strelka_vs_muTect/%s.png", id1)
  png(ff1, width=9, height=6, res=300, units="in")
  par(mfrow=c(2,3), bty="n", cex=0.8)
  plot(log10(v3$n_q20_count), log10(v3$rdNormalFiltered),  xlab="MuTect",
      ylab="Strelka", main="log10(rd) @ normal", col="grey")
  abline(0, 1, lwd=2, col="darkgrey")
      
  points(log10(v3$n_q20_count)[w2kp], log10(v3$rdNormalFiltered)[w2kp], 
         pch=20, col=col1)
  
  plot(log10(v3$t_q20_count), log10(v3$rdTumorFiltered),  xlab="MuTect",
      ylab="Strelka", main="log10(rd) @ tumor", col="grey")
  abline(0, 1, lwd=2, col="darkgrey")
  
  points(log10(v3$t_q20_count)[w2kp], log10(v3$rdTumorFiltered)[w2kp], 
         pch=20, col=col1)
  
  plot(log10(v3$QSS_NT), log10(v3$t_lod_fstar), col="grey", main=id1, 
      xlab="log10(Strelka QSS_NT)", ylab="log10(MuTect LOD)")
  points(log10(v3$QSS_NT)[w2kp], log10(v3$t_lod_fstar)[w2kp], pch=20, col=col1)
  abline(h = log10(10), col="darkgrey")
  abline(v = log10(20), col="darkgrey")

  normal_f_strelka = v3$nAltNormal/(v3$nAltNormal + v3$nRefNormal)
  plot(normal_f_strelka, v3$normal_f, col="grey", main="VAF @ normal",
      xlab="Strelka", ylab="MuTect")
  abline(0, 1, lwd=2, col="darkgrey")
  
  points(normal_f_strelka[w2kp], v3$normal_f[w2kp], pch=20, col=col1)

  tumor_f_strelka = v3$nAltTumor/(v3$nAltTumor + v3$nRefTumor)
  plot(tumor_f_strelka, v3$tumor_f, col="grey", main="VAF @ tumor",
    xlab="Strelka", ylab="MuTect")
  abline(0, 1, lwd=2, col="darkgrey")
  
  points(tumor_f_strelka[w2kp], v3$tumor_f[w2kp], pch=20, col=col1)
  
  plot(log10(v3$QSS_NT), log10(tumor_f_strelka), col="grey",
      xlab="log10(Strelka QSS_NT)", ylab="log10(tumor VAF)")
  abline(h = -1, col="skyblue")
  abline(v = log10(20), col="darkgrey")

  points(log10(v3$QSS_NT)[w2kp], log10(tumor_f_strelka)[w2kp], pch=20, col=col1)

  dev.off()
  
  snvs = rbind(snvs, v3)
}

dim(snvs)
snvs[1:2,]

dim(indels)
indels[1:2,]

# ----------------------------------------------------------------------
# write out intersections
# ----------------------------------------------------------------------

write.table(snvs, "snvs_by_strelka_and_Mutect.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, 
            col.names = TRUE)

write.table(indels, "indels_by_strelka.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = TRUE, 
            col.names = TRUE)

# ----------------------------------------------------------------------
# compare VAF of strelka and muTect
# ----------------------------------------------------------------------

tumor_f_strelka = snvs$nAltTumor/(snvs$nAltTumor + snvs$nRefTumor)
tumor_f = snvs$tumor_f

summary(tumor_f)
summary(tumor_f_strelka)

table(snvs$judgement, snvs$failure_reasons)

w2kp = w2kp1(snvs)

col1 = rgb(0.8, 0.2, 0.2, 0.5)

png("figures/VAF_strelka_vs_muTect_mutations_pass_by_both.png", width=7, 
    height=3.5, units="in", res=300)
par(mar=c(5,4,3,1), mfrow=c(1,2), bty="n")
plot(tumor_f_strelka, tumor_f, xlab="Strelka", ylab="MuTect",
  main="VAF", col="grey", cex=0.5)
points(tumor_f_strelka[w2kp], tumor_f[w2kp], pch=20, col=col1, cex=0.5)

plot(log10(tumor_f_strelka), log10(tumor_f), xlab="Strelka",
     ylab="MuTect", main="log10(VAF)", col="grey", cex=0.5)
points(log10(tumor_f_strelka)[w2kp], log10(tumor_f)[w2kp], pch=20, 
       col=col1, cex=0.5)
dev.off()

# ----------------------------------------------------------------------
# additional summary of read-depth
# ----------------------------------------------------------------------

table(snvs$n_q20_count >= 20, snvs$rdNormalFiltered >= 20)
table(snvs$t_q20_count >= 20, snvs$rdTumorFiltered >= 20)

table(snvs$n_ref_count + snvs$n_alt_count >= 20)
table(snvs$nAltNormal + snvs$nRefNormal >= 20)

table(snvs$t_ref_count + snvs$t_alt_count >= 20)
table(snvs$nAltTumor  + snvs$nRefTumor >= 20)

table(snvs$t_alt_count >= 5, snvs$nAltTumor >= 5)

summary(snvs$rdNormalFiltered - (snvs$nAltNormal + snvs$nRefNormal))
summary(snvs$rdTumorFiltered  - (snvs$nAltTumor  + snvs$nRefTumor))

summary(snvs$t_q20_count  - (snvs$t_ref_count + snvs$t_alt_count))
summary(snvs$n_q20_count  - (snvs$n_ref_count + snvs$n_alt_count))

summary(snvs$t_q20_count / (snvs$nAltTumor  + snvs$nRefTumor))
summary(snvs$n_q20_count / (snvs$nAltNormal + snvs$nRefNormal))

q(save="no")

