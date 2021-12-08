
#--------------------------------------------------------------------
# Step 7: Import neoAg prediction
#--------------------------------------------------------------------

library(magrittr)
library(Biostrings)
library(ggplot2)
library(plyr)
library(data.table)

# Sample to check 
sample_ID = "Pt10"
#--------------------------------------------------------------------
# 1. Import MHC-I prediction
#--------------------------------------------------------------------

# check the relation between affinity and rank for one sample
dat = fread(paste("for_netmhci/_output_tsv/", sample_ID, "_sm9.tsv", sep = ""))
dim(dat)
dat[1:2,]
table(dat$HLA)

summary(dat$affinity)
summary(dat$rank)

dat$t.affinity = 1 - log(dat$affinity, 50000)
summary(dat$t.affinity)
dat$log10.rank = log10(dat$rank/100)
dim(dat)
dat[1:2,]


pdf(paste("figures/affinity_rank_hist_", sample_ID, "_MHCi.pdf", sep = ""), width=7, height=3)
par(mar=c(5,4,1,1), mfrow=c(1,2))
hist(dat$t.affinity, xlab="1 - log(affinity, 50k)", ylab="freq", main="")
abline(v=1-log(500, 50000), col = "black")
abline(v=1-log(50, 50000), col = "red", lty=2)

hist(dat$log10.rank, xlab="log10(rank)", ylab="freq", main="")
abline(v=log10(0.02), col = "black")
abline(v=log10(0.005), col = "red", lty=2)
dev.off()

table(dat$rank < 0.1)
table(dat$rank < 0.5)
table(dat$rank < 2)

table(dat$affinity < 50)
table(dat$affinity < 500)


g1 = ggplot(dat, aes(x=t.affinity, y=log10.rank, shape=HLA, color=HLA)) + 
  geom_point() + labs(x="1 - log(affinity, 50k)", y = "log10(rank)")
g1 = g1 + theme_classic()

g1 = g1 + geom_hline(yintercept=log10(0.02), color = "black")
g1 = g1 + geom_hline(yintercept=log10(0.005), linetype="dashed", color = "red")

g1 = g1 + geom_vline(xintercept=1-log(500, 50000), color = "black")
g1 = g1 + geom_vline(xintercept=1-log(50, 50000), 
                     linetype="dashed", color = "red")
ggsave(paste("figures/affinity_vs_rank_", sample_ID, "_MHCi.png", sep =""), plot = g1, device = "png", 
       width = 4.5, height = 3, units = "in", dpi=300)

# read in MHC-I binding rank information

mhci.files = list.files("for_netmhci/_output_tsv", full.names=TRUE)
length(mhci.files)
mhci.files[1:5]

sampleID = gsub("for_netmhci/_output_tsv/", "", mhci.files, fixed=TRUE)
sampleID = gsub(".tsv", "", sampleID, fixed=TRUE)
sampleID[1:5]

mhci = list()

for(i in 1:length(mhci.files)){
  fnm = mhci.files[i]
  dat = fread(fnm)
  
  dmin = dat[dat[ , .I[which.min(rank)], by = .(HLA, identity)]$V1]

  mhci[[sampleID[i]]] = dmin
}

#--------------------------------------------------------------------
# 2. Import MHC-II prediction
#--------------------------------------------------------------------

# check the relation between affinity and rank for one sample
dat = fread(paste("for_netmhciipan/_output_tsv/", sample_ID, "_hlaii_sm15.tsv", sep=""))
dim(dat)
dat[1:2,]
table(dat$allele)

summary(dat$affinity)
summary(dat$rank)

dat$t.affinity = 1 - log(dat$affinity, 50000)
summary(dat$t.affinity)
dat$log10.rank = log10(dat$rank/100)
dim(dat)
dat[1:2,]

pdf(paste("figures/affinity_rank_hist_",sample_ID, "_MHCiipan.pdf", sep=""), width=7, height=3)
par(mar=c(5,4,1,1), mfrow=c(1,2))
hist(dat$t.affinity, xlab="1 - log(affinity, 50k)", ylab="freq", main="")
abline(v=1-log(500, 50000), col = "black")
abline(v=1-log(50, 50000), col = "red", lty=2)

hist(dat$log10.rank, xlab="log10(rank)", ylab="freq", main="")
abline(v=log10(0.02), col = "black")
abline(v=log10(0.005), col = "red", lty=2)
dev.off()

table(dat$rank < 0.1)
table(dat$rank < 0.5)
table(dat$rank < 2)

table(dat$affinity < 50)
table(dat$affinity < 500)

g1 = ggplot(dat, aes(x=t.affinity, y=log10.rank, color=allele)) + 
  geom_point() + labs(x="1 - log(affinity, 50k)", y = "log10(rank)")
g1 = g1 + theme_classic()

g1 = g1 + geom_hline(yintercept=log10(0.02), color = "black")
g1 = g1 + geom_hline(yintercept=log10(0.005), linetype="dashed", color = "red")

g1 = g1 + geom_vline(xintercept=1-log(500, 50000), color = "black")
g1 = g1 + geom_vline(xintercept=1-log(50, 50000), 
                     linetype="dashed", color = "red")
ggsave(paste("figures/affinity_vs_rank_", sample_ID, "_MHCiipan.png", sep = ""), plot = g1, device = "png", 
       width = 5.5, height = 3, units = "in", dpi=300)

# read in MHC-IIpan binding rank information

mhciipan.files = list.files("for_netmhciipan/_output_tsv", full.names=TRUE)
length(mhciipan.files)
mhciipan.files[1:5]

sampleID = gsub("for_netmhciipan/_output_tsv/", "", mhciipan.files, fixed=TRUE)
sampleID = gsub(".tsv", "", sampleID, fixed=TRUE)
sampleID[1:5]

mhciipan = list()

for(i in 1:length(mhciipan.files)){
  fnm = mhciipan.files[i]
  dat = fread(fnm)
  
  dmin = dat[dat[ , .I[which.min(rank)], by = .(allele, identity)]$V1]
  
  mhciipan[[sampleID[i]]] = dmin
}

#--------------------------------------------------------------------
# 3. re-organize results for HLA-I
#--------------------------------------------------------------------

which.ref = grep("_sm9r", names(mhci), fixed=TRUE)
nms.m = names(mhci)[-which.ref]
nms.r = names(mhci)[which.ref]

length(nms.m)
length(nms.r)

table(nms.r == paste0(nms.m, "r"))

neoAg.mhci = NULL

for(i in 1:length(nms.m)){
  print(i)
  m1 = mhci[[nms.m[i]]]
  r1 = mhci[[nms.r[i]]]
  print(m1)
  print(r1)
  
  sample = gsub("_sm9", "", nms.m[i])
  
  if(nrow(m1) != nrow(r1)){ stop("nrow do not match\n") }
  if(any(m1$HLA != r1$HLA)){ stop("HLA do not match\n") }
  if(any(m1$identity != r1$identity)){ stop("identity do not match\n") }
  
  names(r1)[-c(2,4)] = paste(names(r1)[-c(2,4)], "ref", sep=".")
  
  m1  = cbind(m1, r1[,-c(2,4)])
  df1 = data.frame(sample=rep(sample, nrow(m1)), m1, stringsAsFactors=FALSE)
  
  neoAg.mhci = rbind(neoAg.mhci, df1)
}

dim(neoAg.mhci)
head(neoAg.mhci)

table(neoAg.mhci$pos == neoAg.mhci$pos.ref)

png("figures/compare_rank_mutation_vs_reference_mhci.png", 
    width=8, height=4, units="in", res=300)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
plot(log10(neoAg.mhci$rank.ref), log10(neoAg.mhci$rank), pch=20, cex=0.5, 
     col="darkgrey", xlab="log10(rank.reference)", ylab="log10(rank.mutation)")
abline(0, 1, lwd=2, col='red')

rank.diff = neoAg.mhci$rank - neoAg.mhci$rank.ref

hist(rank.diff[rank.diff !=0 ], breaks=100, main="", 
     xlab="rank.mutation - rank.reference")
abline(v=0, lwd=2, col='red')
dev.off()

#--------------------------------------------------------------------
# 4. re-organize results for HLA-II
#--------------------------------------------------------------------

which.ref = grep("_sm15r", names(mhciipan), fixed=TRUE)
nms.m = names(mhciipan)[-which.ref]
nms.r = names(mhciipan)[which.ref]

length(nms.m)
length(nms.r)
nms.m[1:5]
nms.r[1:5]

table(nms.r == paste0(nms.m, "r"))

neoAg.mhciipan = NULL

for(i in 1:length(nms.m)){
  print(i)
  print(nms.m[i])
  print(nms.r[i])
  m1 = mhciipan[[nms.m[i]]]
  r1 = mhciipan[[nms.r[i]]]
  
  sample = gsub("_hlaii_sm15", "", nms.m[i])
  
  if(nrow(m1) != nrow(r1)){ stop("nrow do not match\n") }
  if(any(m1$allele != r1$allele)){ stop("HLA do not match\n") }
  if(any(m1$identity != r1$identity)){ stop("identity do not match\n") }
  
  names(r1)[-c(2,4)] = paste(names(r1)[-c(2,4)], "ref", sep=".")
  
  m1  = cbind(m1, r1[,-c(2,4)])
  df1 = data.frame(sample=rep(sample, nrow(m1)), m1, stringsAsFactors=FALSE)
  
  neoAg.mhciipan = rbind(neoAg.mhciipan, df1)
}

dim(neoAg.mhciipan)
head(neoAg.mhciipan)

table(neoAg.mhciipan$seq == neoAg.mhciipan$seq.ref)

png("figures/compare_rank_mutation_vs_reference_mhciipan.png", 
    width=8, height=4, units="in", res=300)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
plot(log10(neoAg.mhciipan$rank.ref), log10(neoAg.mhciipan$rank), pch=20, cex=0.5, 
     col="darkgrey", xlab="log10(rank.reference)", ylab="log10(rank.mutation)")
abline(0, 1, lwd=2, col='red')

rank.diff = neoAg.mhciipan$rank - neoAg.mhciipan$rank.ref
hist(rank.diff[rank.diff !=0 ], breaks=100, main="", 
     xlab="rank.mutation - rank.reference")
abline(v=0, lwd=2, col='red')
dev.off()

table(neoAg.mhci$rank == neoAg.mhci$rank.ref)
table(neoAg.mhciipan$rank == neoAg.mhciipan$rank.ref)

table(neoAg.mhci$rank - neoAg.mhci$rank.ref < 0)
table(neoAg.mhciipan$rank - neoAg.mhciipan$rank.ref < 0)

#--------------------------------------------------------------------
# 5. save the neoAg
#--------------------------------------------------------------------

length(mhci)
mhci[[1]]

length(mhciipan)
mhciipan[[1]]

saveRDS(neoAg.mhci, file="output/mhci_min_rank_per_HLA_mutation.rds")
saveRDS(neoAg.mhciipan, file="output/mhciipan_min_rank_per_HLA_mutation.rds")

sessionInfo()
q(save="no")
