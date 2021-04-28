
#--------------------------------------------------------------------
# Step 7: Process neoAg prediction
#--------------------------------------------------------------------

library(magrittr)
library(Biostrings)
library(ggplot2)
library(plyr)
library(data.table)
library(tidyr)
library(ggpubr)
library(scales)

theme_set(theme_classic())

#--------------------------------------------------------------------
# 1. Import MHC-I prediction
#--------------------------------------------------------------------

# check the relation between score and rank for one sample

hla_i = readRDS("../data/netMHCpan4_1_results.rds")
length(hla_i)

dat = as.data.frame(hla_i[["Pt58_pre_hlai_mut"]])
dim(dat)
dat[1:2,]

dat_l1 = pivot_longer(dat[,c(1:3,seq(4,ncol(dat)-2,by=2))], 
                      values_to="EL_score", cols=ends_with("_EL-score"))
dat_l2 = pivot_longer(dat[,c(1:3,seq(5,ncol(dat)-2,by=2))], 
                      values_to="EL_rank", cols=ends_with("_EL_Rank"))

dim(dat_l1)
dat_l1[1:2,]

dim(dat_l2)
dat_l2[1:2,]

dat_l1$HLA = gsub("_EL-score", "", dat_l1$name)
dat_l2$HLA = gsub("_EL_Rank",  "", dat_l2$name)

table(dat_l1$Pos  == dat_l2$Pos)
table(dat_l1$ID   == dat_l2$ID)
table(dat_l1$name == dat_l2$name)
table(dat_l1$Peptide == dat_l2$Peptide)

mhc_i = merge(dat_l1[,-4], dat_l2[,-4])
dim(mhc_i)
mhc_i[1:2,]

table(mhc_i$EL_score == 0)
table(mhc_i$EL_score < 0.001)
table(mhc_i$EL_score < 0.002)
table(mhc_i$EL_score < 0.005)

summary(mhc_i$EL_score)
summary(mhc_i$EL_rank)
summary(mhc_i$EL_rank[which(mhc_i$EL_score <= 0.001)])
summary(mhc_i$EL_rank[which(mhc_i$EL_score > 0.001)])

mhc_i = mhc_i[which(mhc_i$EL_score > 0.001),]
dim(mhc_i)

g1 = ggplot(mhc_i, aes(x=EL_score)) + 
  geom_histogram(color="darkblue", fill="lightblue")

g2 = ggplot(mhc_i, aes(x=EL_rank)) + 
  geom_histogram(color="darkblue", fill="lightblue")

g3 = ggplot(mhc_i, aes(x=EL_score, y=EL_rank, shape=HLA, color=HLA)) + 
  geom_point() + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

gg1 = ggarrange(
  ggarrange(g1, g2, ncol = 2, labels = c("A", "B")), g3, 
  nrow = 2, labels = c("", "C")
) 

pdf("../figures/netMHCpan4.1_score_vs_rank_Pt58_pre.pdf", width=6, height=6)
gg1
dev.off()

#--------------------------------------------------------------------
# 2. Import MHC-II prediction
#--------------------------------------------------------------------

hla_ii = readRDS("../data/netMHCIIpan4_0_results.rds")
length(hla_ii)

dat = as.data.frame(hla_ii[["Pt58_pre_hlaii_mut"]])
dim(dat)
dat[1:2,]

table(dat$Target)

dat_l1 = pivot_longer(dat[,c(1:3,seq(5,ncol(dat)-2,by=2))], 
                      values_to="EL_score", cols=ends_with("_Score"))
dat_l2 = pivot_longer(dat[,c(1:3,seq(6,ncol(dat)-2,by=2))], 
                      values_to="EL_rank", cols=ends_with("_Rank"))

dim(dat_l1)
dat_l1[1:2,]

dim(dat_l2)
dat_l2[1:2,]

dat_l1$HLA = gsub("_Score", "", dat_l1$name)
dat_l2$HLA = gsub("_Rank",  "", dat_l2$name)

table(dat_l1$Peptide == dat_l2$Peptide)
table(dat_l1$Pos == dat_l2$Pos)
table(dat_l1$ID  == dat_l2$ID)
table(dat_l1$HLA == dat_l2$HLA)

mhc_ii = merge(dat_l1[,-4], dat_l2[,-4])
dim(mhc_ii)
mhc_ii[1:2,]

summary(mhc_ii$EL_score)
summary(mhc_ii$EL_rank)

summary(mhc_ii$EL_rank[which(mhc_ii$EL_rank > 0)])
table(mhc_ii$EL_rank == 0)
mhc_ii$EL_rank[which(mhc_ii$EL_rank == 0)] = 0.001


table(mhc_ii$EL_score == 0)
table(mhc_ii$EL_score < 0.001)
table(mhc_ii$EL_score < 0.002)
table(mhc_ii$EL_score < 0.005)
summary(mhc_ii$EL_rank[which(mhc_ii$EL_score <= 0.001)])

mhc_ii = mhc_ii[which(mhc_ii$EL_score > 0.001),]
dim(mhc_ii)
mhc_ii[1:2,]

g1 = ggplot(mhc_ii, aes(x=EL_score)) +
  geom_histogram(color="darkblue", fill="lightblue")

g2 = ggplot(mhc_ii, aes(x=EL_rank)) +
  geom_histogram(color="darkblue", fill="lightblue")

g3 = ggplot(mhc_ii, aes(x=EL_score, y=EL_rank, color=HLA)) + 
  geom_point() + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

gg1 = ggarrange(
  ggarrange(g1, g2, ncol = 2, labels = c("A", "B")), g3, 
  nrow = 2, labels = c("", "C")
) 

pdf("../figures/netMHCIIpan4.0_score_vs_rank_Pt58_pre.pdf", width=6, height=6)
gg1
dev.off()

#--------------------------------------------------------------------
# 3. re-organize results for HLA-I
#--------------------------------------------------------------------

# first transform to the long format
hla_i_long = list()

for(k in names(hla_i)){
  dat = as.data.frame(hla_i[[k]])

  dat_l1 = pivot_longer(dat[,c(1:3,seq(4,ncol(dat)-2,by=2))], 
                        values_to="EL_score", cols=ends_with("_EL-score"))
  dat_l2 = pivot_longer(dat[,c(1:3,seq(5,ncol(dat)-2,by=2))], 
                        values_to="EL_rank", cols=ends_with("_EL_Rank"))
  
  dat_l1$HLA = gsub("_EL-score", "", dat_l1$name)
  dat_l2$HLA = gsub("_EL_Rank",  "", dat_l2$name)
  
  dim(dat_l1)
  dat_l1[1:2,]
  
  dim(dat_l2)
  dat_l2[1:2,]
  
  stopifnot(all(dat_l1$Pos  == dat_l2$Pos))
  stopifnot(all(dat_l1$ID   == dat_l2$ID))
  stopifnot(all(dat_l1$HLA  == dat_l2$HLA))
  stopifnot(all(dat_l1$Peptide  == dat_l2$Peptide))

  mhc_i = merge(dat_l1[,-4], dat_l2[,-4])
  hla_i_long[[k]] = mhc_i
}

length(hla_i_long)
sapply(hla_i_long[1:5], dim)

names(hla_i_long)[1:5]

# next combine _mut and _ref
which.ref = grep("_ref", names(hla_i_long), fixed=TRUE)
nms.m = names(hla_i_long)[-which.ref]
nms.r = names(hla_i_long)[which.ref]

length(nms.m)
length(nms.r)

table(gsub("_ref", "", nms.r) == gsub("_mut", "", nms.m))

neoAg.mhci = NULL

for(i in 1:length(nms.m)){
  m1 = hla_i_long[[nms.m[i]]]
  r1 = hla_i_long[[nms.r[i]]]
  
  sample = gsub("_mut", "", nms.m[i])
  
  if(nrow(m1) != nrow(r1)){ stop("nrow do not match\n") }
  stopifnot(all(m1$Pos  == r1$Pos))
  stopifnot(all(m1$ID   == r1$ID))
  stopifnot(all(m1$HLA  == r1$HLA))
  stopifnot(all(m1$Peptide  == r1$Peptide))
  
  df1 = merge(m1, r1, by = c("Pos", "ID", "HLA"), 
              suffixes = c("_mut", "_ref"))
  
  dim(df1)
  df1[1:5,]
  
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
