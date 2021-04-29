
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
library(ggpointdensity)

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

# next combine _mut and _ref
which.ref = grep("_ref", names(hla_i_long), fixed=TRUE)
nms.m = names(hla_i_long)[-which.ref]
nms.r = names(hla_i_long)[which.ref]

length(nms.m)
length(nms.r)

table(gsub("_ref", "", nms.r) == gsub("_mut", "", nms.m))

neoAg_mhc_i = NULL
cols2align = c("Pos", "ID", "HLA")

for(i in 1:length(nms.m)){
  m1 = hla_i_long[[nms.m[i]]]
  r1 = hla_i_long[[nms.r[i]]]
  
  # remove those with ID like "Un_GL000219v1_8"
  w2rm = grep("^Un_", m1$ID)
  if(length(w2rm) > 0){m1 = m1[-w2rm,]}
  
  w2rm = grep("^Un_", r1$ID)
  if(length(w2rm) > 0){r1 = r1[-w2rm,]}
  
  sample = gsub("_mut", "", nms.m[i])
  
  if(nrow(m1) != nrow(r1)){ stop("nrow do not match\n") }
  
  stopifnot(setequal(apply(m1[,cols2align], 1, paste, collapse=":"), 
                     apply(r1[,cols2align], 1, paste, collapse=":")))

  df1 = merge(m1, r1, by = cols2align, 
              suffixes = c("_mut", "_ref"))
  stopifnot(nrow(df1) == nrow(m1))
  
  df1$sample = rep(sample, nrow(df1))
  dim(df1)
  df1[1:5,]
  
  neoAg_mhc_i = rbind(neoAg_mhc_i, df1)
}

dim(neoAg_mhc_i)
head(neoAg_mhc_i)

rank.diff = neoAg_mhc_i$EL_rank_mut - neoAg_mhc_i$EL_rank_ref
table(rank.diff == 0)
table(abs(rank.diff) > 1)
table(abs(rank.diff) > 5)
table(abs(rank.diff) > 10)

table(rank.diff > 5 & rank.diff <= 10)
table(rank.diff < -5 & rank.diff >= -10)

table(rank.diff > 10)
table(rank.diff < -10)

pdf("../figures/compare_rank_mutation_vs_reference_mhc_i_hist1.pdf", 
    width=3, height=3)
par(mar=c(5,4,1,1))
hist(rank.diff[rank.diff !=0 ], breaks=100, main="", 
     xlab="rank.mutation - rank.reference", xlim=c(-20,20))
abline(v=0, lwd=1, col='red')
dev.off()

pdf("../figures/compare_rank_mutation_vs_reference_mhc_i_hist2.pdf", 
    width=5, height=3)
par(mar=c(5,4,1,1))
hist(rank.diff[abs(rank.diff) > 10 ], breaks=100, main="", 
     xlab="rank.mutation - rank.reference")
abline(v=0, lwd=1, col='red')

dev.off()

neoAg_mhc_i_sub = neoAg_mhc_i[which(neoAg_mhc_i$EL_rank_mut < 2),]
dim(neoAg_mhc_i_sub)
neoAg_mhc_i_sub[1:2,]

min(neoAg_mhc_i_sub$EL_rank_mut)
min(neoAg_mhc_i_sub$EL_rank_ref)


gs1 = ggplot(neoAg_mhc_i_sub, aes(x = log10(EL_rank_ref), 
                                 y = log10(EL_rank_mut))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x = "log10(rank_reference)", y = "log10(rank_mutation)") + 
  geom_abline(intercept = 0, slope = 1)

png("../figures/compare_rank_mutation_vs_reference_mhc_i.png", 
    width=6, height=4, units="in", res=300)
gs1
dev.off()

#--------------------------------------------------------------------
# 4. re-organize results for HLA-II
#--------------------------------------------------------------------

# first transform to the long format
hla_ii_long = list()

for(k in names(hla_ii)){
  dat = as.data.frame(hla_ii[[k]])
  
  dat_l1 = pivot_longer(dat[,c(1:3,seq(5,ncol(dat)-2,by=2))], 
                        values_to="EL_score", cols=ends_with("_Score"))
  
  dat_l2 = pivot_longer(dat[,c(1:3,seq(6,ncol(dat)-2,by=2))], 
                        values_to="EL_rank", cols=ends_with("_Rank"))
  
  dat_l1$HLA = gsub("_Score", "", dat_l1$name)
  dat_l2$HLA = gsub("_Rank",  "", dat_l2$name)
  
  dim(dat_l1)
  dat_l1[1:2,]
  
  dim(dat_l2)
  dat_l2[1:2,]
  
  stopifnot(all(dat_l1$Pos  == dat_l2$Pos))
  stopifnot(all(dat_l1$ID   == dat_l2$ID))
  stopifnot(all(dat_l1$HLA  == dat_l2$HLA))
  stopifnot(all(dat_l1$Peptide  == dat_l2$Peptide))
  
  mhc_ii = merge(dat_l1[,-4], dat_l2[,-4])
  hla_ii_long[[k]] = mhc_ii
}

length(hla_ii_long)
sapply(hla_ii_long[1:5], dim)

# next combine _mut and _ref
which.ref = grep("_ref", names(hla_ii_long), fixed=TRUE)
nms.m = names(hla_ii_long)[-which.ref]
nms.r = names(hla_ii_long)[which.ref]

length(nms.m)
length(nms.r)

table(gsub("_ref", "", nms.r) == gsub("_mut", "", nms.m))

neoAg_mhc_ii = NULL
cols2align = c("Pos", "ID", "HLA")

for(i in 1:length(nms.m)){
  m1 = hla_ii_long[[nms.m[i]]]
  r1 = hla_ii_long[[nms.r[i]]]
  
  # remove those with ID like "Un_GL000219v1_8"
  w2rm = grep("^Un_", m1$ID)
  if(length(w2rm) > 0){m1 = m1[-w2rm,]}
  
  w2rm = grep("^Un_", r1$ID)
  if(length(w2rm) > 0){r1 = r1[-w2rm,]}
  
  sample = gsub("_mut", "", nms.m[i])
  
  if(nrow(m1) != nrow(r1)){ stop("nrow do not match\n") }
  
  stopifnot(setequal(apply(m1[,cols2align], 1, paste, collapse=":"), 
                     apply(r1[,cols2align], 1, paste, collapse=":")))
  
  df1 = merge(m1, r1, by = cols2align, 
              suffixes = c("_mut", "_ref"))
  stopifnot(nrow(df1) == nrow(m1))
  
  df1$sample = rep(sample, nrow(df1))
  dim(df1)
  df1[1:5,]
  
  neoAg_mhc_ii = rbind(neoAg_mhc_ii, df1)
}

dim(neoAg_mhc_ii)
head(neoAg_mhc_ii)

rank.diff = neoAg_mhc_ii$EL_rank_mut - neoAg_mhc_ii$EL_rank_ref

table(rank.diff == 0)
table(abs(rank.diff) > 1)
table(abs(rank.diff) > 5)
table(abs(rank.diff) > 10)

table(rank.diff > 0)
table(rank.diff > 5 & rank.diff <= 10)
table(rank.diff < -5 & rank.diff >= -10)

table(rank.diff > 10)
table(rank.diff < -10)

pdf("../figures/compare_rank_mutation_vs_reference_mhc_ii_hist1.pdf", 
    width=3, height=3)
par(mar=c(5,4,1,1))
hist(rank.diff[rank.diff !=0 ], breaks=100, main="", 
     xlab="rank.mutation - rank.reference", xlim=c(-20,20))
abline(v=0, lwd=1, col='red')
dev.off()

pdf("../figures/compare_rank_mutation_vs_reference_mhc_ii_hist2.pdf", 
    width=5, height=3)
par(mar=c(5,4,1,1))
hist(rank.diff[abs(rank.diff) > 10 ], breaks=100, main="", 
     xlab="rank.mutation - rank.reference")
abline(v=0, lwd=1, col='red')

dev.off()

neoAg_mhc_ii_sub = neoAg_mhc_ii[which(neoAg_mhc_ii$EL_rank_mut < 2),]
dim(neoAg_mhc_ii_sub)
neoAg_mhc_ii_sub[1:2,]

min(neoAg_mhc_ii_sub$EL_rank_mut)
min(neoAg_mhc_ii_sub$EL_rank_mut[neoAg_mhc_ii_sub$EL_rank_mut > 0])
neoAg_mhc_ii_sub$EL_rank_mut[which(neoAg_mhc_ii_sub$EL_rank_mut < 1e-2)] = 5e-3

min(neoAg_mhc_ii_sub$EL_rank_ref)
min(neoAg_mhc_ii_sub$EL_rank_ref[neoAg_mhc_ii_sub$EL_rank_ref > 0])
neoAg_mhc_ii_sub$EL_rank_ref[which(neoAg_mhc_ii_sub$EL_rank_ref < 1e-2)] = 5e-3

gs1 = ggplot(neoAg_mhc_ii_sub, aes(x = log10(EL_rank_ref), 
                                 y = log10(EL_rank_mut))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x = "log10(rank_reference)", y = "log10(rank_mutation)") + 
  geom_abline(intercept = 0, slope = 1)

png("../figures/compare_rank_mutation_vs_reference_mhc_ii.png", 
    width=6, height=4, units="in", res=300)
gs1
dev.off()

#--------------------------------------------------------------------
# 5. save the neoAg
#--------------------------------------------------------------------

saveRDS(neoAg_mhc_i,  file="../data/neoAg_netMHCpan4_1.rds")
saveRDS(neoAg_mhc_ii, file="../data/neoAg_netMHCIIpan4_0.rds")

sessionInfo()
q(save="no")
