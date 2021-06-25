library(data.table)
#--------------------------------------------------------------------
# 1. read in neoantigen estimation
#--------------------------------------------------------------------

mhci = readRDS("../data/neoAg_netMHCpan4_1.rds")
mhciipan = readRDS("../data/neoAg_netMHCIIpan4_0.rds")

colnames(mhci)
dim(mhci)
head(mhci)

dim(mhciipan)
head(mhciipan)

#--------------------------------------------------------------------
# 2. compare reference and mutated peptides for HLA-I
#--------------------------------------------------------------------

rank.diff.cutoffs = c(seq(0.5, 5, by=0.5), seq(6, 20, by=1), seq(25, 50, by=5))

cts.tbl.mhci = matrix(NA, nrow=length(rank.diff.cutoffs), ncol=4)

for(i in 1:length(rank.diff.cutoffs)){
  rc1 = rank.diff.cutoffs[i]
  
  if(i < length(rank.diff.cutoffs)){
    rc2 = rank.diff.cutoffs[i+1]
    wi = mhci$EL_rank_mut - mhci$EL_rank_ref < -rc1
    wi = wi & mhci$EL_rank_mut - mhci$EL_rank_ref >= -rc2
    cts.tbl.mhci[i,1] = sum(wi)
    cts.tbl.mhci[i,3] = sum(wi & mhci$EL_rank_mut < 2)
    
    wi = mhci$EL_rank_ref - mhci$EL_rank_mut < -rc1
    wi = wi & mhci$EL_rank_ref - mhci$EL_rank_mut >= -rc2
    cts.tbl.mhci[i,2] = sum(wi)
    cts.tbl.mhci[i,4] = sum(wi & mhci$EL_rank_ref < 2)
  }else{
    wi = mhci$EL_rank_mut - mhci$EL_rank_ref < -rc1
    cts.tbl.mhci[i,1] = sum(wi)
    cts.tbl.mhci[i,3] = sum(wi & mhci$EL_rank_mut < 2)
    
    wi = mhci$EL_rank_ref - mhci$EL_rank_mut < -rc1
    cts.tbl.mhci[i,2] = sum(wi)
    cts.tbl.mhci[i,4] = sum(wi & mhci$EL_rank_ref < 2)
  }
}

cts.tbl.mhci

df1.mhci = data.frame(cutoff = rank.diff.cutoffs, cts.tbl.mhci, 
                      1-cts.tbl.mhci[,2]/cts.tbl.mhci[,1], 
                      1-cts.tbl.mhci[,4]/cts.tbl.mhci[,3],
                      stringsAsFactors = FALSE)

names(df1.mhci)[2:3] = c("mut.lower.rank", "ref.lower.rank")
names(df1.mhci)[4:5] = c("mut.lower.rank.abs.le.2", "ref.lower.rank.abs.le.2")
names(df1.mhci)[6:7] = c("TDR", "TDR.abs.le.2")

options(width=150)
dim(df1.mhci)
df1.mhci


#--------------------------------------------------------------------
# 3. compare reference and mutated peptides for HLA-II
#--------------------------------------------------------------------

cts.tbl.mhciipan = matrix(NA, nrow=length(rank.diff.cutoffs), ncol=4)

for(i in 1:length(rank.diff.cutoffs)){
  rc1 = rank.diff.cutoffs[i]
  
  if(i < length(rank.diff.cutoffs)){
    rc2 = rank.diff.cutoffs[i+1]
    wi = mhciipan$EL_rank_mut - mhciipan$EL_rank_ref < -rc1
    wi = wi & mhciipan$EL_rank_mut - mhciipan$EL_rank_ref >= -rc2
    cts.tbl.mhciipan[i,1] = sum(wi)
    cts.tbl.mhciipan[i,3] = sum(wi & mhciipan$EL_rank_mut < 2)
    
    wi = mhciipan$EL_rank_ref - mhciipan$EL_rank_mut < -rc1
    wi = wi & mhciipan$EL_rank_ref - mhciipan$EL_rank_mut >= -rc2
    cts.tbl.mhciipan[i,2] = sum(wi)
    cts.tbl.mhciipan[i,4] = sum(wi & mhciipan$EL_rank_ref < 2)
  }else{
    wi = mhciipan$EL_rank_mut - mhciipan$EL_rank_ref < -rc1
    cts.tbl.mhciipan[i,1] = sum(wi)
    cts.tbl.mhciipan[i,3] = sum(wi & mhciipan$EL_rank_mut < 2)
    
    wi = mhciipan$EL_rank_ref - mhciipan$EL_rank_mut < -rc1
    cts.tbl.mhciipan[i,2] = sum(wi)
    cts.tbl.mhciipan[i,3] = sum(wi & mhciipan$EL_rank_ref < 2)
  }
}

df1.mhciipan = data.frame(cutoff = rank.diff.cutoffs, cts.tbl.mhciipan, 
                          1-cts.tbl.mhciipan[,2]/cts.tbl.mhciipan[,1], 
                          1-cts.tbl.mhciipan[,4]/cts.tbl.mhciipan[,3], 
                          stringsAsFactors = FALSE)

names(df1.mhciipan)[2:3] = c("mut.lower.rank", "ref.lower.rank")
names(df1.mhciipan)[4:5] = c("mut.lower.rank.abs.le.2", "ref.lower.rank.abs.le.2")
names(df1.mhciipan)[6:7] = c("TDR", "TDR.abs.le.2")

dim(df1.mhciipan)
df1.mhciipan

#--------------------------------------------------------------------
# 4. fit the relation between cutoffs and TDR
#--------------------------------------------------------------------

pdf("../figures/rank_diff_vs_rank_diff_cond_rank<2.pdf", width=6, height=3)
par(mfrow=c(1,2), mar=c(5,4,2,1), bty="n")

plot(rank.diff.cutoffs, df1.mhci$TDR, ylim=c(0,1), ylab="TDR", 
     pch=20, col="darkgrey", main="MHC-I")
lines(rank.diff.cutoffs, df1.mhci$TDR, lwd=1, col="darkgrey")
points(rank.diff.cutoffs, df1.mhci$TDR.abs.le.2, pch=20, col="purple")
lines(rank.diff.cutoffs, df1.mhci$TDR.abs.le.2, lwd=1.5, col="purple")
legend("bottomright", legend=c("all", "rank < 2%"), pch=c(20,20), 
       col=c("darkgrey", "purple"), bty="n")


plot(rank.diff.cutoffs, df1.mhciipan$TDR, ylim=c(0,1), ylab="TDR", 
     pch=20, col="darkgrey", main="MHC-IIpan")
lines(rank.diff.cutoffs, df1.mhciipan$TDR, lwd=1, col="darkgrey")
points(rank.diff.cutoffs, df1.mhciipan$TDR.abs.le.2, pch=20, col="purple")
lines(rank.diff.cutoffs, df1.mhciipan$TDR.abs.le.2, lwd=1.5, col="purple")
legend("bottomright", legend=c("all", "rank < 2%"), pch=c(20,20), 
       col=c("darkgrey", "purple"), bty="n")

dev.off()

pdf("../figures/rank_diff_vs_TDR_lowess.pdf", width=6, height=3)
par(mfrow=c(1,2), mar=c(5,4,2,1), bty="n")

plot(rank.diff.cutoffs, df1.mhci$TDR, ylab="TDR", 
     pch=20, col="darkred", main="MHC-I",ylim = c(0,1))
l1 = lowess(rank.diff.cutoffs, df1.mhci$TDR, f=1/3)
lines(l1, lwd=1.5, col="orange")

plot(rank.diff.cutoffs, df1.mhciipan$TDR, ylab="TDR", 
     pch=20, col="darkblue", main="MHC-IIpan",ylim = c(0,1))
l2 = lowess(rank.diff.cutoffs, df1.mhciipan$TDR, f=1/3)
lines(l2, lwd=1.5, col="deepskyblue")
dev.off()

#--------------------------------------------------------------------
# 5. add weights
#--------------------------------------------------------------------

l1 = as.data.frame(l1)
l2 = as.data.frame(l2)

l1
l2

#mhci$rank.diff = mhci$rank - mhci$rank.ref
mhci$rank.diff = mhci$EL_rank_mut - mhci$EL_rank_ref
mhci$weight = rep(0, nrow(mhci))

for(k in 1:nrow(l1)){
  rc1 = l1$x[k]
  
  if(k < nrow(l1)){
    rc2 = l1$x[k + 1]
    w2update = which(mhci$rank.diff < -rc1 & mhci$rank.diff >= -rc2)
    mhci$weight[w2update] = l1$y[k]
  }else{
    w2update = which(mhci$rank.diff < -rc1)
    mhci$weight[w2update] = l1$y[k]
  }
}


#mhciipan$rank.diff = mhciipan$rank - mhciipan$rank.ref
mhciipan$rank.diff = mhciipan$EL_rank_mut - mhciipan$EL_rank_ref
mhciipan$weight = rep(0, nrow(mhciipan))

for(k in 1:nrow(l2)){
  rc1 = l2$x[k]
  
  if(k < nrow(l2)){
    rc2 = l2$x[k + 1]
    w2update = which(mhciipan$rank.diff < -rc1 & mhciipan$rank.diff >= -rc2)
    mhciipan$weight[w2update] = l2$y[k]
  }else{
    w2update = which(mhciipan$rank.diff < -rc1)
    mhciipan$weight[w2update] = l2$y[k]
  }
}

dim(mhci)
head(mhci)

dim(mhciipan)
head(mhciipan)

#--------------------------------------------------------------------
# 6. estimate neoAg burden
#--------------------------------------------------------------------

mhci = as.data.table(mhci)
mhciipan = as.data.table(mhciipan)
mhci
mhciipan

length(unique(mhci$sample))
length(unique(mhciipan$sample))

nb = list()
colnames(mhci)
nb$mhci = mhci[EL_rank_mut < 2, .N, by=sample]
mhci[EL_rank_mut < 2, .N, by=sample]
nb$mhci.lt.ref = mhci[EL_rank_mut < 2 & rank.diff < -0.5, .N, by=sample]
nb$mhci.lt.ref.weighted = mhci[EL_rank_mut < 2 & rank.diff < -0.5, sum(weight), 
                               by=sample]

nb$mhciipan = mhciipan[EL_rank_mut < 2, .N, by=sample]
nb$mhciipan.lt.ref = mhciipan[EL_rank_mut < 2 & rank.diff < -0.5, .N, by=sample]
nb$mhciipan.lt.ref.weighted = mhciipan[EL_rank_mut < 2 & rank.diff < -0.5, sum(weight), 
                                       by=sample]
length(nb)
samples = unique(mhci$sample)
samples
samples_i = unique(mhci$sample)
samples_ii = unique(mhciipan$sample)
nb.mx = matrix(0, nrow=length(samples), ncol=6)
rownames(nb.mx) = samples_i

for(i in 1:3){
  w2update = match(nb[[i]]$sample, samples_i)
  nb.mx[w2update,i] = unlist(nb[[i]][,2])
}

for(i in 4:6){
  w2update = match(nb[[i]]$sample, samples_ii)
  nb.mx[w2update,i] = unlist(nb[[i]][,2])
}



colnames(nb.mx) = names(nb)
dim(nb.mx)
nb.mx[1:6,]

summary(nb.mx)
pairs(log10(nb.mx+1))

plot(log10(nb.mx[,"mhci.lt.ref"] + 1), 
     log10(nb.mx[,"mhci.lt.ref.weighted"] + 1))

# check the weighted summation of mutaiton burdern for a few samples

wts.test1 = wts.test2 = rep(NA, length(samples))
for(i in 1:length(samples)){
  sam1 = samples[i]
  wsam = which(mhci$sample == sam1 & mhci$rank.diff < -0.5 & mhci$EL_rank_mut < 2)
  wts.test1[i] = sum(mhci$weight[wsam])
  
  wsam = which(mhciipan$sample == sam1 & mhciipan$rank.diff < -0.5 & mhciipan$rank < 2)
  wts.test2[i] = sum(mhciipan$weight[wsam])
}

summary(wts.test1 - nb.mx[,"mhci.lt.ref.weighted"])
summary(wts.test2 - nb.mx[,"mhciipan.lt.ref.weighted"])

#--------------------------------------------------------------------
# Step 7. Save datasets 
#--------------------------------------------------------------------

nb.df = data.frame(sample = rownames(nb.mx), nb.mx, stringsAsFactors = FALSE)
dim(nb.df)
nb.df[1:6,]

write.table(nb.df, file = "../output/neoAg_burden.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

sessionInfo()
q(save="no")
