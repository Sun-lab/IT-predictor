
library(data.table)
library(ggplot2)
library(ggpointdensity)
theme_set(theme_bw())

#--------------------------------------------------------------------
# 1. read in neoantigen estimation
#--------------------------------------------------------------------

mhci  = readRDS("../data/neoAg_netMHCpan4_1.rds")
mhcii = readRDS("../data/neoAg_netMHCIIpan4_0.rds")

colnames(mhci)
dim(mhci)
head(mhci)

dim(mhcii)
head(mhcii)


t1 = table(mhci$ID)
table(t1)

t1[t1 > 500]

mhci_check = mhci[which(mhci$ID == "7_140753336_A_T"),]
dim(mhci_check)
mhci_check[1:5,]

table(mhci_check$sample)
tc = table(apply(mhci_check[,c("Pos", "ID", "sample")], 1, 
                 paste, collapse=":"))
table(tc)

t2 = table(mhcii$ID)
table(t2)
t2[t2 > 1000]

mhcii_check = mhcii[which(mhcii$ID == "7:140753336:A:T"),]
dim(mhcii_check)
mhcii_check[1:5,]

table(mhcii_check$sample)
tc = table(apply(mhcii_check[,c("Pos", "ID", "sample")], 1, 
                 paste, collapse=":"))
table(tc)

#--------------------------------------------------------------------
# 1b. compare score vs. rank
#--------------------------------------------------------------------

summary(mhci$EL_rank_mut)
summary(mhci$EL_score_mut)

# a random sample of 10,000 peptides
set.seed(2021)
mhci_sam = mhci[sample(nrow(mhci), 10000),]
dim(mhci_sam)

g1 = ggplot(mhci_sam, aes(x=EL_score_mut, y=log10(EL_rank_mut))) + 
  geom_pointdensity() + scale_color_viridis_c() + 
  geom_hline(yintercept=log10(2), col="grey") + ggtitle("MHC-I")

png("../figures/step8_mhc_i_EL_score_mut_vs_EL_rank_mut.png", 
    width=5, height=4, units="in", res=300)
g1
dev.off()

summary(mhcii$EL_rank_mut)
summary(mhcii$EL_score_mut)
min0 = min(mhcii$EL_rank_mut[mhcii$EL_rank_mut > 0])
min0
mhcii$EL_rank_mut[which(mhcii$EL_rank_mut == 0)] = 0.5*min0

# a random sample of 10,000 peptides
set.seed(2021)
mhcii_sam = mhcii[sample(nrow(mhcii), 10000),]
dim(mhcii_sam)

g2 = ggplot(mhcii_sam, aes(x=EL_score_mut, y=log10(EL_rank_mut))) + 
  geom_pointdensity() + scale_color_viridis_c() + 
  geom_hline(yintercept=log10(2), col="grey") + ggtitle("MHC-II")

png("../figures/step8_mhc_ii_EL_score_mut_vs_EL_rank_mut.png", 
    width=5, height=4, units="in", res=300)
g2
dev.off()

#--------------------------------------------------------------------
# 2. compare reference and mutated peptides for HLA-I
#--------------------------------------------------------------------

score.diff.cutoffs = seq(0.05, 0.8, by=0.05)

cts.tbl.mhci = matrix(NA, nrow=length(score.diff.cutoffs), ncol=4)

for(i in 1:length(score.diff.cutoffs)){
  rc1 = score.diff.cutoffs[i]
  
  if(i < length(score.diff.cutoffs)){
    rc2 = score.diff.cutoffs[i+1]
    wi = mhci$EL_score_mut - mhci$EL_score_ref > rc1
    wi = wi & mhci$EL_score_mut - mhci$EL_score_ref <= rc2
    cts.tbl.mhci[i,1] = sum(wi)
    cts.tbl.mhci[i,3] = sum(wi & mhci$EL_score_mut > 0.2)
    
    wi = mhci$EL_score_ref - mhci$EL_score_mut > rc1
    wi = wi & mhci$EL_score_ref - mhci$EL_score_mut <= rc2
    cts.tbl.mhci[i,2] = sum(wi)
    cts.tbl.mhci[i,4] = sum(wi & mhci$EL_score_ref > 0.2)
  }else{
    wi = mhci$EL_score_mut - mhci$EL_score_ref > rc1
    cts.tbl.mhci[i,1] = sum(wi)
    cts.tbl.mhci[i,3] = sum(wi & mhci$EL_score_mut > 0.2)
    
    wi = mhci$EL_score_ref - mhci$EL_score_mut > rc1
    cts.tbl.mhci[i,2] = sum(wi)
    cts.tbl.mhci[i,4] = sum(wi & mhci$EL_score_ref > 0.2)
  }
}

cts.tbl.mhci

df1.mhci = data.frame(cutoff = score.diff.cutoffs, cts.tbl.mhci, 
                      round(1-cts.tbl.mhci[,2]/cts.tbl.mhci[,1], 3), 
                      round(1-cts.tbl.mhci[,4]/cts.tbl.mhci[,3], 3), 
                      stringsAsFactors = FALSE)

names(df1.mhci)[2:3] = c("mut_H", "ref_H")
names(df1.mhci)[4:5] = c("mut_H_mut_ge_0.2", "ref_H_ref_ge_0.2")
names(df1.mhci)[6:7] = c("TDR", "TDR_ge_0.2")

dim(df1.mhci)
df1.mhci

#--------------------------------------------------------------------
# 3. compare reference and mutated peptides for HLA-II
#--------------------------------------------------------------------

cts.tbl.mhcii = matrix(NA, nrow=length(score.diff.cutoffs), ncol=4)

for(i in 1:length(score.diff.cutoffs)){
  rc1 = score.diff.cutoffs[i]
  
  if(i < length(score.diff.cutoffs)){
    rc2 = score.diff.cutoffs[i+1]
    wi = mhcii$EL_score_mut - mhcii$EL_score_ref > rc1
    wi = wi & mhcii$EL_score_mut - mhcii$EL_score_ref <= rc2
    cts.tbl.mhcii[i,1] = sum(wi)
    cts.tbl.mhcii[i,3] = sum(wi & mhcii$EL_score_mut > 0.2)
    
    wi = mhcii$EL_score_ref - mhcii$EL_score_mut > rc1
    wi = wi & mhcii$EL_score_ref - mhcii$EL_score_mut <= rc2
    cts.tbl.mhcii[i,2] = sum(wi)
    cts.tbl.mhcii[i,4] = sum(wi & mhcii$EL_score_ref > 0.2)
  }else{
    wi = mhcii$EL_score_mut - mhcii$EL_score_ref > rc1
    cts.tbl.mhcii[i,1] = sum(wi)
    cts.tbl.mhcii[i,3] = sum(wi & mhcii$EL_score_mut > 0.2)
    
    wi = mhcii$EL_score_ref - mhcii$EL_score_mut > rc1
    cts.tbl.mhcii[i,2] = sum(wi)
    cts.tbl.mhcii[i,3] = sum(wi & mhcii$EL_score_ref > 0.2)
  }
}

df1.mhcii = data.frame(cutoff = score.diff.cutoffs, cts.tbl.mhcii, 
                       round(1-cts.tbl.mhcii[,2]/cts.tbl.mhcii[,1], 3),
                       round(1-cts.tbl.mhcii[,4]/cts.tbl.mhcii[,3], 3), 
                       stringsAsFactors = FALSE)

names(df1.mhcii)[2:3] = c("mut_H", "ref_H")
names(df1.mhcii)[4:5] = c("mut_H_mut_ge_0.2", "ref_H_ref_ge_0.2")
names(df1.mhcii)[6:7] = c("TDR", "TDR_ge_0.2")

dim(df1.mhcii)
df1.mhcii

#--------------------------------------------------------------------
# 4. estimate neoAg burden
#--------------------------------------------------------------------

mhci  = as.data.table(mhci)
mhcii = as.data.table(mhcii)
mhci
mhcii

length(unique(mhci$sample))
length(unique(mhcii$sample))

# take maximum across position in the peptide and HLA alleles
mhci_max = mhci[,.(EL_score_mut_max = max(EL_score_mut)), by = .(ID, sample)]
mhci_max

mhcii_max = mhcii[,.(EL_score_mut_max = max(EL_score_mut)), by = .(ID, sample)]
mhcii_max

dim(mhci)
dim(mhcii)
dim(mhci_max)
dim(mhcii_max)

summary(mhci$EL_score_mut)
summary(mhcii$EL_score_mut)

summary(mhci_max$EL_score_mut)
summary(mhcii_max$EL_score_mut)

quantile(mhci$EL_score_mut, c(0.95, 0.98, 0.99))
quantile(mhcii$EL_score_mut, c(0.95, 0.98, 0.99))

quantile(mhci_max$EL_score_mut, c(0.95, 0.98, 0.99))
quantile(mhcii_max$EL_score_mut, c(0.95, 0.98, 0.99))

table(mhci$EL_score_mut > 0.187)/nrow(mhci)
table(mhcii$EL_score_mut > 0.40)/nrow(mhcii)

table(mhci_max$EL_score_mut_max > 0.9154)/nrow(mhci_max)
table(mhcii_max$EL_score_mut_max > 0.9524)/nrow(mhcii_max)

nb = list()
colnames(mhci)

## count one somatic mutation multiple times if multple position of it
## is bound to multiple HLA allele
## or count one somatic mutation at most once 

nb$mhci  = mhci[EL_score_mut > 0.187, .N, by=sample]
nb$mhci_max  = mhci_max[EL_score_mut_max > 0.9154, .N, by=sample]

nb$mhcii = mhcii[EL_score_mut > 0.40, .N, by=sample]
nb$mhcii_max = mhcii_max[EL_score_mut_max > 0.9524, .N, by=sample]

length(nb)
lapply(nb, dim)
lapply(nb, head)

samples_i  = unique(mhci$sample)
samples_ii = unique(mhcii$sample)

length(samples_i)
length(samples_ii)

stopifnot(all(gsub("_hlai",  "", samples_i) == 
                gsub("_hlaii", "", samples_ii)))

samples = gsub("_hlai",  "", samples_i)

nb.mx = matrix(0, nrow=length(samples), ncol=4)
rownames(nb.mx) = samples

for(i in 1:2){
  w2update = match(nb[[i]]$sample, samples_i)
  stopifnot(! any(is.na(w2update)))
  nb.mx[w2update,i] = unlist(nb[[i]][,2])
}

for(i in 3:4){
  w2update = match(nb[[i]]$sample, samples_ii)
  stopifnot(! any(is.na(w2update)))
  nb.mx[w2update,i] = unlist(nb[[i]][,2])
}


colnames(nb.mx) = c("mhci_all", "mhci_max", "mhcii_all", "mhcii_max")
dim(nb.mx)
nb.mx[1:6,]

summary(nb.mx)

pdf(file="../figures/step8_compare_nb.pdf", width=8, height=8)
pairs(log10(nb.mx+1))
dev.off()

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
