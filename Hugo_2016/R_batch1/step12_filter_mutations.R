
# ----------------------------------------------------------------------
# read in mutation calls
# ----------------------------------------------------------------------

setwd("~/research/Immuno/data/Hugo_2016")

muts = read.table("snvs_indels_by_strelka_and_Mutect_with_anno.txt", 
                  header=TRUE, sep="\t", quote="", as.is=TRUE)
dim(muts)
muts[1:2,]

tb.id = table(muts$id)
tb.id.snv   = table(muts$id[which(muts$type=="SNV")])
tb.id.indel = table(muts$id[which(muts$type=="indel")])

sort(tb.id, decreasing=TRUE)[1:10]
sort(tb.id.snv, decreasing=TRUE)[1:10]
sort(tb.id.indel, decreasing=TRUE)[1:10]

# ----------------------------------------------------------------------
# filter by depth
# ----------------------------------------------------------------------

table(muts$nAltTumor + muts$nRefTumor >= 20)
table(muts$nAltNormal + muts$nRefNormal >= 20)
table(muts$vafTumor >= 0.05, muts$nAltTumor >= 5)

## note that indels do not have muTect information
table(muts$t_ref_count + muts$t_alt_count >= 20)
table(muts$n_ref_count + muts$n_alt_count >= 20)
table(muts$tumor_f >= 0.05, muts$t_alt_count >= 5)

w2kp = which(muts$nAltTumor + muts$nRefTumor >= 20)
length(w2kp)

w2kp = setdiff(w2kp, which(muts$nAltNormal + muts$nRefNormal < 20))
length(w2kp)

w2kp = setdiff(w2kp, which(muts$t_ref_count + muts$t_alt_count < 20))
length(w2kp)

w2kp = setdiff(w2kp, which(muts$n_ref_count + muts$n_alt_count < 20))
length(w2kp)

w2kp = setdiff(w2kp, which(muts$nAltTumor < 5))
length(w2kp)

w2kp = setdiff(w2kp, which(muts$t_alt_count < 5))
length(w2kp)

muts = muts[w2kp,]
dim(muts)

# ----------------------------------------------------------------------
# filter by VAF
# ----------------------------------------------------------------------

summary(muts$vafTumor)
summary(muts$tumor_f)

table(muts$vafTumor >= 0.05, muts$tumor_f >= 0.05, useNA = "ifany")

muts = muts[-which(muts$vafTumor < 0.05 | muts$tumor_f < 0.05),]
dim(muts)

# ----------------------------------------------------------------------
# filter by ExAC VAF
# ----------------------------------------------------------------------

table(muts$ExAC_nontcga_ALL > 1e-4, useNA="ifany")
table(muts$ExAC_nontcga_ALL > 1e-3, useNA="ifany")

muts = muts[-which(muts$ExAC_nontcga_ALL > 1e-4),]
dim(muts)

table(muts$Func.ensGene, useNA="ifany")
table(muts$ExonicFunc.ensGene, useNA="ifany")

nonSynonymous = rep(FALSE, nrow(muts))
efuns = c("nonsynonymous SNV", "stopgain", "stoploss")
efuns = c(efuns, "frameshift substitution")
efuns = c(efuns, "nonframeshift substitution")
funs = c("exonic;splicing", "ncRNA_exonic;splicing", "splicing")
nonSynonymous[which(muts$ExonicFunc.ensGene %in% efuns)] = TRUE
nonSynonymous[which(muts$Func.ensGene %in% funs)] = TRUE

table(nonSynonymous)
table(nonSynonymous, muts$Func.ensGene)
table(nonSynonymous, muts$ExonicFunc.ensGene)

muts$nonSynonymous = nonSynonymous

# ----------------------------------------------------------------------
# write out filtered mutation data
# ----------------------------------------------------------------------

write.table(muts, file = "snvs_indels_by_strelka_and_Mutect_with_anno_filtered.txt", 
            append = FALSE, quote = FALSE, sep = "\t",eol = "\n", 
            row.names = FALSE, col.names = TRUE)

# ----------------------------------------------------------------------
# summarize mutation by sample
# ----------------------------------------------------------------------

mut1 = muts[which(muts$nonSynonymous),]
dim(mut1)

tb.id = sort(table(mut1$id), decreasing=TRUE)
tb.id.snv  = table(mut1$id[which(mut1$type=="SNV")])
tb.id.indel = table(mut1$id[which(mut1$type=="indel")])

tb.id[1:10]
sort(tb.id.snv, decreasing=TRUE)[1:10]
sort(tb.id.indel, decreasing=TRUE)[1:10]

tb.id.snv   = tb.id.snv[match(names(tb.id), names(tb.id.snv))]
tb.id.indel = tb.id.indel[match(names(tb.id), names(tb.id.indel))]
tb.id.indel[which(is.na(tb.id.indel))] = 0

pdf("figures/n_non_synonymous_muts.pdf", width=3, height=3)
par(mar=c(5,4,1,1), bty="n")
plot(log10(as.numeric(tb.id.snv)), as.numeric(tb.id.indel), 
     xlab="log10(# of SNVs)", ylab="# of indels", 
     pch=20, col=rgb(0.8, 0.2, 0.2, 0.5))
dev.off()

# ----------------------------------------------------------------------
# summarize mutation by gene
# ----------------------------------------------------------------------

gene.id = unlist(tapply(mut1$Gene.ensGene, mut1$id, unique))
length(gene.id)

tgene = sort(table(gene.id), decreasing=TRUE)
tgene[1:50]
table(tgene > 2)
table(tgene > 3)
table(tgene > 4)

gene2kp = names(tgene)[which(tgene > 2)]
length(gene2kp)
gene2kp[1:2]

dat1 = matrix(NA, nrow=length(tb.id), ncol=length(gene2kp))
rownames(dat1) = names(tb.id)
colnames(dat1) = gene2kp
dim(dat1)
dat1[1:3,1:3]

for(i in 1:nrow(dat1)){
  id.i = rownames(dat1)[i]
  for(j in 1:ncol(dat1)){
    gene.j = colnames(dat1)[j]
    dat1[i,j] = length(which(mut1$id==id.i & mut1$Gene.ensGene==gene.j))
  }
}
sum(dat1)
table(mut1$Gene.ensGene %in% gene2kp)

dat1 = as.data.frame(dat1)
dat0 = data.frame(rownames(dat1), cbind(tb.id, tb.id.snv, tb.id.indel), 
                  stringsAsFactors=FALSE)
dim(dat0)
dat0[1:5,]

dat1 = cbind(dat0, dat1)
names(dat1)[1:4] = c("id", "mb", "mb_SNV", "mb_indel")
dim(dat1)
dat1[1:5,1:6]

# ----------------------------------------------------------------------
# write out summary
# ----------------------------------------------------------------------

write.table(dat1, file = "nonsyn_mutation_sample_level.txt", 
            append = FALSE, quote = FALSE, sep = "\t",eol = "\n", 
            row.names = FALSE, col.names = TRUE)

q(save="no")
