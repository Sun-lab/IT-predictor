
# ----------------------------------------------------------------------
# read in mutation calls
# ----------------------------------------------------------------------

setwd("~/research/Immuno/R_batch5/step9")

snvs = read.table("snvs_by_strelka_and_Mutect.txt", header=TRUE, 
                  sep="\t", as.is=TRUE)
dim(snvs)
snvs[1:2,]

indels = read.table("indels_by_strelka.txt", header=TRUE, 
                    sep="\t", as.is=TRUE)
dim(indels)
indels[1:2,]

sort(table(snvs$id), decreasing=TRUE)[1:10]
sort(table(indels$id), decreasing=TRUE)[1:10]

cols = c("seqnames", "start", "end", "REF", "ALT", "id")

input = rbind(snvs[,cols], indels[,cols])
dim(input)
input[1:2,]

# ----------------------------------------------------------------------
# read output 
# ----------------------------------------------------------------------

setwd("~/research/Immuno/R_batch5/step10")
f1 = "snvs_indels_by_strelka_and_Mutect.hg38_multianno.txt"
annotated = read.table(f1, sep="\t", header=TRUE, as.is=TRUE, na.string=".", 
                       quote="")
names(annotated)

dim(annotated)
annotated[1:2,1:8]

# ----------------------------------------------------------------------
# make sure input and annotated ones are one-to-one matach 
# ----------------------------------------------------------------------

ids.input = apply(input[,1:5], 1, paste, collapse=":")
ids.annot = apply(annotated[,1:5], 1, paste, collapse=":")
table(ids.input == ids.annot)

# ----------------------------------------------------------------------
# table a subset of the many columns of annoations 
# ----------------------------------------------------------------------

setwd("~/research/Immuno/R_batch5/step11")

anno1 = annotated[,c(1:13,25,28,31,34,39,42,75:77,85:90)]
dim(anno1)
anno1[1:2,]

summary(anno1$ExAC_nontcga_ALL)
table(anno1$ExAC_nontcga_ALL >= 1e-4)
table(anno1$ExAC_nontcga_ALL >= 5e-4)
table(anno1$ExAC_nontcga_ALL >= 0.001)
table(anno1$ExAC_nontcga_ALL >= 0.005)
table(anno1$ExAC_nontcga_ALL >= 0.010)

pdf("figures/ExAC_VAF.pdf", width=5, height=4)
hist(log10(anno1$ExAC_nontcga_ALL), main="", xlab="ExAC_nontcga_ALL")
dev.off()

table(is.na(anno1$avsnp150))
summary(anno1$ExAC_nontcga_ALL[which(! is.na(anno1$avsnp150))])

get1 <- function(v){ 
  if(is.na(v)){ res = NA }else{
    res = strsplit(v, split=";")[[1]][1]
  }
  res
}

table(anno1$Func.ensGene)
table(anno1$ExonicFunc.ensGene)

table(anno1$Func.ensGene=="exonic", is.na(anno1$ExonicFunc.ensGene))

w2check = which(anno1$Func.ensGene!="exonic"& !is.na(anno1$ExonicFunc.ensGene))
table(anno1$Func.ensGene[w2check])
table(anno1$ExonicFunc.ensGene[w2check])

# ----------------------------------------------------------------------
# generate a new file with mutation call status from Mutect
# as well as some annoation columns
# ----------------------------------------------------------------------

dim(snvs)
snvs[1:2,]

dim(indels)
indels[1:2,]

match(names(indels), names(snvs))

indels.fill = matrix(NA, nrow=nrow(indels), ncol=ncol(snvs)-ncol(indels))
indels = cbind(indels, indels.fill)
names(indels) = names(snvs)

snvs.indels = rbind(snvs, indels)
dim(snvs.indels)

snvs.indels = cbind(snvs.indels, anno1[,-(1:5)])

dim(snvs.indels)
snvs.indels[1:2,]

write.table(snvs.indels, 
            file = "snvs_indels_by_strelka_and_Mutect_with_anno.txt", 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", 
            row.names = FALSE, col.names = TRUE)

sessionInfo()
q(save="no")

