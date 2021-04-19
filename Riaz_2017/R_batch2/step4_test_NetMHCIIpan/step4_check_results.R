
# check whether the results of NetMHCIIpan is complete

library(seqinr)
library(data.table)

pep_file = "../../data/netmhcii/Pt89_on_mut.txt"

hla = c("DRB1_0403", "DRB1_0901", 
        "HLA-DQA10301-DQB10303", "HLA-DQA10302-DQB10302", 
        "HLA-DQA10301-DQB10302", "HLA-DQA10302-DQB10303", 
        "HLA-DPA10103-DPB10201", "HLA-DPA10103-DPB10401", 
        "HLA-DPA10103-DPB10401", "HLA-DPA10103-DPB10201")

pep = read.fasta(pep_file, seqtype="AA", as.string=TRUE)
is.list(pep)
length(pep)
pep[1]

res0 = fread("Pt89_on_hlaii_mut_v0.txt", fill=TRUE)
dim(res0)
res0[1:2,]

table(res0$allele)

res0_sub = res0[which(res0$allele == "HLA-DPA10103-DPB10201"),]
dim(res0_sub)

dim(unique(res0_sub))

paste1 <- function(x,y){
  for(k in 2:length(x)){
    if(x[k] == "" && y[k] == "Rank"){ x[k] = x[k-1]}
  }
  x[x!=""] = paste0(x[x!=""], "_")
  paste0(x, y)
}

header = fread("Pt89_on_hlaii_mut.txt", nrows=2, fill=TRUE)
header
header = as.matrix(header)
header = paste1(header[1,], header[2,])
header

res1 = fread("Pt89_on_hlaii_mut.txt", skip=2, header=FALSE)
names(res1) = header
dim(res1)
res1[1:2,]

sessionInfo()
q(save="no")



