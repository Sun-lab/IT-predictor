
#--------------------------------------------------------------------
# Step 5: check the ouptut of MHC-I or MHC-II predictions
#--------------------------------------------------------------------

library(stringr)
library(data.table)

#--------------------------------------------------------------------
# all the output files
#--------------------------------------------------------------------

outs = list.files("../output", pattern="txt.gz", 
                  full.names="TRUE")

length(outs)
outs[1:3]

output_info = gsub("../output/", "", outs)
output_info = gsub(".txt.gz", "", output_info)

output_info = strsplit(output_info, "_hla")

table(sapply(output_info, length))

output_info = data.frame(matrix(unlist(output_info), byrow = TRUE, ncol=2))
names(output_info) = c("sample", "label")
dim(output_info)
output_info[1:2,]

table(table(output_info$sample))
table(output_info$label)

#--------------------------------------------------------------------
# read in HLA-I results
#--------------------------------------------------------------------

outs_hla_i = outs[grep("_hlai_", outs)]
length(outs_hla_i)

hla_i_nms = gsub("../output/|.txt.gz", "", outs_hla_i)
hla_i_nms[1:4]

paste1 <- function(x,y){
  for(k in 2:length(x)){
    if(x[k] == "" && y[k] != "Ave"){ x[k] = x[k-1]}
  }
  x[x!=""] = paste0(x[x!=""], "_")
  paste0(x, y)
}

dat_hla_i = list()

for(i in 1:length(outs_hla_i)){
  input_i  = outs_hla_i[i]
  header_i = fread(input_i, nrows=2)
  header_i = as.matrix(header_i)
  header_i = paste1(header_i[1,], header_i[2,])
  
  dat_i = fread(input_i, skip=2)
  names(dat_i) = header_i
  dim(dat_i)
  dat_i[1:2,]
  
  cols2check = grep("_core|_icore", header_i)
  for(c1 in cols2check){
    if(any(dat_i$Peptide != dat_i[[c1]])){
      stop("core or icore does not match with peptide\n")
    }
  }
  
  dat_i = dat_i[,(names(dat_i)[cols2check]):=NULL]
  
  dat_hla_i[[i]] = dat_i
}

names(dat_hla_i) = hla_i_nms

length(dat_hla_i)
summary(sapply(dat_hla_i, nrow))
table(sapply(dat_hla_i, ncol))

saveRDS(dat_hla_i, file = "../data/netMHCpan4_1_results.rds")

#--------------------------------------------------------------------
# read in HLA-II results
#--------------------------------------------------------------------

outs_hla_ii = outs[grep("_hlaii_", outs)]
length(outs_hla_ii)

hla_ii_nms = gsub("../output/|.txt.gz", "", outs_hla_ii)
hla_ii_nms[1:4]

paste2 <- function(x,y){
  for(k in 2:length(x)){
    if(x[k] == "" && y[k] == "Rank"){ x[k] = x[k-1]}
  }
  x[x!=""] = paste0(x[x!=""], "_")
  paste0(x, y)
}

dat_hla_ii = list()


for(i in 1:length(outs_hla_ii)){
  input_i  = outs_hla_ii[i]
  header_i = fread(input_i, nrows=2, fill=TRUE)
  header_i = as.matrix(header_i)
  header_i = paste2(header_i[1,], header_i[2,])
  
  dat_i = fread(input_i, skip=2)
  names(dat_i) = header_i
  dim(dat_i)
  dat_i[1:2,]
  
  dat_hla_ii[[i]] = dat_i
}

names(dat_hla_ii) = hla_ii_nms

length(dat_hla_ii)

summary(sapply(dat_hla_ii, nrow))
table(sapply(dat_hla_ii, ncol))

saveRDS(dat_hla_ii, file = "../data/netMHCIIpan4_0_results.rds")

sessionInfo()
q(save="no")
