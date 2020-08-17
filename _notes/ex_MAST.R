# check here: https://www.bioconductor.org/packages/release/bioc/html/MAST.html

library(MAST)

# sim_matrix is a big matris of count data, with one row per gene and one column per cell

sim_matrix_log = log2(1 + sim_matrix) #log transformed data

dim(sim_matrix_log)
sim_matrix_log[1:3, 1:4]
cell_id = colnames(sim_matrix_log)   # get the cell id from the data
gene_id = rownames(sim_matrix_log)   # get the gene id from the data

# replace the diagnosis with cell type label. 
diagnosis = as.character(meta$phenotype) #
diagnosis[diagnosis == 1] = "Case"
diagnosis[diagnosis == 0] = "Control"

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey = cell_id)
length(fData)
length(cData)

dim(meta)
meta[1:2,]

dim(meta_ind)
meta_ind[1:2,]

sca = FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta$CDR)
colData(sca)$diagnosis = as.factor(diagnosis)
colData(sca)$ind = as.factor(meta$individual)
colData(sca)$RIN = meta_ind$RIN[match(meta$individual, meta_ind$individual)]

colData(sca)

getOption("mc.cores")

date()
b0 = zlm(formula = ~ diagnosis + cngeneson + RIN, sca = sca, parallel = TRUE)
date()

lrt0 = lrTest(b0, "diagnosis")

dim(lrt0)
lrt0[1,,]
mast_pval_glm = apply(lrt0, 1, function(x){x[3,3]})

