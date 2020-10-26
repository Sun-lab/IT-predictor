library(data.table)
library(R.utils)  

# --------------------- read in TPM data ---------------------  
fnm1 = "GSE120575_Sade_Feldman_melanoma_single_cells_nolabel_1.txt.gz"
fnm2 = "GSE120575_Sade_Feldman_melanoma_single_cells_2.txt.gz"
sf_tpm1 = fread(paste0("../scRNAseq/Sade_Feldman_2018/", fnm1), fill = TRUE, header = TRUE)
sf_tpm2 = fread(paste0("../scRNAseq/Sade_Feldman_2018/", fnm2), fill = TRUE, drop = 16293, col.names = colnames(sf_tpm1))
ls_tpm = list(sf_tpm1,sf_tpm2)
rm(sf_tpm1)
rm(sf_tpm2)
sf_tpm = rbindlist(ls_tpm)
rm(ls_tpm)

# --------------------- read in cell type data ---------------------  
cell_type = read.table("cell_type.txt",header = TRUE,sep = "\t",as.is=TRUE)
dim(cell_type)
cell_type[1:2,]
table(cell_type$type)
clustered_cells = cell_type$cell[cell_type$type != "unclustered"]
selected.cols = c("V1", clustered_cells)
cluster_sf_tpm = sf_tpm[, ..selected.cols] # 55737 x 14778
cluster_sf_tpm[1:5,1:5]
rm(sf_tpm)
rm(clustered_cells)
rm(selected.cols)

# --------------------- Select ~10,000 genes expressed in at least a proportion of cells ---------------------
gene.grouping = rep(1:28, each = 2000)[1:dim(cluster_sf_tpm)[1]]
cells.per.gene = NULL
for (i in 1:28) {
  temp.counts = rowSums(cluster_sf_tpm[gene.grouping==i,-1] > 0)
  cells.per.gene = c(cells.per.gene,temp.counts)
}
large.express.genes.all = cluster_sf_tpm$GeneSymbol[which(cells.per.gene > 500)] # kept 12341 genes
rm(gene.grouping)
clustered_tpm_LEgenes = cluster_sf_tpm[cluster_sf_tpm$GeneSymbol %in% large.express.genes.all, ] 
clustered_tpm_LEgenes[1:5,1:5]

# --------------------- Make Reference sample file ---------------------
typesOFcells = cell_type$type[match(colnames(clustered_tpm_LEgenes)[-1], cell_type$cell)]
write.table(clustered_tpm_LEgenes, "CIBERSORTx_ref_sample.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = c("GeneSymbol", typesOFcells))


