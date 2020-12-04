# Reanalysis of single cell RNA-seq data from Sade-Feldman et al. (2018)

We downloaded the scRNA-seq data in TPM format from NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575, including two files

* GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz	(120.9 Mb)
* GSE120575_patient_ID_single_cells.txt.gz	(81.1 Kb)

## 1_SF2018_CD8Tcells: Refine the CD8 T cell clustering.

Re-clustering the CD8 T cells using a simple approach (PCA followed by K-means) and compare the results versus Sade-Feldman et al. (2018). We take the intersection of the genes annotated to the two cluster CD8_G and CD8_B by Sade-Feldman et al. (2018), and to the 3 clusters of our K-means results. Here are the description of the CD8_G and CD8_B clusters from Sade-Feldman et al. (2018):

> Clustering all CD8+ T cells (n = 6,350) revealed two major cell states: CD8_G with increased expression> of genes linked to memory, activation, and cell survival (e.g., IL7R, TCF7, REL, FOXP1, FOSL2, and STAT4)
> and reduced expression of co-inhibitory molecules; and CD8_B enriched for genes linked to cell exhaustion 
(e.g., CD38, HAVCR2, ENTPD1, PDCD1, BATF, LAG3, CTLA4, and PTPN6)

## 2_Sade_Feldman_2018_hvg: Refine the clustering of all the cells.

Re-clustering the CD8 T cells using a simple approach (PCA followed by K-means) and compare the results versus Sade-Feldman et al. (2018). We also combine the definition of CD8_G and CD8_B cells from ```1_SF2018_CD8Tcells```, and seek to re-clustering NK cells separately. 


## Reference

Sade-Feldman et al. (2018). "Defining T cell states associated with response to checkpoint immunotherapy in melanoma." Cell 175.4 (2018): 998-1013.
