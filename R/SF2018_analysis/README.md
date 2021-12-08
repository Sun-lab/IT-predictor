# Reanalysis of single cell RNA-seq data from Sade-Feldman et al. (2018)

We downloaded the scRNA-seq data in TPM format from NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575, including two files

* GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz	(120.9 Mb)
* GSE120575_patient_ID_single_cells.txt.gz	(81.1 Kb)

## 1_SF2018_CD8Tcells: Refine the CD8 T cell clustering.

Re-clustering the CD8 T cells using a simple approach (PCA followed by K-means) and compare the results versus Sade-Feldman et al. (2018). We take the intersection of the genes annotated to the two cluster CD8_G and CD8_B by Sade-Feldman et al. (2018), and to the 3 clusters of our K-means results. Here are the description of the CD8_G and CD8_B clusters from Sade-Feldman et al. (2018):

> Clustering all CD8+ T cells (n = 6,350) revealed two major cell states: CD8_G with increased expression
> of genes linked to memory, activation, and cell survival (e.g., IL7R, TCF7, REL, FOXP1, FOSL2, and STAT4)
> and reduced expression of co-inhibitory molecules; and CD8_B enriched for genes linked to cell exhaustion 
(e.g., CD38, HAVCR2, ENTPD1, PDCD1, BATF, LAG3, CTLA4, and PTPN6)

## 2_Sade_Feldman_2018_hvg: Refine the clustering of all the cells.

Re-clustering the CD8 T cells using a simple approach (PCA followed by K-means) and compare the results versus Sade-Feldman et al. (2018). We also combine the definition of CD8_G and CD8_B cells from ```1_SF2018_CD8Tcells```, and seek to re-clustering NK cells separately. 

## 3_SF2018_cell_type_correlation: check the correlation between SF reference and LM22 reference used by CIBERSORT

Since CIBERSORT (Newman et al. 2015) is a popular software package for cell type composition estimation and its signature matrix LM22, the gene expression of 22 immune cell types have been used in many earlier studies, we compare the two signature matrix. Note here we did not select the best marker gene set from SF 2018 data and just use the genes in the LM22 matrix (around 300 genes after excluding those that are lowly expressed in SF 2018 data). We do observed overall consistent gene expression pattern, i.e., gene expression of the same or similar cell types tend to have higher expression. 

## 4a_SF2018_DEanalysis /  4b_SF2018_DEanalysis_filtered_cells

Differential expression analysis to compare one cell type vs. all the other cells, using all the cells or top 200 cells with largest proportion of expressed genes (```_filtered_cells```). 

## 5a_SF2018_signature_genes /  5b_SF2018_signature_genes_filtered_cells

Select signature genes for each cell type based on the DE results from step 4. We selected around 60 genes with smallest p-value, and largest fold changes. Technically, we define fold changes as the average expression from other cell types versus the average expression of the cell type of interest, and thus selected those with smallest fold changes. 

## Reference

Sade-Feldman et al. (2018). "Defining T cell states associated with response to checkpoint immunotherapy in melanoma." Cell 175.4 (2018): 998-1013.

Newman et al. (2015) "Robust enumeration of cell subsets from tissue expression profiles." Nature methods 12.5 (2015): 453-457.
