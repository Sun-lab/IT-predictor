
#--------------------------------------------------------------------
# Process annotated Hugo data to extract sequences for making peptides 
# for somatic mutatins, same for Hugo et al data in 'mhcpeptide_hugo_mutdata3.R'
#
# Data: snvs_indels_by_strelka_and_Mutect_with_anno_filtered.txt
# Steps: 
  #1. Import data
  #2. Filter dataset
  #3. Extract the mutation  data on the amino acid level and sequence level
  #4. Extract ensembletrans id
  #5. save dataset as riaz_hugo_mutdata3.RData
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Step 1: Import Data
#--------------------------------------------------------------------

library(readr)
library(stringr)

input = "snvs_indels_by_strelka_and_Mutect_with_anno_filtered_corrected.txt"
#164831
hugo_mutdata = read_delim(input, delim="\t", col_names=TRUE)
# 4 parsing errors: for Pt17 observation, variable GTEx_V6_gene is not NA
hugo_mutdata[which(hugo_mutdata$id=="Pt17" & hugo_mutdata$QSS==20 & !is.na(hugo_mutdata$GTEx_V6_gene)), "GTEx_V6_gene"]
# I think the program just replaces it with NA .i.e. if even if I manually replace with NA it is still 164831
# only two subjects anyway


dim(hugo_mutdata)
hugo_mutdata[1:10, 1:5]
hugo_mutdata[1:10, c(1,2,78,97)]

# there maybe some duplicated rows due to file processing in prevoius steps
du = duplicated(hugo_mutdata)
which(du)
names(hugo_mutdata)

length(unique(hugo_mutdata$id))
tb1 = table(hugo_mutdata$id)
length(tb1)
sort(tb1)

patient_info = data.frame(tb1)
names(patient_info) = c("sample", "total_mutations")
dim(patient_info)
patient_info[1:2,]

#--------------------------------------------------------------------
# Step 2: Filter dataset 
#--------------------------------------------------------------------
#--------subset to where nonsynomous== TRUE

hugo_mutdata0 = hugo_mutdata[which(hugo_mutdata$nonSynonymous==TRUE), ]
dim(hugo_mutdata0)

hugo_mutdata0[1:10, c(1,2,78,97)]

# check NA in AAChange.ensGene
table(is.na(hugo_mutdata0$AAChange.ensGene))

# observe that these all are  where Func.ensGene == "splicing"
table(hugo_mutdata0$Func.ensGene[is.na(hugo_mutdata0$AAChange.ensGene)])

hugo_mutdata1 = hugo_mutdata0[!is.na(hugo_mutdata0$AAChange.ensGene),]
dim(hugo_mutdata1)

hugo_mutdata1[1:10, c(1,2,78,97)]
hugo_mutdata1 = as.data.frame(hugo_mutdata1)

#--------------------------------------------------------------------
# Step 3: Extract mutation data 
#--------------------------------------------------------------------

# get the wildtype amino acid, mutation position, variant amino acid
# by default of str_match, get the first match only
mutamino = str_match(hugo_mutdata1$AAChange.ensGene, "(p.)([A-Z])(\\d+)([A-Z])")
mutseq   = str_match(hugo_mutdata1$AAChange.ensGene, "(c.)([A-Z])(\\d+)([A-Z])")
mut = cbind(mutamino, mutseq)
mut = as.data.frame(mut)
dim(mut)
mut[1:5,]

table(is.na(mut$V1))
table(is.na(mut$V6))

# the amino acid changes were not annotated for indels
# this is default of ANNOVA, thoug there is a way around it. check 
# http://annovar.openbioinformatics.org/en/latest/user-guide/gene/
table(is.na(mut$V1), hugo_mutdata1$type, useNA="ifany")

mut           = mut[,c(3,4,5,8,9,10)]
colnames(mut) = c("W", "pos", "M", "Wseq", "seqpos", "Mseq")
mut$pos       = as.numeric(as.character(mut$pos))
mut$seqpos    = as.numeric(as.character(mut$seqpos))
hugo_mutdata2 = cbind(hugo_mutdata1, mut)
hugo_mutdata2 = hugo_mutdata2[which(hugo_mutdata1$type == "SNV"),]
dim(hugo_mutdata2)
hugo_mutdata2[1:10, c(1,2,97:103)]

#--------------------------------------------------------------------
# Step 4: Extract ENST
#--------------------------------------------------------------------

mutensembl = str_match(hugo_mutdata2$AAChange.ensGene, "(ENST)(\\d+)")
mutensembl = mutensembl[,1]
mutensembl = data.frame(EnsembleID=mutensembl, stringsAsFactors=FALSE)
hugo_mutdata3 = cbind(hugo_mutdata2, mutensembl)

dim(hugo_mutdata3)
hugo_mutdata3[1:10, c(1,2,78,97:104)]

#--------------------------------------------------------------------
# Step 5. Save RData 
#--------------------------------------------------------------------

length(unique(hugo_mutdata3$id))
tb2 = table(hugo_mutdata3$id)
length(tb2)
sort(tb2)

patient_info$n_mutations_neoAg = rep(0, nrow(patient_info))
match1 = match(names(tb2), patient_info$sample)
table(is.na(match1))

patient_info$n_mutations_neoAg[match1] = tb2
dim(patient_info)
patient_info[1:2,]

write.table(patient_info, file = "hugo_patient_mb_info.txt", 
            quote=FALSE, col.names= TRUE, row.names = FALSE)

save(hugo_mutdata3, file="hugo_mutdata3.RData")

sessionInfo()
q(save="no")

