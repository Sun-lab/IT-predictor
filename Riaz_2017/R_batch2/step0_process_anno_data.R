
#--------------------------------------------------------------------
# Process annotated Riaz data to extract sequences for making peptides 
# for somatic mutatins
#
# Data: snvs_indels_by_strelka_and_Mutect_with_anno_filtered.txt
# Steps: 
  #1. Import data
  #2. Filter dataset
  #3. Extract the mutation data on the amino acid and DNA level
  #4. Extract ensembletrans id
  #5. save dataset as riaz_mutdata.RData
#--------------------------------------------------------------------

aas = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

#--------------------------------------------------------------------
# Step 1: Import Data
#--------------------------------------------------------------------

library(data.table)
library(stringr)

input = "../data/snvs_indels_by_strelka_and_Mutect_with_anno_filtered.txt"
riaz_mutdata = fread(input)

dim(riaz_mutdata)
names(riaz_mutdata)

riaz_mutdata[1:10, 1:5]
riaz_mutdata[1:10, c(1,2,78,97)]

length(unique(riaz_mutdata$id))
tb1 = table(riaz_mutdata$id)
length(tb1)
sort(tb1)

patient_info = data.frame(tb1)
names(patient_info) = c("sample", "total_mutations")
dim(patient_info)
patient_info[1:2,]

patient_info$patient = str_extract(patient_info$sample, "\\S+(?=_)")
patient_info$time    = str_extract(patient_info$sample, "(?<=_)\\S+")
dim(patient_info)
patient_info[1:2,]

table(patient_info$time)
length(unique(patient_info$patient[which(patient_info$time == "pre")]))
length(unique(patient_info$patient[which(patient_info$time == "on")]))

#--------------------------------------------------------------------
# Step 2: Filter dataset 
#--------------------------------------------------------------------

# chose nonsynomous mutations
riaz_mutdata0 = riaz_mutdata[which(riaz_mutdata$nonSynonymous==TRUE), ]
dim(riaz_mutdata0)

riaz_mutdata0[1:2, c(1,2,78,97)]

# check NA in AAChange.ensGene
table(is.na(riaz_mutdata0$AAChange.ensGene))

# observe that these are the cases where Func.ensGene == "splicing"
table(riaz_mutdata0$Func.ensGene[is.na(riaz_mutdata0$AAChange.ensGene)])

riaz_mutdata1 = riaz_mutdata0[!is.na(riaz_mutdata0$AAChange.ensGene),]
dim(riaz_mutdata1)

riaz_mutdata1[1:2, c(1,2,78,97)]
riaz_mutdata1 = as.data.frame(riaz_mutdata1)

table(riaz_mutdata1$id)

#--------------------------------------------------------------------
# Step 3: Extract mutation data 
#--------------------------------------------------------------------

# get the wildtype amino acid, mutation position, variant amino acid
# by default of str_match, get the first match only
mutamino = str_match(riaz_mutdata1$AAChange.ensGene, "(p.)([A-Z])(\\d+)([A-Z])")
mutseq   = str_match(riaz_mutdata1$AAChange.ensGene, "(c.)([A-Z])(\\d+)([A-Z])")
mut = cbind(mutamino, mutseq)
mut = as.data.frame(mut)
dim(mut)
mut[1:5,]

table(is.na(mut$V1))
table(is.na(mut$V6))

# the amino acid changes were not annotated for indels
# this is default of ANNOVA, thoug there is a way around it. check 
# http://annovar.openbioinformatics.org/en/latest/user-guide/gene/
table(is.na(mut$V1), riaz_mutdata1$type, useNA="ifany")

mut           = mut[,c(3,4,5,8,9,10)]
colnames(mut) = c("W", "pos", "M", "Wseq", "seqpos", "Mseq")
mut$pos       = as.numeric(as.character(mut$pos))
mut$seqpos    = as.numeric(as.character(mut$seqpos))
riaz_mutdata2 = cbind(riaz_mutdata1, mut)
riaz_mutdata2 = riaz_mutdata2[which(riaz_mutdata1$type == "SNV"),]
dim(riaz_mutdata2)
riaz_mutdata2[1:2, c(1,2,97:103)]

#--------------------------------------------------------------------
# Step 4: Extract ENST
#--------------------------------------------------------------------

mutensembl = str_match(riaz_mutdata2$AAChange.ensGene, "(ENST)(\\d+)")
mutensembl = mutensembl[,1]
mutensembl = data.frame(EnsembleID=mutensembl, stringsAsFactors=FALSE)
riaz_mutdata = cbind(riaz_mutdata2, mutensembl)

dim(riaz_mutdata)
riaz_mutdata[1:2, c(1,2,78,97:104)]

#--------------------------------------------------------------------
# Step 5: check amino acid annoation, remove those 
# where the whild type amino acid is "X"
#--------------------------------------------------------------------

table(riaz_mutdata$W %in% aas)
table(riaz_mutdata$M %in% aas)
table(riaz_mutdata$M[!riaz_mutdata$M %in% aas])

riaz_mutdata[! riaz_mutdata$W %in% aas, c(1,9:13,78,97:104)]
riaz_mutdata[! riaz_mutdata$M %in% aas, c(1,9:13,78,97:104)][1:2,]

riaz_mutdata = riaz_mutdata[which(riaz_mutdata$W %in% aas),]
dim(riaz_mutdata)

#--------------------------------------------------------------------
# Step 6. Save RData 
#--------------------------------------------------------------------

length(unique(riaz_mutdata$id))
tb2 = table(riaz_mutdata$id)
length(tb2)
sort(tb2)

patient_info$n_mutations_neoAg = rep(0, nrow(patient_info))
match1 = match(names(tb2), patient_info$sample)
table(is.na(match1))

patient_info$n_mutations_neoAg[match1] = tb2
dim(patient_info)
patient_info[1:2,]

write.table(patient_info, file = "../data/riaz_patient_mb_info.txt", 
            quote=FALSE, col.names= TRUE, row.names = FALSE)

save(riaz_mutdata, file="../data/riaz_mutdata.RData")

sessionInfo()
q(save="no")

