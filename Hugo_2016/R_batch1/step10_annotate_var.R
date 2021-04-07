
# ----------------------------------------------------------------------
# read in mutation calls
# ----------------------------------------------------------------------

setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

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

write.table(input, file = "snvs_indels_by_strelka_and_Mutect.avinput", 
            append = FALSE, quote = FALSE, sep = "\t",eol = "\n", 
            row.names = FALSE, col.names = FALSE)

# ----------------------------------------------------------------------
# run it 
# ----------------------------------------------------------------------

cmd = "/fh/fast/sun_w/bin/annovar/table_annovar.pl"
cmd = paste0(cmd, " snvs_indels_by_strelka_and_Mutect.avinput")
cmd = paste0(cmd, " /fh/fast/sun_w/bin/annovar/humandb/ -buildver hg38")
cmd = paste0(cmd, " -out snvs_indels_by_strelka_and_Mutect")
cmd = paste0(cmd, "  -remove -protocol ensGene,dbnsfp33a,exac03nontcga")
cmd = paste0(cmd, ",avsnp150,clinvar_20170905")
cmd = paste0(cmd, " -operation g,f,f,f,f -nastring .")

cmd

date()
system(cmd)
date()

q(save="no")

