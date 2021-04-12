
sams = scan("/pine/scr/l/y/lyzhou/riaz_bam/fastq/SRRList.txt", 
            what=character())
length(sams)

sams = sort(sams)
sams[1:5]

## ---------------------------------------------------------------------
## Organize HLA-HD result
## ---------------------------------------------------------------------

setwd('/pine/scr/w/e/weisun/riaz_bam/hlahd_output')
hlahd_out = list()
hla_comp = NULL

for(sami in sams){
  filei = list.files(path = sami, pattern = '_final.result.txt', recursive = T)  
  filei
  filenmn = sprintf("%s/%s", sami, filei)
  ncol = max(count.fields(filenmn, sep = "\t"))
  hlai = read.table(filenmn, header=F, sep='\t', row.names = 1, as.is = T,
                    fill=T, col.names=paste0('V', seq_len(ncol)))
  
  hlai[,2] = apply(hlai, 1, function(x) if(x[2] == '-') x[1] else x[2])
  hlahd_out[[sami]] = hlai
}

sort(table(unlist(hlahd_out)), decreasing=TRUE)

save(hlahd_out, file = 'hlahd_result.Rdata')

sessionInfo()

q(save ='no')
