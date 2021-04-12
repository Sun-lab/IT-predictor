
setwd('/fh/fast/sun_w/licai/HLA/')

sams = read.csv('_hugo_samples.csv', as.is =T)$x
sams

## ---------------------------------------------------------------------
## Organize optitype result
## ---------------------------------------------------------------------

setwd('./optitype_out')
optitype_out = NULL
sam_null = NULL
for(sami in sams){
  filei = list.files(path = sami, pattern = '_result.tsv', recursive = T)  
  filei
  if(length(filei) ==0L){
    print(sami)
    sam_null = c(sam_null, sami)
    next
  }
  hlai = read.table(sprintf("%s/%s", sami, filei), as.is = T, header = T)
  rownames(hlai) = sami
  optitype_out = rbind(optitype_out, hlai)
}
# sapply(1:6, function(x) hla_out[,x] <- gsub('[A-C][*]', '', hla_out[,x]))
write.table(optitype_out, file = 'optitype_result.txt', sep = '\t', 
            col.names = T, row.names = T, quote = F )
sam_null

## ---------------------------------------------------------------------
## Organize HLA-HD result
## ---------------------------------------------------------------------

setwd('../hlahd.1.2.0.1/estimation')
hlahd_out = list()
hla_comp = NULL
for(sami in sams){
  filei = list.files(path = sami, pattern = '_final.result.txt', recursive = T)  
  filei
  filenmn = sprintf("%s/%s", sami, filei)
  ncol <- max(count.fields(filenmn, sep = "\t"))
  hlai = read.table(filenmn, header=F, sep='\t', row.names = 1, as.is = T,
                    fill=T, col.names=paste0('V', seq_len(ncol)))
  
  hlai[,2] = apply(hlai, 1, function(x) if(x[2] == '-') x[1] else x[2])
  hlahd_out[[sami]] = hlai
}

save(hlahd_out, file = 'hlahd_result.Rdata')

## ---------------------------------------------------------------------
## Compare two result: 4 digit
## ---------------------------------------------------------------------

hlahd_out2 = t(sapply(hlahd_out, function(x) 
  substr(c(x[1, 1:2], x[2, 1:2], x[3, 1:2]),7,11)))
colnames(hlahd_out2) = paste0(rep(LETTERS[1:3],each=2), 1:2)

optitype_out2 = sapply(optitype_out[,1:6], substr, 3, 7)
rownames(optitype_out2) = rownames(optitype_out)
hlahd_out2[1:5,]
optitype_out2[1:5,]

table(rownames(hlahd_out2) == rownames(optitype_out2))

for(coli in 0:2){
  indi = which(hlahd_out2[,coli*2+1] == optitype_out2[,coli*2+2] & 
               hlahd_out2[,coli*2+2] == optitype_out2[,coli*2+1] )
  hlahd_out2[indi, coli*2+1:2] = hlahd_out2[indi, coli*2+2:1]
}

compare_res = (hlahd_out2 == optitype_out2)
apply(compare_res, 2, sum)[1:3*2]

ind2 = which(apply(compare_res, 1, all) == F)
hlahd_out2[ind2,]
optitype_out2[ind2,]

q(save ='no')
