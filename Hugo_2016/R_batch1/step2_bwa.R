
setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

srr1 = scan(file="SRR_Acc_List.txt", what=character(0))
srr2 = scan(file="SRR_Acc_List_Vanderbilt.txt", what=character(0))

srr1[1:5]
srr2[1:5]

srr = c(srr1, srr2)
nn = length(srr)

setwd("/fh/fast/sun_w/research/Immuno/R_batch2")

idxf = "/fh/fast/sun_w/research/data/human/hg38/Homo_sapiens_assembly38.fasta "

cmd1 = "/fh/fast/sun_w/bin/bwa-0.7.17/bwa"
errs = paste0("--error=bwa_", srr, ".txt")
outs = paste0("--output=bwa_", srr, ".txt")
path = "/fh/fast/sun_w/research/data/ncbi_sra/public/fastq"
inp1 = sprintf("%s/%s_1.fastq.gz", path, srr)
inp2 = sprintf("%s/%s_2.fastq.gz", path, srr)
oupt = sprintf("/fh/fast/sun_w/research/Immuno/data/Hugo_2016/sams/%s.sam", srr)

cmds = "sbatch --constraint x10sle --comment=\"naccesspolicy:singlejob\""
cmds = sprintf("%s -n 1 -c 4 %s %s --wrap=\"%s", cmds, errs, outs, cmd1)
cmds = sprintf("%s mem -t 4 %s %s %s > %s\"", cmds, idxf, inp1, inp2, oupt)

cat(cmds[1])

cat(cmds, file="step2_bwa_batch.sh", sep="\n")

q(save="no")
