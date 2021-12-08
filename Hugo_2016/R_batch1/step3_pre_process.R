
setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

srr1 = scan(file="SRR_Acc_List.txt", what=character(0))
srr2 = scan(file="SRR_Acc_List_Vanderbilt.txt", what=character(0))

srr1[1:5]
srr2[1:5]

srr = c(srr1, srr2)

nn = length(srr)

setwd("/fh/fast/sun_w/research/Immuno/R_batch2")

codes = scan("step3_pre_process_ex.sh", what=character(), sep="\n", 
             blank.lines.skip=FALSE)
length(codes)
codes[1:5]

shells=paste0("step3/step3_pre_process_", srr, ".sh")

for(i in 1:nn){
  srr1 = srr[i]
  codes[3] = sprintf("sampleName=\"%s\"", srr1)
  cat(codes, file=shells[i], sep="\n")
  cat("\n", file=shells[i], append=TRUE)
}

logFile = paste0("prep_", srr, ".txt")
cmds = "sbatch -n 1 -c 2 --mem=16192 --error="
cmds = paste0(cmds, logFile, " --output=", logFile, " ", shells)

cat(cmds, file="step3_pre_process_batch.sh", sep="\n")

q(save="no")
