
setwd("/fh/fast/sun_w/research/Immuno/data/Hugo_2016")

srr1 = scan(file="SRR_Acc_List.txt", what=character(0))
srr2 = scan(file="SRR_Acc_List_Vanderbilt.txt", what=character(0))

srr1[1:5]
srr2[1:5]

srr = c(srr1, srr2)

# seems sratoolkit on sever is too old, use my own version located at
# /fh/fast/sun_w/bin/sratoolkit.2.8.2-1-ubuntu64/bin/
# first run /fh/fast/sun_w/bin/sratoolkit.2.8.2-1-ubuntu64/bin/vdb-config -i 
# to configure, set working directory, and add dbGap repository key
# RUN THE SCRIPT IN THE WORKING DIRECTORY

setwd("/fh/fast/sun_w/research/Immuno/R_batch2")

cmd1 = "/fh/fast/sun_w/bin/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump"
cmds = paste0("sbatch --wrap=\"", cmd1, " --split-files ", srr, "\"")

cat(cmds, file="step1_download_files_Hugo_2016.sh", sep="\n")

q(save="no")
