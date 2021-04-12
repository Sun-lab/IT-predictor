
## ---------------------------------------------------------------------
## get sample names for Riaz data
## ---------------------------------------------------------------------


srr = scan("/pine/scr/l/y/lyzhou/riaz_bam/fastq/SRRList.txt", 
           what=character())
length(srr)

srr = sort(srr)
srr[1:5]

setwd('/pine/scr/l/y/lyzhou/riaz_bam/fastq/')
fqfiles1 = system('ls SRR*_1.fastq', intern=TRUE)
fqfiles2 = system('ls SRR*_2.fastq', intern=TRUE)
length(fqfiles1)
length(fqfiles2)

fq.srr1 = gsub("_1.fastq", "", fqfiles1, fixed=TRUE)
fq.srr2 = gsub("_2.fastq", "", fqfiles2, fixed=TRUE)
fq.srr1[1:5]

table(fq.srr1 == fq.srr2)
table(fq.srr1 == srr)

## ---------------------------------------------------------------------
## HLA-HD:
##
## 1. download hla-hd from https://www.genome.med.kyoto-u.ac.jp/HLA-HD/ 
## 2. $ module load bowtie2 
##    (bowtie2/2.2.5 being used)
## 3. Uncompress the downloaded tar.gz file by: 
##    $ tar -zxvf hlahd.version.tar.gz
## 4. Then, move to the uncompressed directory and type: sh install.sh 
## 5. After the installation, add the current directory to your PATH (important):
##    $ export PATH=$PATH:/path_to_HLA-HD_install_directory/bin
## 6. Before running the HLA-HD, check the value of open files on your computer by typing:
##    $ ulimit -Sa
##    If open files are less than 1024, please type:
##    $ ulimit -n 1024
## 7. Run following R and submit jobs 
##
## ---------------------------------------------------------------------

setwd('/pine/scr/w/e/weisun/riaz_bam/')

# export PATH=$PATH:/pine/scr/w/e/weisun/riaz_bam/hlahd.1.2.0.1/bin

fqdir = '/pine/scr/l/y/lyzhou/riaz_bam/fastq/'
outputdir = '/pine/scr/l/y/lyzhou/riaz_bam/hlahd_output/'
shfile = "step13_hlahd.sh"
cat("", file=shfile)
cmd1 = "hlahd.sh -t 4 -m 100 -c 0.95 -f freq_data/"

for(sami in srr){
  cmd2 = sprintf("%s %s%s_1.fastq %s%s_2.fastq", cmd1, fqdir, sami, fqdir, sami)
  cat(sprintf('sbatch -N 1 -n 4 --wrap=\"%s HLA_gene.split.txt dictionary/ %s %s\"\n', cmd2, sami, outputdir), 
      file=shfile, append = TRUE)
}

## ---------------------------------------------------------------------
## Optitype:
## 
## 1. download 
##      OptiType: https://github.com/FRED-2/OptiType
##      Razer3 (64bit): http://packages.seqan.de/razers3/
## 2. module load software (Version being used):
##      Python             (Python/3.6.5-foss-2016b-fh3)
##      SAMtools           (samtools/1.0)
##      HDF5               (HDF5/1.8.18-foss-2016b)
##      GLPK               (GLPK/4.61-foss-2016b)
## 3. create and activate virture env for python
##      $ python3 -m venv /path/to/new/virtual/environment
##      $ source /path/to/new/virtual/environment/bin/activate
##      $ pip install numpy
##      $ pip install pyomo
##      $ pip install pysam
##      $ pip install matplotlib
##      $ pip install tables
##      $ pip install pandas
##      $ pip install future
## 4. create a configuration file config.ini in the same direcoty as OptiType
##      can be found at /fh/fast/sun_w/licai/HLA/OptiType/config.ini
## 5. Run following R and submit jobs 
##
## ---------------------------------------------------------------------

# cd /fh/fast/sun_w/licai/HLA/

for(sami in sams){
  fnm = paste0(sami, '.sh')
  # create shell script for each sample
  write('#! /bin/bash',file = fnm )
  write('source ~/mypy/bin/activate', file = fnm, append = T) # activate viture env
  write(sprintf('python ./OptiType/OptiTypePipeline.py -i %s%s_1.fastq %s%s_2.fastq --dna -v -o ./optitype_out/%s \n',
                fqdir, sami, fqdir, sami, sami), 
        file = fnm, append = T)
  # cat(sprintf('sbatch --exclusive %s.sh \n', sami))
  cat(sprintf('sbatch --mem=768000 --partition=largenode -c 6 -t 0-12 %s.sh \n', sami))
}
