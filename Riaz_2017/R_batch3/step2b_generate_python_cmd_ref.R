
library(stringr)

models = list.files("../../../PEPPRMINT/results/PEPPRMINT", pattern = ".h5")

models

n1 = as.numeric(str_extract(models, "(?<=MA_)\\d+"))
s1 = as.numeric(str_extract(models, "(?<=split)\\d+"))

n1
s1

fnm = "../../../PEPPRMINT/python/run_files/Riaz_2017_ref.sh"
cat("", file=fnm)

for(i in 1:length(models)){
  m1 = models[i]
  mtag = sprintf("%d_split%d", n1[i], s1[i])
  log_file = sprintf("logfiles/PEPPRMINT_Riaz_2017_%s_ref.log", mtag)
  
  code = "python3 _prediction_PEPPRMINT.py \\\n"
  code = paste0(code, "--input_test_pred riaz_peptide_ref.txt \\\n")
  code = paste0(code, "--input_alist riaz_HLA_I_allele_list.txt \\\n")
  code = paste0(code, "--test_data_name PEPPRMINT_Riaz_2017 \\\n")
  code = paste0(code, "--data_dir ../data/Riaz_2017 \\\n")
  code = paste0(code, "--results_dir ../results/Riaz_2017 \\\n")
  code = paste0(code, "--test_data_name PEPPRMINT_Riaz_2017_ref \\\n")
  code = paste0(code, "--m_tag ", mtag, " \\\n")
  code = paste0(code, "--model ", m1, " \\\n")
  code = paste0(code, "--neoantigen \\\n")
  code = paste0(code, "> ", log_file, " \n\n")
  cat(code, file = fnm, append = TRUE)
}


sessionInfo()
q(save="no")



