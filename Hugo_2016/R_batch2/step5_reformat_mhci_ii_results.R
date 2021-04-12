
#--------------------------------------------------------------------
# Step 5: reformat the ouptut of MHC-I or MHC-II predictions
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# MHC-I
#--------------------------------------------------------------------

outs = list.files("for_netmhci/_output", pattern="txt", 
                  full.names="TRUE")

length(outs)
outs[1:3]

for(i in 1:length(outs)){
  print(i)
  input.i  = outs[i]
  output.i = gsub("txt$", "tsv", input.i)
  print(output.i)
  output.i = gsub("_output", "_output_tsv", output.i)
  cmd = sprintf("perl _summarize_mhci_results.pl %s %s", input.i, output.i)
  system(cmd)
}

#--------------------------------------------------------------------
# MHC-II
#--------------------------------------------------------------------

outs = list.files("for_netmhciipan/_output", pattern="txt", 
                  full.names="TRUE")

length(outs)
outs[1:3]

for(i in 1:length(outs)){
  print(i)
  input.i  = outs[i]
  output.i = gsub("txt$", "tsv", input.i)
  print(output.i)
  output.i = gsub("_output", "_output_tsv", output.i)
  cmd = sprintf("perl _summarize_mhciipan_results.pl %s %s", input.i, output.i)
  system(cmd)
}

sessionInfo()
q(save="no")
