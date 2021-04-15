#!/bin/bash

matchID=$1
HLAI=$2

# matchID=Pt10_on
# HLAI=HLA-A02:01,HLA-A30:02,HLA-B18:01,HLA-B44:02,HLA-C05:01,HLA-C05:01 

~/bin/NetMHC/netMHCpan-4.1/netMHCpan -xls -xlsfile ../output/${matchID}_hlai_mut.txt -f ../data/netmhci/${matchID}_mut.txt -a ${HLAI} -l 9 > ../output/${matchID}_hlai_mut.log
~/bin/NetMHC/netMHCpan-4.1/netMHCpan -xls -xlsfile ../output/${matchID}_hlai_ref.txt -f ../data/netmhci/${matchID}_ref.txt -a ${HLAI} -l 9 > ../output/${matchID}_hlai_ref.log
