#!/bin/bash

matchID=$1
HLAII=$2

# matchID=Pt10_on
# HLAII=DRB1_0101,DRB1_0301,HLA-DQA10101-DQB10501,HLA-DQA10501-DQB10201

~/bin/NetMHC/netMHCIIpan-4.0/netMHCIIpan -xls -xlsfile ../output/${matchID}_hlaii_mut.txt -context 1 -f ../data/netmhcii/${matchID}_mut.txt -a ${HLAII} -length 15 > ../output/${matchID}_hlaii_mut.log
~/bin/NetMHC/netMHCIIpan-4.0/netMHCIIpan -xls -xlsfile ../output/${matchID}_hlaii_ref.txt -context 1 -f ../data/netmhcii/${matchID}_ref.txt -a ${HLAII} -length 15 > ../output/${matchID}_hlaii_ref.log
