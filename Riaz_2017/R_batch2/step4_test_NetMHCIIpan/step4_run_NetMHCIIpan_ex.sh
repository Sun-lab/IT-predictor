#!/bin/bash

matchID=Pt89_on
HLAII=DRB1_0403,DRB1_0901,HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10401,HLA-DQA10301-DQB10302,HLA-DQA10301-DQB10303,HLA-DQA10302-DQB10302,HLA-DQA10302-DQB10303


~/bin/NetMHC/netMHCIIpan-4.0/netMHCIIpan -xls -xlsfile ${matchID}_hlaii_mut.txt -context 1 -f ../../data/netmhcii/${matchID}_mut.txt -a ${HLAII} -length 15 > ${matchID}_hlaii_mut.log
