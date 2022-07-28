#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_info_rna_edit_mark.py
#         Author:   yujie
#    Description:   cp_info_rna_edit_mark.py
#        Version:   1.0
#           Time:   2022/07/12 16:16:51
#  Last Modified:   2022/07/12 16:16:51
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import re
import time
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   cp_info_rna_edit_mark.py\n\
\npython3   RNA编辑的基因末尾加标识符\n\
\n\
0.需要-i1 -i2 -o\n\
\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i1', '--infile1', metavar='[infile]', help='final_gene_annotaion.info1', type=str, default='final_gene_annotaion.info1', required=False)
optional.add_argument(
    '-i2', '--infile2', metavar='[infile]', help='rna编辑基因列表,一行一个', type=str, default='tmp_gene', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='final_gene_annotaion.info', type=str, default='final_gene_annotaion.info', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

with open('F:\\4326\\final_gene_annotaion.info', 'r') as info_handle, open('F:\\4326\\tmp_gene', 'r') as name_handle, open('F:\\4326\\out.info', 'wb') as out_handle:
    name_list = []
    exist_list = []
    [name_list.append(i.strip()) for i in name_handle]
    n, m = 0, 0
    for i in info_handle:
        n += 1
        gene_name = i.split()[2]
        if gene_name not in name_list:
            out_handle.write(i.encode())
        else:
            m += 1
            i = i.rstrip()+'\t1\n'
            out_handle.write(i.encode())
            exist_list.append(gene_name)
