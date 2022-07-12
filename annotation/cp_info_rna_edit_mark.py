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
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
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
