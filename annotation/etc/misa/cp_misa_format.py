#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_misa_format.py
#         Author:   yujie
#    Description:   cp_misa_format.py
#        Version:   1.0
#           Time:   2022/05/09 10:44:49
#  Last Modified:   2022/05/09 10:44:49
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
import pandas as pd

parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   cp_misa_format.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='F:/4286/Camellia_sinensis_L_O_Kuntze_cv_Xillian_1.misa', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_file(infile):
    list_ssr = []
    list_size = []
    test_content = pd.read_table(infile)  # 有标题,根据标题来取值
    list_SSR = test_content["SSR"]
    list_SIZE = test_content["size"]
    for i in list_SSR:
        list_ssr.append(i)
    for i in list_SIZE:
        list_size.append(i)
    return list_ssr, list_size


list_ssr, list_size = read_file(args.infile)
print(list_ssr)
