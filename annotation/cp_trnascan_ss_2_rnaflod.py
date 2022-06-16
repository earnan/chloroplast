#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_trnascan_ss_2_rnaflod.py
#         Author:   yujie
#    Description:   cp_trnascan_ss_2_rnaflod.py
#        Version:   1.0
#           Time:   2022/06/16 10:01:53
#  Last Modified:   2022/06/16 10:01:53
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

from Bio import SeqIO
from Bio.Seq import Seq
#from icecream import ic
import argparse
import linecache
import os
import re
import time
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   cp_trnascan_ss_2_rnaflod.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infasta', metavar='[infasta]', help='fasta', type=str, required=False)
optional.add_argument(
    '-ss', '--infile', metavar='[infile]', help='infile', type=str,  required=False)
optional.add_argument(
    '-n', '--table', metavar='[codon table]', help='默认11', type=int, default=11, required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[outdir]', help='outdir', type=str,  required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()
# ##################################################################################名字映射


def name_mapping(s, table=11):  # 名字映射 2脊椎动物    5无脊椎动物
    if table == 11:
        amino_acid_1 = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                        '*', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', '*', 'W', 'R', 'S', 'R', 'G']
        amino_acid_2 = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                        'Tyr', 'Ter', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Ter', 'Trp', 'Arg', 'Ser', 'Arg', 'Gly']

    if s in amino_acid_2:
        i = amino_acid_2.index(s)
        c = amino_acid_1[i]
    else:
        if s == 'Sup':
            print('Sup may be W')
            c = 'W'
        else:
            c = input(s+': ')
    return c


# ##########################################################################################运行tRNAscan-SE
if args.infasta:
    abs_path = os.path.abspath(args.infasta)
    indir_path = os.path.dirname(abs_path)
    if args.table == 11:
        cmd = "/share/nas1/yuj/software/miniconda3/envs/trnascan/bin/tRNAscan-SE -o {0}/tRNA.out -f {0}/tRNA.ss -m {0}/tRNA.stats {1} -O".format(
            indir_path, args.infasta)  # -O 细胞器
    os.system(cmd)
# ###############################################################################把ss文件存进不同的字典里,需要用到local()函数
"""
createvar = locals()
createvar['list_'+1] = []
"""
# os.path.dirname(path)	              返回文件路径
# os.path.abspath(path)	              返回绝对路径(包含文件名的全路径)

if args.infile:
    infile_path = args.infile
else:
    infile_path = os.path.join(indir_path, 'tRNA.ss')

with open(infile_path, 'r') as fi_handle:
    createvar = locals()
    line_n = 0  # 行号
    multiple = 1  # 倍数
    for line in fi_handle:
        line_n += 1
        while line_n > 6*multiple:
            multiple += 1
        if (line_n % 6) == 1:  # 需要加上括号 避免出bug
            createvar['dict_'+str(multiple)] = {}
            name = line.strip().split('\t')[0].rstrip(
                ')').replace(".", "-").replace(" (", "-")
            createvar['dict_'+str(multiple)]['name'] = name
        elif line.startswith('Type:'):
            s = line.lstrip('Type:').strip().split('\t')[0]
            createvar['dict_'+str(multiple)]['type'] = name_mapping(s)
        elif line.startswith('Seq:'):
            createvar['dict_'+str(multiple)
                      ]['seq'] = line.lstrip('Seq:').strip()
        elif line.startswith('Str:'):  # 改成括号存进去
            createvar['dict_'+str(multiple)
                      ]['str'] = line.lstrip('Str:').strip().replace(">", "(").replace("<", ")")
print(multiple)
# ########################################################################################写入文件
# os.path.exists(path)                  判断一个目录是否存在
# os.makedirs(path)                     多层创建目录
# os.mkdir(path)                        创建目录
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

for i in range(1, multiple+1):
    name = createvar['dict_'+str(i)]['name']
    rnatype = createvar['dict_'+str(i)]['type']
    seq = createvar['dict_'+str(i)]['seq']
    rnastr = createvar['dict_'+str(i)]['str']

    old = name.split('-')[-3]
    new = name.split('-')[-3].split('a')[0]+rnatype
    file_prefix = name.replace(old, new)
    with open(os.path.join(args.outdir, file_prefix+'.fold'), 'w') as fo_handle:
        all_str = '>{}\n{}\n{}\n'.format(file_prefix, seq, rnastr)
        fo_handle.write(all_str)
# ############################################################################################rnaplot画图
s = "cd trna.structure/trn && for i in *.fold;do echo $i;RNAplot -o svg < $i;done && rename _ss.svg .svg *.svg && rm *.fa *.fold && cd ../../ && rm *.ps"
print(s)
