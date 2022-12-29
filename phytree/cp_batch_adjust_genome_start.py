#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_batch_adjust_genome_start.py
#         Author:   yujie
#    Description:   cp_batch_adjust_genome_start.py
#        Version:   1.0
#           Time:   2022/12/19 13:38:47
#  Last Modified:   2022/12/20 13:38:47
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################

#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Blast import NCBIWWW
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
import time
# import copy  # 深度拷贝
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
python3 cp_batch_adjust_genome_start.py -i1 [fasta dir] [-o1 [fa out dir]] -i2 [gbk dir] \n\
python3 cp_batch_adjust_genome_start.py -i1 [fasta dir] [-o1 [fa out dir]]\n\
python3 cp_batch_adjust_genome_start.py -i2 [gbk dir] \n\
Take lsc as the starting point to adjust the chloroplast genome.')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-i1', '--input1', help='fasta dir path', type=str)
optional.add_argument('-o1', '--output1', help='fa out path', type=str)
required.add_argument('-i2', '--input2', help='gbk dir path', type=str)
#optional.add_argument('-o2', '--output2', help='gbk out path', type=str)
optional.add_argument(
    '-info', '--info', help='show update log and exit', action='store_true')
optional.add_argument('-h', '--help', action='help',
                      help='show this help message and exit')
args = parser.parse_args()

if args.info:
    print('\nUpdate log:')
    print('\t2022/12/19 🎉init(all): 以LSC为起点,调整叶绿体基因组序列')
    print('\t2022/12/29 ✨feat(main): 生成对gbk文件的调整命令')
    sys.exit(0)

# fasta
infasta_dir = os.path.abspath(args.input1)
if not args.output1:
    outfasta_dir = os.path.join(os.path.dirname(infasta_dir), 'complete2')
else:
    outfasta_dir = args.output1
if not os.path.exists(outfasta_dir):
    os.makedirs(outfasta_dir)

# gbk
ingbk_dir = os.path.abspath(args.input2)


if 'ir.log' not in os.listdir(infasta_dir):
    print('log file does not exist')
else:
    log_path = os.path.join(infasta_dir, 'ir.log')
    with open(log_path, 'r') as log_handle:
        dict_seq_start = {}
        dict_gbk_start = {}
        for line in log_handle:
            if line.startswith('-'*6):
                seq_id = line.split(' ')[1]
                dict_seq_start[seq_id] = ''
                dict_gbk_start[seq_id] = ''
            else:
                if line.startswith('LSC:'):
                    start = line.strip().lstrip(
                        'LSC:').split('-')[0]
                    # fa
                    cmd1 = 'python3 /share/nas1/yuj/script/mitochondrion/annotation/mt_move_gene_pos.py -fa "{}" -s {} -o "{}" &'.format(
                        os.path.join(args.input1, seq_id), start, os.path.join(outfasta_dir, seq_id))
                    dict_seq_start[seq_id] = cmd1
                    # gbk
                    if args.input2:
                        cmd2 = 'perl /share/nas6/xul/program/chloroplast/bin/cp_format_gbk.pl -f "{}.gbk" -s {} &'.format(
                            os.path.join(args.input2, seq_id).rstrip('.fasta'), start)
                        dict_gbk_start[seq_id] = cmd2


commands1_sh_path = os.path.join(infasta_dir, 'commands1.sh')
commands2_sh_path = os.path.join(ingbk_dir, 'commands2.sh')
with open(commands1_sh_path, 'w') as cmd_handle:
    for k, v in dict_seq_start.items():
        cmd_handle.write(v+'\n')
with open(commands2_sh_path, 'w') as cmd_handle:
    for k, v in dict_gbk_start.items():
        cmd_handle.write(v+'\n')

os.system('sh {}'.format(commands1_sh_path))
os.system('sh {}'.format(commands2_sh_path))
