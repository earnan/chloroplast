#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   ir_fasta.py
#         Author:   yujie
#    Description:   ir_fasta.py
#        Version:   2.0
#           Time:   2022/04/18 17:22:12
#  Last Modified:   2022/12/07 10:27:03
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
# from Bio import SeqIO
# from Bio.Seq import Seq
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
# import time
# import copy  # 深度拷贝
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()  # 括号里不加东西，程序会自动补出帮助信息
# 自定的必选参数，要加上required=True才是真正的必选项
required = parser.add_argument_group('required arguments')
optional1 = parser.add_argument_group(
    'optional arguments group 1')  # 自定的可选参数，默认状态都是可选
optional2 = parser.add_argument_group('optional arguments group 2')

# help='input fasta file or seq string',
required.add_argument(
    '-i', '--input', metavar='.fasta/seq', type=str, required=True)
# help='output  reverse-complementary fasta file',
optional1.add_argument('-o', '--output', metavar='ir.fasta', type=str)
optional1.add_argument(
    '-l', '--lenth',  help='display sequence length', action='store_false')
optional1.add_argument(
    '-gc', '--gccount', help='display sequence gc_count', action='store_true')

optional2.add_argument('-s1', '--seq1', help='ATG→CAT', action='store_true')
optional2.add_argument('-s2', '--seq2', help='CAU→AUG', action='store_true')
optional2.add_argument('-s3', '--seq3', help='CAU→ATG', action='store_true')
optional2.add_argument('-s4', '--seq4', help='U→T', action='store_true')
optional2.add_argument(
    '-s5', '--seq5', help='str.upper()', action='store_true')
optional2.add_argument(
    '-s6', '--seq6', help='str.lower()', action='store_true')
args = parser.parse_args()


def ir1(s):  # DNA反向互补
    re = s[::-1]  # 字符串反向
    c = ""  # 定义字符串c接收互补序列
    for i in re:
        if i == 'A':
            c = c + 'T'
        elif i == 'G':
            c = c + 'C'
        elif i == 'T':
            c = c + 'A'
        elif i == 'C':
            c = c + 'G'
        elif i == 'N':
            c = c + 'N'
    return c


def ir2(s):  # RNA反向互补
    re = s[::-1]  # 字符串反向
    c = ""  # 定义字符串c接收互补序列
    for i in re:
        if i == 'A':
            c = c + 'U'
        elif i == 'G':
            c = c + 'C'
        elif i == 'U':
            c = c + 'A'
        elif i == 'C':
            c = c + 'G'
        elif i == 'N':
            c = c + 'N'
    return c


def readfasta(fasta_path):  # 将fa文件读取为字典,自带判断是不是fasta格式文件
    seq_dict = {}
    first_content = linecache.getline(fasta_path, 1)
    if first_content.startswith('>'):
        with open(fasta_path, 'r') as fa_handel:
            seq_id = ''
            seq_dict = {}
            for line in fa_handel:
                if line.startswith('>'):
                    seq_id = line.strip()
                    seq_dict[seq_id] = ''
                else:
                    seq_dict[seq_id] += line.strip()
    elif not first_content.startswith('>'):
        seq_dict[os.path.basename(fasta_path)] = first_content.strip()
    return seq_dict


def judgment_input_type(input_str):  # 判断输入的是文件还是粘贴的序列
    seq_dict = {}
    if os.path.isfile(input_str):
        abs_path = os.path.abspath(input_str)
        abs_dir = os.path.dirname(abs_path)
        seq_dict = readfasta(abs_path)
    else:
        seq_dict['undefined'] = input_str.strip()
        abs_dir = os.getcwd()
    return seq_dict, abs_dir


def gccount(seq):
    seq = seq.upper()
    BaseSum = 0  # 碱基总个数初始化
    no_c, no_g, no_a, no_t, no_n = 0, 0, 0, 0, 0  # 各碱基数量
    BaseSum = len(seq)
    no_c = seq.count('C')
    no_g = seq.count('G')
    no_a = seq.count('A')
    no_t = seq.count('T')
    no_n = seq.count('N')
    s = 'Total_base_number:{}bp\n\
        A_percentage:{:.1f}%\n\
            T_percentage:{:.1f}%\n\
                C_percentage:{:.1f}%\n\
                    G_percentage:{:.1f}%\n\
                        N_percentage:{:.1f}%\n\
                            GC_content:{:.1f}%\n'.format(
        BaseSum, no_a*100/BaseSum, no_t*100/BaseSum, no_c*100/BaseSum, no_g*100/BaseSum, no_n*100/BaseSum, (no_g+no_c)*100/BaseSum)
    return s


# 必选参数
seq_dict, abs_dir = judgment_input_type(args.input)
# 可选参数组1
if args.output and (not args.seq1) and (not args.seq2) and (not args.seq3) and (not args.seq4) and (not args.seq5) and (not args.seq6):
    with open(args.output, 'w') as output_handle:
        for seq_id, seq in seq_dict.items():
            ir_seq = ir1(seq)
            output_handle.write(seq_id+'\n')
            output_handle.write(ir_seq+'\n')
    print('done')
for seq_id, seq in seq_dict.items():
    if args.gccount:
        s = gccount(seq)
        print(seq_id)
        print(s)
    if args.lenth:
        print(len(seq), seq_id)

# 可选参数组2
n = 0
for seq_id, seq in seq_dict.items():
    n += 1
    if args.seq1:
        print(ir1(seq))
        if args.output:
            file_name = args.output
        else:
            file_name = seq_id+'s1.fa'
        out_path = os.path.join(abs_dir, file_name)  # 生成绝对路径
        print(out_path)
        with open(out_path, 'w') as output_file:
            # print(os.path.abspath("s1.fa"))
            # path = os.path.abspath(out_file)#获得绝对路径
            # print(os.path.dirname(path))#绝对路径刨去文件名
            output_file.write(ir1(seq)+'\n')
    if args.seq2:
        print(ir2(seq))
    if args.seq3:
        s = seq
        if s.find('T') >= 0:
            c = ir1(s)
            c = c.replace('T', 'U')
        if s.find('U') >= 0:
            c = ir2(s)
            c = c.replace('U', 'T')
        print(c)
    if args.seq4:
        print(seq.replace('U', 'T'))
    if args.seq5:
        print(seq.upper())
    if args.seq6:
        print(seq.lower())
