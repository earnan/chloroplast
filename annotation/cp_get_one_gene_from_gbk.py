#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_get_one_gene_from_gbk.py
#         Author:   yujie
#    Description:   cp_get_one_gene_from_gbk.py
#        Version:   1.0
#           Time:   2022/12/30 11:24:55
#  Last Modified:   2022/12/30 11:24:55
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################


from Bio import SeqIO
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
python3 cp_get_one_gene_from_gbk.py -i [gbk file] -s [gene name]\n\
')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-i', '--input', help='gbk file path',
                      default="F:\\4923\\2\\ref_adv\\gbk\\Carya_kweichowensis_MH121170.1.gbk", type=str)
required.add_argument('-s', '--search_name',
                      help=' gene name', default='ycf1', type=str)
optional.add_argument(
    '-info', '--info', help='show update log and exit', action='store_true')
optional.add_argument('-h', '--help', action='help',
                      help='show this help message and exit')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t2023/01/04  🎉init(all): 从gbk中取出特定基因序列')
    sys.exit(0)


def ir(s):  # 反向互补
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
    return c


def merge_sequence(ele, complete_seq):  # 合并获取到的序列
    gene_seq = ""
    tmp_list = []  # 位置列表
    for ele1 in ele.location.parts:
        if ele1.strand == (-1):
            # print('minus')
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际起点,从end中取不用+1
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际终点,从start取+1
            gene_seq += ir(complete_seq[ele1.start:ele1.end])
        elif ele1.strand == (1):
            # print('plus')
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际起点,要+1
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际终点,不用+1
            # 切片没问题,索引从start到end-1,也就是对应start+1到end的序列
            gene_seq += complete_seq[ele1.start:ele1.end]
    return tmp_list, gene_seq


gbk_file_path = args.input
search_name = args.search_name
seq_record = SeqIO.read(gbk_file_path, "genbank")
complete_seq = str(seq_record.seq)
dict_gene_seq = {}
for ele in seq_record.features:
    # and ('pseudo' not in ele.qualifiers.keys()):
    if ele.type == "CDS" and 'gene' in ele.qualifiers.keys():
        if ele.qualifiers['gene'][0] == search_name:
            tmp_list, gene_seq = merge_sequence(ele, complete_seq)
            dict_gene_seq[len(gene_seq)] = gene_seq

print(dict_gene_seq[max(dict_gene_seq.keys())].strip())
