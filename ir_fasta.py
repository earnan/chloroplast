#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   ir_fasta.py
#         Author:   yujie
#    Description:   ir_fasta.py
#        Version:   1.0
#           Time:   2022/04/18 17:22:12
#  Last Modified:   2022/04/18 17:22:12
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

parser = argparse.ArgumentParser(add_help=False, usage='\npython3   将fa序列反向互补')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[xxx.fasta]', help='输入fa文件', type=str, required=False)
optional.add_argument(
    '-l', '--lenth', metavar='[基因序列长度]', type=bool, help="有-i 参数即可输出", default='1', required=False)
optional.add_argument('-gc', '--gccount', help='GC?默认否,运行则-gc ',
                      action='store_true', required=False)

optional.add_argument('-s1', '--seq1',
                      metavar='[ATG→CAT]', help='DNA反向互补', type=str, required=False)
optional.add_argument('-s2', '--seq2',
                      metavar='[CAU→AUG]', help='RNA反向互补', type=str, required=False)
optional.add_argument('-s3', '--seq3',
                      metavar='[CAU→ATG]', help='RNA与DNA间反向互补', type=str, required=False)
optional.add_argument('-s4', '--seq4',
                      metavar='[U→T]', help='不反向不互补仅替换', type=str, required=False)
optional.add_argument('-s5', '--seq5',
                      metavar='[str.upper()]', help='字符串大写', type=str, required=False)
optional.add_argument('-s6', '--seq6',
                      metavar='[str.lower()]', help='字符串小写', type=str, required=False)

optional.add_argument('-o', '--output',
                      metavar='[ir_xxx.fasta]', help='输出反向后的fa文件', type=str, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def ir1(s):  # 反向互补
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


def ir2(s):  # 反向互补
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
    return c


def readfasta(input_file):  # fa将文件读取为字典及列表
    seq_id = ''
    seq_dict = {}
    seqlist = []
    for line in input_file:
        if line.startswith('>'):
            seq_id = line.strip('\n')
            seq_dict[seq_id] = ''
        else:
            seq_dict[seq_id] += line.strip('\n')

    return seq_dict


def judgment_input_type(string):
    if os.path.isfile(string):
        abs_path = os.path.abspath(string)
        abs_dir = os.path.dirname(abs_path)
        input_file = open(string, 'r')
        for line in input_file:
            if len(line.strip('\n')):  # and (not line.startswith('>')):
                seq = line.strip('\n')
    else:  # type(string) == type("a"):
        seq = string
        abs_dir = os.getcwd()
    return seq, abs_dir


def gccount(seq):
    BaseSum = 0  # 碱基总个数初始化
    no_c, no_g, no_a, no_t, no_n = 0, 0, 0, 0, 0  # 各碱基数量
    BaseSum = len(seq)
    no_c = seq.count('C')
    no_g = seq.count('G')
    no_a = seq.count('A')
    no_t = seq.count('T')
    no_n = seq.count('N')
    # print(BaseSum)
    # print('序列个数:'+str(len(Length)))
    # print('最短的序列长度是:'+str(min(Length)))
    # print('最长的序列长度是:'+str(max(Length)))
    print("总碱基数目:", BaseSum)
    print("A百分比", "%.1f" % ((float(no_a*100))/BaseSum), '%')
    print("T百分比", "%.1f" % ((float(no_t*100))/BaseSum), '%')
    print("C百分比", "%.1f" % ((float(no_c*100))/BaseSum), '%')
    print("G百分比", "%.1f" % ((float(no_g*100))/BaseSum), '%')
    print("N百分比", "%.1f" % ((float(no_n*100))/BaseSum), '%')
    print("GC含量是", "%.1f" % ((float(no_g*100+no_c*100))/BaseSum), '%')
    return 0


if args.input:
    input_file = open(args.input, 'r')
    seq_dict = readfasta(input_file)
    seq_id = list(seq_dict.keys())[0]
    seq = list(seq_dict.values())[0]
    if args.gccount:
        gccount(seq)
    if args.lenth:
        print(len(seq))
    ir_seq = ir1(seq)
    # print(len(ir_seq))
    if args.output:
        output_file = open(args.output, 'w')
        output_file.write(str(seq_id)+'\n')
        output_file.write(ir_seq+'\n')
        output_file.close()
    input_file.close()
    print('done')

if args.seq1:
    (seq, abs_dir) = judgment_input_type(args.seq1)
    if args.output:
        file_name = args.output
    else:
        file_name = "s1.fa"
        print(ir1(seq))
    out_path = os.path.join(abs_dir, file_name)  # 生成绝对路径
    print(out_path)
    output_file = open(out_path, 'w')
    # print(os.path.abspath("s1.fa"))
    # path = os.path.abspath(out_file)#获得绝对路径
    # print(os.path.dirname(path))#绝对路径刨去文件名
    output_file.write(ir1(seq)+'\n')
    output_file.close()

if args.seq2:
    print(ir2(args.seq2))
if args.seq3:
    s = args.seq3
    if s.find('T') >= 0:
        c = ir1(s)
        c = c.replace('T', 'U')
    if s.find('U') >= 0:
        c = ir2(s)
        c = c.replace('U', 'T')
    print(c)
if args.seq4:
    print(args.seq4.replace('U', 'T'))
if args.seq5:
    print(args.seq5.upper())
if args.seq6:
    print(args.seq6.lower())
