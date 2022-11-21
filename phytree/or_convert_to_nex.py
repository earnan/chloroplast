#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   or_convert_to_nex.py
#         Author:   yujie
#    Description:   or_convert_to_nex.py
#        Version:   1.0
#           Time:   2022/10/14 14:25:00
#  Last Modified:   2022/10/14 14:25:00
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
import time
import copy  # 深度拷贝
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   or_convert_to_nex.py\n\
#         Author:   yujie\n\
#    Description:   or_convert_to_nex.py\n\
#        Version:   1.0\n\
#           Time:   2022/10/14 14:38:06\n\
#  Last Modified:   2022/10/14 14:38:06\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
#\n\
##########################################################\n\
\n\
\npython3   or_convert_to_nex.py\n\
根据比对好的fasta文件和模型预测结果生成.nex格式\n\
Function:\n\
1.常规使用\n\
1.1 -i []  -m [] -o [] \n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\chloroplast\phytree\or_convert_to_nex.py\n\
Path: /share/nas1/yuj/script/chloroplast/phytree/or_convert_to_nex.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--inaln', metavar='[aln.fa]', help='inaln', type=str, default='E:/Server/3463-3/fullbi/3463-3fullml.trim.aln.fasta', required=False)
optional.add_argument(
    '-m', '--modeltest', metavar='[jmodeltest file]', help='jmodeltest file', type=str, default='E:/Server/3463-3/fullbi/jmodeltest', required=False)
optional.add_argument(
    '-o', '--outnex', metavar='[.nex]', help='.nex', type=str, default='F:/3463.nex', required=False)
#optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',action='store_false', required=False)
#optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',action='store_true', required=False)
optional.add_argument('-info', help='更新日志,使用时-info',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[help information]')
args = parser.parse_args()


if args.info:
    print('\n更新日志:')
    print('\t20221117  模型预测结果解析完成 第一版ok')
    print('\n')
    sys.exit(0)


#################################################################
# 格式化成2016-03-20 11: 45: 39形式
begin_time = time.time()
start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
print('Start Time : {}'.format(start_time))
#################################################################

# 文件路径
in_aln_fa_path = args.inaln  # i
in_model_path = args.modeltest  # m
out_nex_path = args.outnex  # o


def format_fasta(seq, num):  # 格式化字符串
    number_of_effective_rows = 1
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        if (index + 1) % num == 0:  # 可以用来换行
            format_seq += "\n"
            number_of_effective_rows += 1
    return format_seq, number_of_effective_rows

# 读取比对文件,返回formatted_accession_seq_dict


def get_formatted_seq_dict(in_aln_fa_path):  # 贼耗时
    with open(in_aln_fa_path, 'r') as fa_handle:
        # nex存放具体序列的部分,在开头
        all_seq_dict = {}  # 直接读取文件  :  序列,无换行
        formatted_accession_seq_dict = {}  # 登录号 :  序列,换行
        all_seq_id_dict = {}  # id : 名字长度
        previous_current_note = ''  # 上一个id
        previous_final_note = ''  # 上一个id
        all_seq_dict[previous_current_note] = ''  # 赋空值,第一次要写的序列为空
        format_seq = ''
        int_ntax = 0  # 序列总个数

        for line in fa_handle:
            if line.startswith('>'):
                '''写入序列'''
                int_ntax += 1
                format_seq, number_of_effective_rows = format_fasta(
                    all_seq_dict[previous_current_note], 60)
                if format_seq != '':  # 由于是先写序列,后写id,直接写入,那么第一行会空出一行
                    # 第一次要写的序列为空,正好可以判断为空就不写入
                    formatted_accession_seq_dict[previous_final_note] = format_seq
                '''读取成字典'''
                current_note = line.strip().lstrip('>').replace('.', '_')  # 改为下划线的登录号
                all_seq_id_dict[current_note] = len(current_note)
                previous_current_note = current_note
                all_seq_dict[current_note] = ''
                previous_final_note = current_note
                formatted_accession_seq_dict[previous_final_note] = ''
            else:
                all_seq_dict[current_note] += line.strip().upper()

        str_seq_len = len(all_seq_dict[previous_current_note])
        #  最后一个序列,要单独拿出来存入字典
        format_seq, number_of_effective_rows = format_fasta(
            all_seq_dict[previous_current_note], 60)
        formatted_accession_seq_dict[previous_final_note] = format_seq
    return all_seq_id_dict, formatted_accession_seq_dict, int_ntax, str_seq_len, number_of_effective_rows


def get_model_parameters(in_model_path):
    dict_nst = {'JC': '1', 'F81': '1', 'K80': '2', 'HKY': '2', 'GTR': '6'}
    dict_rates = {'G': 'gamma', 'I': 'propinv', 'I+G': 'invgamma'}
    str_nst, str_rates, str_ngammacat, str_statefreqpr, str_revmat, str_pinvar, str_shapepr = '', '', '5', '', '', '', ''
    bic_number_of_row = 0
    bic = ''
    with open(in_model_path, 'r') as model_handle:
        number_of_row = 0
        for line in model_handle:
            number_of_row += 1
            if line.startswith('*             BAYESIAN INFORMATION CRITERION (BIC) '):
                bic_number_of_row = number_of_row
            # 22代表bic这一部分占了22行，可能不同文件不一样
            if bic_number_of_row != 0 and (bic_number_of_row <= number_of_row <= bic_number_of_row+22):
                bic += line
                # model
                if line.startswith('   Model'):
                    str_model = line.split()[-1]
                    str_nst = dict_nst[str_model.split('+')[0]]
                    if len(str_model.split('+')) == 1:
                        flag_str_rates = input(
                            "please enter 1 or 2\n>1:equal\n>2:adgamma\n: ")
                        if flag_str_rates == '1':
                            str_rates = 'equal'
                        elif flag_str_rates == '2':
                            str_rates = 'adgamma'
                    elif len(str_model.split('+')) == 2:
                        str_rates = dict_rates[str_model.split('+')[-1]]
                    elif len(str_model.split('+')) == 3:
                        str_rates = dict_rates['+'.join(
                            str_model.split('+')[1:])]
                # statefreqpr
                elif line.strip() != '' and line.split()[0] in ['freqA', 'freqC', 'freqG', 'freqT']:
                    str_statefreqpr += (line.split()[-1]+',')
                # str_revmat
                elif line.strip() != '' and line.split()[0] in ['R(a)', 'R(b)', 'R(c)', 'R(d)', 'R(e)', 'R(f)']:
                    str_revmat += (line.split()[-1]+',')
                # str_pinvar
                elif line.startswith('   p-inv'):
                    str_pinvar = line.split()[-1]
                # str_shapepr
                elif line.startswith('   gamma shape'):
                    str_shapepr = line.split()[-1]
            if bic_number_of_row != 0 and (number_of_row > bic_number_of_row+22):
                break
    bic = "---------------------------------------------------------------\
         \n*                                                             *\n"+bic
    str_statefreqpr = str_statefreqpr.rstrip(',')
    str_revmat = str_revmat.rstrip(',')
    return bic, str_nst, str_rates, str_ngammacat, str_statefreqpr, str_revmat, str_pinvar, str_shapepr


bic, str_nst, str_rates, str_ngammacat, str_statefreqpr, str_revmat, str_pinvar, str_shapepr = get_model_parameters(
    in_model_path)
all_seq_id_dict, formatted_accession_seq_dict, int_ntax, str_seq_len, number_of_effective_rows = get_formatted_seq_dict(
    in_aln_fa_path)
all_seq_id_dict = copy.deepcopy(dict(sorted(all_seq_id_dict.items(
), key=lambda x: x[0], reverse=False)))  # 字典按照键的大小排序  再转回字典深拷贝   进而完成对初始字典的排序

print(bic)
print(str_nst, str_rates, str_ngammacat, str_statefreqpr,
      str_revmat, str_pinvar, str_shapepr)
print(int_ntax, str_seq_len)


with open(out_nex_path, 'w') as nex_handle:

    s = "#NEXUS\n\
BEGIN DATA;\n\
dimensions ntax={0} nchar={1};\n\
format missing=?\n\
datatype=DNA gap= - interleave;\n\
\n\
matrix\n\
".format(int_ntax, str_seq_len)
    nex_handle.write(s)

    for i in range(number_of_effective_rows):
        for j in all_seq_id_dict.keys():
            nex_handle.write(j.ljust(max(all_seq_id_dict.values())+1) +
                             formatted_accession_seq_dict[j][61*i:61*(i+1)-1]+'\n')
        if i < number_of_effective_rows-1:
            nex_handle.write('\n')

    s = ";\n\
END;\n\
begin mrbayes;\n\
log start filename = log.txt;\n\
lset nst={} rates={} Ngammacat={};\n\
prset statefreqpr = fixed({}) revmat = fixed({}) pinvar=fixed({}) shapepr = fixed({});\n\
mcmcp ngen=2000000 printfreq=1000 samplefreq=100 nchains=4 nruns=2 savebrlens=yes checkpoint=yes checkfreq=5000;\n\
mcmc;\n\
sumt conformat=Simple contype=Halfcompat relburnin=yes burninfrac=0.25;\n\
sump relburnin=yes burninfrac=0.25;\n\
end;\
".format(str_nst, str_rates, str_ngammacat, str_statefreqpr,
         str_revmat, str_pinvar, str_shapepr)
    nex_handle.write(s)


###############################################################
end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
print('Done')
###############################################################

'''
#NEXUS\n\
BEGIN DATA;\n\
dimensions ntax={$int_ntax} nchar={$str_seq_len};\n\
format missing=?\n\
datatype=DNA gap= - interleave;\n\
\n\
matrix\n\
{$formatted_accession_seq_dict}\n\
;\n\
END;\n\
begin mrbayes;\n\
log start filename = log.txt;\n\
lset nst={6} rates={invgamma} Ngammacat={5};\n\
prset statefreqpr = fixed({0.3646,0.1334,0.1327,0.3692}) revmat = fixed({0.969,1.634,0.578,0.452,1.581,1.000}) pinvar=fixed({0.12}) shapepr = fixed({0.61});\n\
mcmcp ngen=2000000 printfreq=1000 samplefreq=100 nchains=4 nruns=2 savebrlens=yes checkpoint=yes checkfreq=5000;\n\
mcmc;\n\
sumt conformat=Simple contype=Halfcompat relburnin=yes burninfrac=0.25;\n\
sump relburnin=yes burninfrac=0.25;\n\
end;\n\
'''
