#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_irscope.py
#         Author:   yujie
#    Description:   cp_irscope.py
#        Version:   1.0
#           Time:   2022/12/29 16:49:16
#  Last Modified:   2023/01/04 16:49:16
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
import copy  # 深度拷贝
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
python3 cp_irscope.py -i [gbk dir] [-s {search_name}]\n\
')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-i', '--input', help='gbk dir path',
                      default='/share/nas1/yuj/project/GP-20220919-4923_20221212/2/ref_adv2', type=str)
#required.add_argument('-o', '--output', help='output dir path', type=str)
optional.add_argument('-s', '--search_name',
                      help='search gene', default='ycf1', type=str)
optional.add_argument(
    '-info', '--info', help='show update log and exit', action='store_true')
optional.add_argument('-h', '--help', action='help',
                      help='show this help message and exit')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t2023/01/04 🎉init(all): irscope分析, 对gbk进行处理\n\
\t\t1.修改gbk中物种名称\n\
\t\t2.根据gbk查找基因\n\
\t\t3.修改mauve中显示的名称')
    sys.exit(0)

# #################################################################################子函数


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


# #########################################################初始化文件路径
ref_adv_dir = args.input

adv_gbk_dir = os.path.join(ref_adv_dir, 'gbk')
gbk_path_list = []
for x in os.listdir(adv_gbk_dir):
    if x.endswith('gbk'):
        gbk_path_list.append(os.path.join(adv_gbk_dir, x))
gbk_path_list.sort()

adv_fa_dir = os.path.join(ref_adv_dir, 'fasta')
fa_path_list = []
for x in os.listdir(adv_fa_dir):
    if x.endswith('fasta'):
        fa_path_list.append(os.path.join(adv_fa_dir, x))


#start_time = time.time()

# ##########################################################################1.修改gbk中显示的名称
for i in gbk_path_list:
    edit_flag1 = False
    edit_flag2 = False
    with open(i, 'r') as gbk_handle:
        tmp_str = ''
        for line in gbk_handle:
            if line.startswith('VERSION'):
                accession = line.strip().split()[1]
                tmp_str += line
            elif line.startswith('  ORGANISM'):  # 修改irscope中的显示
                if line.find(accession) >= 0:
                    edit_flag1 = True
                    new_str = line
                else:
                    new_str = line.rstrip()+' '+accession+'\n'
                tmp_str += new_str
            elif line.startswith('                     /organism="'):  # 修改mauve中的显示
                if line.find(accession) >= 0:
                    edit_flag2 = True
                    new_str = line
                else:
                    new_str = line.rstrip().rstrip('"')+' '+accession+'"'+'\n'
                tmp_str += new_str
                break
            else:
                tmp_str += line
        txt_all = gbk_handle.read()
    if edit_flag1 == False or edit_flag2 == False:
        with open(i, 'w') as gbk_handle:
            gbk_handle.write(tmp_str+txt_all)

# #######################################################################2.查找基因,以ycf1为例
search_name = args.search_name
for i in gbk_path_list:
    unique_fa_path = os.path.join(
        adv_fa_dir, os.path.basename(i).rstrip('.gbk')+'.fasta')
    unique_prefix = os.path.basename(i).rstrip('.gbk')
    seq_record = SeqIO.read(i, "genbank")
    complete_seq = str(seq_record.seq)
    dict_gene_seq = {}
    for ele in seq_record.features:
        # and ('pseudo' not in ele.qualifiers.keys()):
        if ele.type == "CDS" and 'gene' in ele.qualifiers.keys():  # 读取gbk中的信息
            if ele.qualifiers['gene'][0] == search_name:
                tmp_list, gene_seq = merge_sequence(ele, complete_seq)
                dict_gene_seq[len(gene_seq)] = gene_seq
    with open(f'{ref_adv_dir}/{search_name}_{unique_prefix}', 'w') as fo:
        fo.write(dict_gene_seq[max(dict_gene_seq.keys())].strip())
    #cmd1 = f'python3 /share/nas1/yuj/script/chloroplast/annotation/cp_get_one_gene_from_gbk.py -i {i} -s {search_name} > {ref_adv_dir}/{search_name}_{unique_prefix}'
    cmd2 = f'blastn -outfmt 6 -query {ref_adv_dir}/{search_name}_{unique_prefix} -subject {unique_fa_path} > {ref_adv_dir}/blastn_{search_name}_{unique_prefix}'
    # os.system(cmd1)
    os.system(cmd2)


# #######################################################################展现ycf1的假基因信息
#search_name = 'ycf1'
gbk_count = 0
for i in gbk_path_list:  # i: gbk_path
    gbk_count += 1
    unique_prefix = os.path.basename(i).rstrip('.gbk')
    #s = ''
    print('\n', gbk_count, i)
    print(f'{"-"*40}{search_name}{"-"*40}')
    with open(f'{ref_adv_dir}/blastn_{search_name}_{unique_prefix}', 'r') as blastn_handle:
        dict_blastn_stat = {}
        ln = 0
        for line in blastn_handle:
            ln += 1
            # blastn_type_note = f'gene{" "*12}'  # cds或gene 标签
            # pseudo_note = f'\n{" "*21}/pseudo'  # /pseudo 标签
            blastn_lenth = int(line.split()[3])
            blastn_start = int(line.split()[-4])
            blastn_end = int(line.split()[-3])
            if blastn_start < blastn_end:
                query_point_1 = blastn_start
                query_point_2 = blastn_end
                query_pos_info = f'{query_point_1}..{query_point_2}'
            elif blastn_start > blastn_end:
                query_point_1 = blastn_end
                query_point_2 = blastn_start
                query_pos_info = f'complement({query_point_1}..{query_point_2})'

            if blastn_lenth not in dict_blastn_stat.keys():
                dict_blastn_stat[blastn_lenth] = []
                dict_blastn_stat[blastn_lenth].append(query_pos_info)
            else:
                dict_blastn_stat[blastn_lenth].append(query_pos_info)

            if search_name == 'ycf1':
                if 700 <= blastn_lenth <= 1300:  # 假基因的长度
                    blastn_type_note = f'gene{" "*12}'  # cds或gene 标签
                    pseudo_note = f'\n{" "*21}/pseudo'  # /pseudo 标签
                    s = f'{" "*5}{blastn_type_note}{query_pos_info}\n{" "*21}/gene ="{search_name}"{pseudo_note}'
                    print(s)
                elif 4000 <= blastn_lenth:  # 正常基因长度
                    pseudo_note = ''
                    blastn_type_note = f'CDS{" "*13}'
                    s = f'{" "*5}{blastn_type_note}{query_pos_info}\n{" "*21}/gene ="{search_name}"{pseudo_note}'
                    print(s)
                else:
                    s = ''

            elif search_name == 'rps19':
                if 240 <= blastn_lenth <= 320:  # 完整的基因
                    pseudo_note = ''
                    blastn_type_note = f'CDS{" "*13}'
                    s = f'{" "*5}{blastn_type_note}{query_pos_info}\n{" "*21}/gene ="{search_name}"{pseudo_note}'
                    print(s)
                else:
                    s = ''

            else:  # 其他基因
                blastn_type_note = f'gene{" "*12}'  # cds或gene 标签
                pseudo_note = f'\n{" "*21}/pseudo'  # /pseudo 标签
                s = f'{" "*5}{blastn_type_note}{query_pos_info}\n{" "*21}/gene ="{search_name}"{pseudo_note}'
                print(s)

            '''
            subject_point_1  gbk中的起点,不论正负
            subject_point_2  gbk中的终点,不论正负
            query_point_1    待查找的起点,不论正负
            query_point_2    待查找的终点,不论正负
            tmp_point_1      上一组起点,不管正负
            tmp_point_2      上一组终点,不管正负
            '''
            tmp_point_1 = 0
            tmp_point_2 = 0
            if s != '':
                with open(i, 'r') as fi_handle:
                    tmp_str = ''
                    n = 0
                    for line in fi_handle:
                        n += 1
                        if line.startswith(f'{" "*5}gene'):
                            subject_point_1 = int(
                                re.findall(r'[0-9]+', line)[0])
                            subject_point_2 = int(
                                re.findall(r'[0-9]+', line)[-1])
                            if query_point_1 == subject_point_1 or query_point_2 == subject_point_2:  # 已经存在
                                print(
                                    f'{n}{line.replace(" ","!").rstrip()}{"-"*5}')
                                tmp_str += line
                            elif n > 30 and query_point_2 < subject_point_1 and query_point_1 > tmp_point_2:  # 查找插入的位置 要在30行外查找,排除第一个基因影响
                                print(
                                    f'{n}{line.replace(" ","-").rstrip()}{"-"*5}')
                                tmp_str += (s+'\n'+line)
                                break
                            else:
                                tmp_str += line

                            tmp_point_1 = subject_point_1  # 存入上一组起点,不管正负
                            tmp_point_2 = subject_point_2  # 存入上一组终点.不管正负
                        else:
                            tmp_str += line
                    txt_all = fi_handle.read()
                #if flag:tmp_str+txt_all

# source
# gene
# CDS tRNA rRNA
    print(dict_blastn_stat, '\n')
