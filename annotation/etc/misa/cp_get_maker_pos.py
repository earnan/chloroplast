#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_get_maker_pos.py
#         Author:   yujie
#    Description:   cp_get_maker_pos.py
#        Version:   1.0
#           Time:   2022/05/09 16:20:32
#  Last Modified:   2022/05/09 16:20:32
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
\npython3   .py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-1', '--indir', metavar='[dir]', help='all fa',
                      type=str, default='F:/4286', required=False)
optional.add_argument('-2', '--infile', metavar='[.maker]', help='maker',
                      type=str, default='F:/4286/ssr20220509.xls', required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[outdir]', help='outdir', type=str, default='F:/4286', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_fasta_to_dic1(infasta1, infasta2, infasta3):  # 最简单,针对普通fasta文件 >物种名
    if os.path.isfile(infasta1):
        infasta = infasta1
    elif os.path.isfile(infasta2):
        infasta = infasta2
    elif os.path.isfile(infasta3):
        infasta = infasta3
    with open(infasta, 'r') as f:
        seq_id = ''
        id_index = []
        dict_seq = {}
        for line in f:
            # 如果是">"开头的，就创建一个key键
            if line.startswith('>'):
                seq_id = line.strip('\n').lstrip('>')  # ID为键
                id_index.append(line.replace(
                    "\n", "").replace(">", ""))  # 顺便创建索引的列表
                dict_seq[seq_id] = ''  # 有key无value
            # 如果不是">"开头的，在这个key键下添加value值
            else:
                dict_seq[seq_id] += line.strip('\n')
    return dict_seq


def get_maker_in_fa_pos(indir, infile):
    # ###############################获取所有分析物种
    list_species = []
    with open(infile, 'r') as f:
        line1 = f.readline().strip()
        content = line1.split('\t')
        for i in range(2, len(content)):
            list_species.append(content[i])
    # ################################取出每个物种的ssr序列
    # 初始化变量
    createvar = locals()
    for i in range(len(list_species)):
        createvar['list_ssr_'+str(i+1)] = []
    # 使用pandas库
    test_content = pd.read_table(infile)  # 有标题,根据标题来取值
    for i in range(len(list_species)):
        createvar['list_SSR_'+str(i+1)] = test_content[list_species[i]]
        for j in createvar['list_SSR_'+str(i+1)]:
            createvar['list_ssr_'+str(i+1)].append(j)
    # ###################################打开对应物种的.fasta完整序列
    # 初始化变量
    for i in range(len(list_species)):
        infasta1 = indir+'/'+list_species[i]+'.fasta'
        infasta2 = indir+'/'+list_species[i]+'.1.fasta'
        infasta3 = indir+'/'+list_species[i].split('_')[0]+'.fasta'
        dict_seq = read_fasta_to_dic1(infasta1, infasta2, infasta3)
        createvar['seq'+str(i+1)] = dict_seq[list_species[i]]
    """
    return createvar['list_ssr_1'], createvar['list_ssr_2'], createvar['list_ssr_3'], createvar['list_ssr_4'], createvar['list_ssr_5'], \
        createvar['list_ssr_6'], createvar['list_ssr_7'], createvar['list_ssr_8'], createvar['list_ssr_9'], createvar['list_ssr_10'], createvar['list_ssr_11'], \
        createvar['seq1'], createvar['seq2'], createvar['seq3'], createvar['seq4'], createvar['seq5'],\
        createvar['seq6'], createvar['seq7'], createvar['seq8'], createvar['seq9'], createvar['seq10'], createvar['seq11'], list_species
createvar = locals()
createvar['list_ssr_1'], createvar['list_ssr_2'], createvar['list_ssr_3'], createvar['list_ssr_4'], createvar['list_ssr_5'], \
    createvar['list_ssr_6'], createvar['list_ssr_7'], createvar['list_ssr_8'], createvar['list_ssr_9'], createvar['list_ssr_10'], createvar['list_ssr_11'], \
    createvar['seq1'], createvar['seq2'], createvar['seq3'], createvar['seq4'], createvar['seq5'],\
    createvar['seq6'], createvar['seq7'], createvar['seq8'], createvar['seq9'], createvar['seq10'], createvar['seq11'], list_species = get_maker_in_fa_pos(
    args.indir, args.infile)
    # list_ssr_1, list_ssr_2, list_ssr_3, list_ssr_4, list_ssr_5, list_ssr_6, list_ssr_7, list_ssr_8, list_ssr_9,  list_ssr_10,  list_ssr_11,\
    # seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9, seq10, seq11 = get_maker_in_fa_pos(
    # args.indir, args.infile)
    """
    # ##########################获取位置
    dict_all_pos_list = {}
    for i in range(len(list_species)):
        tmp_str = ''
        seq = createvar['seq'+str(i+1)]
        list_ssr = createvar['list_ssr_'+str(i+1)]
        tmp_str = list_species[i]+':' + seq
        dict_all_pos_list[tmp_str] = []
        lenth = 0
        start = 0
        for j in range(len(list_ssr)):
            ssr = createvar['list_ssr_'+str(i+1)][j]
            lenth = len(ssr)
            n = seq.find(ssr, start)
            dict_all_pos_list[tmp_str].append([n+1, lenth])
            start = n+lenth
        # print(len(dict_all_pos_list[list_species[i]]))
    return dict_all_pos_list, list_species


"""
PRIMER_SEQUENCE_ID=Camellia_sinensis_L_O_Kuntze_cv_Xillian_1_1
SEQUENCE=GGGCGAACGACGGGAATTG
TARGET=353,18
PRIMER_PRODUCT_SIZE_RANGE=100-280
PRIMER_MAX_END_STABILITY=250
=
"""
"""start-3 size+6"""


def write_p3in(dict_all_pos_list, list_species, outdir):
    fa_dict = {}
    fa_list = []
    for i in dict_all_pos_list.keys():
        fa_dict[i.split(':')[0]] = i.split(':')[1]
        fa_list.append(i)

    for i in range(len(list_species)):
        outfile = outdir+'/'+list_species[i]+'.ssr.p3in'
        ic(outfile)
        with open(outfile, 'wb') as f:
            n = 1
            for j in dict_all_pos_list[fa_list[i]]:
                f.write('PRIMER_SEQUENCE_ID={0}_{1}\n'.format(
                    list_species[i], n).encode())
                f.write('SEQUENCE={}\n'.format(
                    fa_dict[list_species[i]]).encode())
                f.write('TARGET={0},{1}\n'.format(j[0]-3, j[1]+6).encode())
                f.write('PRIMER_PRODUCT_SIZE_RANGE=100-280\n'.encode())
                f.write('PRIMER_MAX_END_STABILITY=250\n'.encode())
                f.write('=\n'.encode())
                n += 1

            print('done')
        #createvar['seq'+str(i+1)] = dict_seq[list_species[i]]
    return fa_list


dict_all_pos_list, list_species = get_maker_in_fa_pos(args.indir, args.infile)
fa_list = write_p3in(dict_all_pos_list, list_species, args.outdir)
for i in dict_all_pos_list[fa_list[0]]:
    print(i[0])
