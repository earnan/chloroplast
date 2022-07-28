#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_add_gene_seq.py
#         Author:   yujie
#    Description:   cp_add_gene_seq.py
#        Version:   1.0
#           Time:   2022/04/19 13:43:58
#  Last Modified:   2022/04/19 13:43:58
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
\npython3   mt_add_gene_seq.py\n\
\n\
1.常规使用\n\
1.1查看密码子 -i -p \n\
\n\
2.递归查找与存储\n\
2.1起始子查找,-m 最大查找次数\n\
2.2存储序列,-sn 基因名\n\
2.3存储蛋白,-sp 基因名\n\
\n\
Path: E:\OneDrive\jshy信息部\Script\chloroplast\annotation\cp_add_gene_seq.py\n\
Version: V1.0'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infasta', metavar='[infasta]', help='输入fasta文件', type=str, default='F:/Epipactis_helleborine_FULLCP.fsa', required=False)
optional.add_argument(
    '-p', '--posstr', metavar='[pos_str]', help="输入位置,形如'124353-124892:-;126001-126552:-'", type=str, default='68847-69098:-;69781-70072:-;71079-71149:-', required=False)


optional.add_argument(
    '-m', '--maxnumber', metavar='[max_number]', help='最大递归查找次数', type=int, default=0, required=False)
optional.add_argument('-trans', '--trans_flag', help='翻译?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-sn', '--nuc_file_name',
                      metavar='[store 2 dna]', help='默认否,值为NULL,存储则输入gene名', type=str,  default='NULL', required=False)
optional.add_argument('-sp', '--pro_file_name',
                      metavar='[store 2 protein]', help='默认否,值为NULL,存储则输入蛋白名', type=str,  default='NULL', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_file(infasta):  # 读取文件
    with open(infasta, 'r') as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.strip('\n')
    return seq


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


def format_pos(pos_str):  # 读取输入的位置为位置列表
    pos_list = []
    content = pos_str.split(';')
    for ele in content:
        if ele.split(':')[-1] == '-':
            tmp = ele.split(':')[0]+':'+'-1'
            pos_list.append(tmp)
        elif ele.split(':')[-1] == '+':
            tmp = ele.split(':')[0]+':'+'1'
            pos_list.append(tmp)
    return pos_list


####################################################################################


def merge_sequence(pos_list, seq):  # 合并获取到的序列,顺便排一下位置顺序

    # 20220728
    """判断是否是trna,返回一个flag"""
    flag_gene_type = 'NULL'
    len_trna_type = 0
    if len(pos_list) == 1:
        start = pos_list[0].split(':')[0].split('-')[0]
        end = pos_list[0].split(':')[0].split('-')[-1]
        len_trna_type = abs(int(end)-int(start))+1
        if 55 <= len_trna_type <= 100:
            # pos_list[0].split(':')[0]   14323-1527
            flag_gene_type = 'trna'

    cds_seq = ""
    if int(pos_list[0].split(':')[-1]) == -1:
        pos_list = pos_list[::-1]

    for ele in pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        start_index = start-1
        end_index = end
        if strand == (-1):
            # ic('minus')
            # seq[start_index:end_index] 角标从start_index到end_index    取的是索引start-1一直到end  取的是start一直到end的碱基
            cds_seq += ir(seq[start_index:end_index])
            # ic(cds_seq)
        elif strand == (1):
            # ic('plus')
            cds_seq += seq[start_index:end_index]
            # ic(cds_seq)
    return cds_seq, pos_list, flag_gene_type, len_trna_type

#######################################################################################################################


def trans2acid(cds_seq):  # 翻译成氨基酸,返回是否正确以及第一个终止子在基因序列上的相对位置
    start_code_table = ['TTG', 'CTG', 'ATT',
                        'ATC', 'ATA', 'ATG', 'GTG']  # 没有rna编辑
    tmp_flag = False
    inter_number = 0
    if len(cds_seq) % 3 == 1:
        print('len(sequence) not a multiple of three! {}=3n+1'.format(len(cds_seq)))
    elif len(cds_seq) % 3 == 2:
        print('len(sequence) not a multiple of three! {}=3n+2'.format(len(cds_seq)))

    coding_dna = Seq(cds_seq)
    acid = coding_dna.translate(table=11)
    print('------------------------------------------------------------')
    print(acid)

    if cds_seq[0:3] not in start_code_table and cds_seq[0:3] != 'ACG':
        print('#####start is wrong!')
    else:
        if cds_seq[0:3] == 'ACG':
            print('-----------------------------------------may be rna edit!')
        if acid.count('*') > 1:
            print('#####interior is wrong!')
            inter_number = acid.find('*')
            print(inter_number)
            print('\n')
        elif acid.count('*') < 1:
            print('#####end is wrong!')
        else:
            if not acid.endswith('*'):
                print('#####interior is wrong!')
                inter_number = acid.find('*')
                print(inter_number)
                print('\n')
            else:
                tmp_flag = True
                print(
                    '------------------------------------------------------------ok')
    return tmp_flag, inter_number, acid

###################################################################################################################
# 如果内部有终止子,则开始尝试返回新的基因位置


def get_new_pos(tmp_pos_list, inter_number):
    # pos_list = []  # 原位置
    # tmp_pos_list = []  # 排序后位置
    # inter_number = 200  # 包括第一个终止子在内的前面所有氨基酸数
    inter_pos = 3*inter_number  # 包括第一个终止子在内的前面所有碱基数

    strand_list = []
    lenth_list = []
    for ele in tmp_pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        lenth = end-start+1
        lenth_list.append(lenth)
        strand_list.append(strand)

    print(lenth_list)
    print(strand_list)
    lenth_sum = 0
    for i in lenth_list:
        lenth_sum += i
    inter = lenth_sum-inter_pos  # 剩余的碱基数
    # ic(inter)

    if inter <= lenth_list[-1]:
        print('lie in [{}]'.format(tmp_pos_list[-1]))
        # 126552-(200-1) 最后一个碱基位置
        new_pos = inter + int(tmp_pos_list[-1].split('-')[0])-3
        print(new_pos)
        print('\n')
    elif inter_pos <= lenth_list[0]+lenth_list[1]:
        new_pos = 124845  # 124892-(600-552-1)最后一个碱基位置
    return 0


#################################################################################################################
# 存储获取到的dna序列或蛋白


def storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq):  # 20220722 新增子函数
    if flag_gene_type == 'trna':  # 20220629   trna 存起来
        print('\nType: tRNA  Len: '+str(len_trna_type)+'\n')
        current_abs_path = os.getcwd()
        if nuc_file_name != 'NULL':
            with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
                f_handle.write(cds_seq+'\n')

    if flag_gene_type == 'NULL':  # 20220722   把 cds 存起来
        #print('\nType: tRNA  Len: '+str(len_trna_type)+'\n')
        current_abs_path = os.getcwd()
        if nuc_file_name != 'NULL':
            with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
                f_handle.write(cds_seq+'\n')
    return 0
#################################################################################################################
# 循环查找


def loop_look(infasta, posstr, trans_flag, loop_count, maxnumber, nuc_file_name, pro_file_name):
    inter_number = False  # 20220728 add  初始值为false

    seq = read_file(infasta)
    pos_list = format_pos(posstr)
    cds_seq, tmp_pos_list, flag_gene_type, len_trna_type = merge_sequence(
        pos_list, seq)
    print('\n'+cds_seq)

    storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq)

    if trans_flag and (flag_gene_type != 'trna'):  # 翻译
        tmp_flag, inter_number, acid = trans2acid(cds_seq)
        current_abs_path = os.getcwd()
        if pro_file_name != 'NULL':
            with open(os.path.join(current_abs_path, pro_file_name+'.acid'), 'w') as f_handle:
                f_handle.write(str(acid)+'\n')

        if tmp_flag == True:
            new_posstr = posstr
            print(new_posstr)
        elif tmp_flag == False:
            n += 1
            print('第{}次查找中'.format(n))
            new_posstr = '124353-124892:-;126001-126552:-'

            if n <= maxnumber:
                loop_look(infasta, new_posstr, trans_flag,
                          loop_count, maxnumber, nuc_file_name, pro_file_name)
            else:
                print('{}次查找未有结果,取消第{}次查找'.format(n-1, n))
    return tmp_pos_list, inter_number


if __name__ == '__main__':
    """
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    begin_time = time.time()
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('Start Time : {}'.format(start_time))
    #################################################################
    """
    loop_count = 0  # 控制递归次数,在loop_look函数外部定义全局变量
    tmp_pos_list, inter_number = loop_look(
        args.infasta, args.posstr, args.trans_flag, loop_count, args.maxnumber, args.nuc_file_name, args.pro_file_name)
    if type(inter_number) == type(1):
        get_new_pos(tmp_pos_list, inter_number)
    """
    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
    """

"""
def trans2acid(codon):  # 翻译成氨基酸
    """"""
    The Bacterial and Plant Plastid Code (11):
    Stnd    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    This    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts      = ---M---------------M------------MMMM---------------M------------
    Base1       = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2       = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3       = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """"""
    genetic_code_number = 11
    acid = ''
    code_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', }
    acid = code_table[codon]
    return acid
"""
