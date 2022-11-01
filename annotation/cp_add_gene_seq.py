#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_add_gene_seq.py
#         Author:   yujie
#    Description:   cp_add_gene_seq.py
#        Version:   2.0
#           Time:   2022/04/19 13:43:58
#  Last Modified:   2022/11/01 13:43:58
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
# import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
# import time
# import copy  # 深度拷贝
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   cp_add_gene_seq.py\n\
#         Author:   yujie\n\
#    Description:   cp_add_gene_seq.py\n\
#        Version:   2.0\n\
#           Time:   2022/04/19 13:43:58\n\
#  Last Modified:   2022/11/01 13:43:58\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
#\n\
##########################################################\n\
\n\
\npython3   cp_add_gene_seq.py\n\
\n\
1.常规使用\n\
1.1查看密码子 -i -p \n\
1.2序列被误判为trna,强制翻译使用 -f \n\
\n\
2.递归查找与存储\n\
2.1起始子/终止子自动查找,-m 最大查找次数\n\
\n\
3.存储为文件\n\
3.1存储序列,-sn 基因名\n\
3.2存储蛋白,-sp 基因名\n\
\n\
Path: E:\OneDrive\jshy信息部\Script\chloroplast\annotation\cp_add_gene_seq.py\n\
Path: /share/nas1/yuj/script/chloroplast/annotation/cp_add_gene_seq.py\n\
Version: V2.0'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infasta', metavar='[infasta]', help='输入fasta文件', type=str, default='F:/Epipactis_helleborine_FULLCP.fsa', required=False)
optional.add_argument(
    '-p', '--posstr', metavar='[pos_str]', help="输入位置,形如'124353-124892:-;126001-126552:-'", type=str, default='68847-69098:-;69781-70072:-;71079-71149:-', required=False)


optional.add_argument(
    '-m', '--maxnumber', metavar='[max_number]', help='最大递归查找次数,默认0,假查找', type=int, default=0, required=False)
optional.add_argument('-trans', '--trans_flag',
                      help='翻译?默认是,不运行则-trans', action='store_false', required=False)
optional.add_argument('-f', '--force_flag',
                      help='强制翻译?默认否,运行则-f', action='store_true', required=False)
optional.add_argument('-sn', '--nuc_file_name',
                      metavar='[store 2 dna]', help='默认否,值为NULL,存储则输入gene名', type=str,  default='NULL', required=False)
optional.add_argument('-sp', '--pro_file_name',
                      metavar='[store 2 protein]', help='默认否,值为NULL,存储则输入蛋白名', type=str,  default='NULL', required=False)
optional.add_argument('-info', help='更新日志,使用时-info',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t20221101  添加终止子错误时的查找 更新一些提示信息')

    print('\n')
    sys.exit(0)

##########################################################################################
# 子函数 功能简单


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
    # pos_list 某基因格式化的位置
    # seq 全长序列
    """
    20220728
    判断是否是trna,返回一个flag"""
    flag_gene_type = 'NULL'
    len_trna_type = 0
    if len(pos_list) == 1:
        start = pos_list[0].split(':')[0].split('-')[0]
        end = pos_list[0].split(':')[0].split('-')[-1]
        len_trna_type = abs(int(end)-int(start))+1
        if 55 <= len_trna_type <= 100:
            flag_gene_type = 'trna'

    # -----------------20221101 参考线粒体  解决跨首尾基因
    seq_len = len(seq)
    if int(pos_list[0].split(':')[-1]) == 1 and int(pos_list[0].split(':')[0].split('-')[0]) > int(pos_list[0].split(':')[0].split('-')[-1]):  # 14323-1527:1
        pos1 = '{0}-{1}:1'.format(pos_list[0].split(':')
                                  [0].split('-')[0], seq_len)
        pos2 = '1-{}:1'.format(pos_list[0].split(':')[0].split('-')[-1])
        pos_list = [pos1, pos2]
    # ------------------
    cds_seq = ""
    tmp_pos_list = pos_list
    if int(pos_list[0].split(':')[-1]) == -1:
        tmp_pos_list = pos_list[::-1]

    for ele in tmp_pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        start_index = start-1
        end_index = end
        if strand == (-1):
            # seq[start_index:end_index] 角标从start_index到end_index    取的是索引start-1一直到end  取的是start一直到end的碱基
            cds_seq += ir(seq[start_index:end_index])
        elif strand == (1):
            cds_seq += seq[start_index:end_index]
    return cds_seq, tmp_pos_list, flag_gene_type, len_trna_type

#######################################################################################################################


def trans2acid(cds_seq):  # 翻译成氨基酸,返回是否正确以及第一个终止子在基因序列上的相对位置

    start_code_list = ['TTG', 'CTG', 'ATT',
                       'ATC', 'ATA', 'ATG', 'GTG']  # 没有rna编辑

    #end_codon_list = ['TAA', 'TAG', 'TGA']
#
#
#
#
#
#
    tmp_flag = 0  # 0为正确,1为起始x,2为内部错,3为末尾错
    inter_number = 0
    if len(cds_seq) % 3 == 1:
        print('len(sequence) not a multiple of three! {}=3n+1'.format(len(cds_seq)))
    elif len(cds_seq) % 3 == 2:
        print('len(sequence) not a multiple of three! {}=3n+2'.format(len(cds_seq)))

    coding_dna = Seq(cds_seq)
    acid = coding_dna.translate(table=11)
    print('------------------------------------------------------------')
    print(acid)
#
#
#
#
#
    if cds_seq[0:3] == 'ACG':
        print('-----------------------------------------may be rna edit!')

    if cds_seq[0:3] not in start_code_list and cds_seq[0:3] != 'ACG':
        print('#####start is wrong!')
        tmp_flag = 1
    else:  # 起始正确,又分3种情况
        if acid.count('*') > 1:  # 终止多于1,意味着提前终止
            print('#####interior is wrong!')
            tmp_flag = 2
            inter_number = acid.find('*')
            print('Index of the first stop codon :{}'.format(inter_number))
            print('\n')
        elif acid.count('*') < 1:  # 终止小于1 未终止
            print('#####end is wrong!')
            tmp_flag = 3
        else:  # 1个终止子,提前终止或者没问题
            if not acid.endswith('*'):
                print('#####interior is wrong!')
                tmp_flag = 2
                inter_number = acid.find('*')
                print('Index of the first stop codon :{}'.format(inter_number))
                print('\n')
            else:
                tmp_flag = 0
                print(
                    '------------------------------------------------------------ok')
    return tmp_flag, inter_number, acid


###################################################################################################################


# 如果内部有终止子,则开始尝试返回新的基因位置，指开头到第一个终止子这一段


def get_current_first_end_pos(tmp_pos_list, inter_number):
    # tmp_pos_list = []  # 排序后位置,起始子序列在列表第一位，终止子在列表最后一位
    # inter_number = 200  # 包括第一个终止子在内的前面所有密码子个数
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

    print('Lenth list:{0}\tStrand list:{1}\tTotal length:{2}bp'.format(
        lenth_list, strand_list, sum(i for i in lenth_list)))  # 序列可能是多段

    if inter_pos == 0:  # 说明序列正确，终止子在末尾
        print('Stop codon lie in [{}]'.format(
            tmp_pos_list[-1]))  # 说明终止子位于最后一段内
    else:
        remaining_bp = sum(i for i in lenth_list)-inter_pos  # 剩余的碱基数,不包括第一个终止子
        if remaining_bp <= lenth_list[-1]:  # 剩余的长度小于最后一段，说明终止子位于最后一段内
            print('Stop codon lie in [{}]'.format(tmp_pos_list[-1]))
            # 20220808 更正第一个终止密码子出现位置的计算公式
            if tmp_pos_list[-1].split(':')[-1] == 1:  # plus链基因
                current_first_end_pos = int(re.findall(
                    r'\d+', tmp_pos_list[-1])[-1])-remaining_bp-1  # 终止子中间那个碱基位置
                print('{}-{}:+'.format(current_first_end_pos -
                      1, current_first_end_pos+1))
                print('\n')
            elif tmp_pos_list[-1].split(':')[-1] == -1:  # minus链基因
                current_first_end_pos = int(re.findall(
                    r'\d+', tmp_pos_list[-1])[0])+remaining_bp+1
                print('{}-{}:-'.format(current_first_end_pos -
                      1, current_first_end_pos+1))
                print('\n')

        elif inter_pos <= lenth_list[0]:  # 包括第一个终止在内的长度小于第一段，说明终止子位于第一段
            print('Stop codon lie in [{}]'.format(tmp_pos_list[0]))
        else:
            print('Stop codon lie in [{}]'.format(tmp_pos_list[1]))
    #print('Already Run {}s'.format(time.time()-begin_time))
    return 0

#################################################################################################################

# 20220722 新增子函数


def storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq):  # 存储获取到的dna序列或蛋白
    if flag_gene_type == 'trna':  # 20220629   trna 存起来
        print('\nType: tRNA  Len: '+str(len_trna_type)+'\n')
    # 20221101 和线粒体一样 精简
    current_abs_path = os.getcwd()
    if nuc_file_name != 'NULL':
        with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
            f_handle.write(cds_seq+'\n')
    '''
    if flag_gene_type == 'NULL':  # 20220722   把 cds 存起来
        #print('\nType: tRNA  Len: '+str(len_trna_type)+'\n')
        current_abs_path = os.getcwd()
        if nuc_file_name != 'NULL':
            with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
                f_handle.write(cds_seq+'\n')'''
    return 0
#################################################################################################################


# 循环查找
def loop_look(infasta, posstr, trans_flag, loop_count, maxnumber, nuc_file_name, pro_file_name):

    start_codon_list = ['TTG', 'CTG', 'ATT',
                        'ATC', 'ATA', 'ATG', 'GTG']  # 没有rna编辑

    end_codon_list = ['TAA', 'TAG', 'TGA']

    inter_number = False  # 20220728 add  初始值为false
#
#
    seq = read_file(infasta)
    pos_list = format_pos(posstr)
    cds_seq, tmp_pos_list, flag_gene_type, len_trna_type = merge_sequence(
        pos_list, seq)
    print('\n'+cds_seq)
    if maxnumber > 0:
        print('current pos:{}'.format(tmp_pos_list))
        print('current pos:{}'.format(posstr))
    storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq)

    if (trans_flag and (flag_gene_type != 'trna')) or args.force_flag:  # 翻译
        tmp_flag, inter_number, acid = trans2acid(cds_seq)
        if maxnumber > 0:
            if tmp_flag != 1:
                print('\n[START CONDON]The correct start codon was found after {} searches / Total times: {}'.format(
                    len(loop_count_flag1), len(loop_count_flag)))
                if tmp_flag != 3:
                    print('\n[STOP CONDON]The correct stop codon was found after {} searches / Total times: {}'.format(
                        len(loop_count_flag3), len(loop_count_flag)))
        current_abs_path = os.getcwd()
        # #########################################################################################################################################
        # 第一层if else
        if tmp_flag == 0:  # 起始ok
            if len(posstr.split(';')) != len(tmp_pos_list):  # ??????
                if len(tmp_pos_list) == 2:
                    new_pos_str = tmp_pos_list[0].replace(
                        ':1', ':+')+';'+tmp_pos_list[1].replace(':1', ':+')  # 考虑跨首尾
                else:
                    new_pos_str = tmp_pos_list
            else:
                new_pos_str = posstr
            print('Correct Position: [{}]'.format(new_pos_str))
            if pro_file_name != 'NULL':
                with open(os.path.join(current_abs_path, pro_file_name+'.acid'), 'w') as f_handle:
                    f_handle.write(str(acid)+'\n')
        # #####################################################################################################################################
        # 第一层if else
        #
        #
        elif tmp_flag == 1:  # 起始错,这个优先要满足的条件就不对
            # posstr 能传到这里  形如 1-7:+;14020-14078:+
            # tmp_pos_list 也能传到这里  形如['1-7:1', '14020-14078:1']
            loop_count += 1
            loop_count_flag.append(000)
            loop_count_flag1.append(1)  # 每有一次查找,列表元素个数就+1
            if maxnumber > 0:
                print(
                    '\n[START CONDON]Start search......Times:{} / Total times:{}'.format(loop_count, len(loop_count_flag)))
            cds_seq = cds_seq[3:]  # 已经判断起始错误了,因此直接把序列剪掉前面3个碱基

            # ##############################################################
            # 定义为第二层if else
            if cds_seq[0:3] not in start_codon_list and maxnumber != 0:  # 20220805  如果为假查找，就不进行下一步了
                start_flag = False
                if maxnumber > 0:
                    print('old pos:{}'.format(tmp_pos_list))
                    print('old pos:{}'.format(posstr))
                # 20220808 以下自动返回位置，也就是开头往后挪6bp
                if posstr.split(':')[-1] == '+':
                    new_pos_str = posstr.replace(posstr.split(
                        '-')[0], str(int(posstr.split('-')[0])+6))
                elif posstr.split(':')[-1] == '-':
                    new_pos_str = posstr.replace(re.findall(
                        r'\d+', posstr)[-1], str(int(re.findall(
                            r'\d+', posstr)[-1])-6))

                if loop_count <= maxnumber:
                    loop_look(infasta, new_pos_str, trans_flag,
                              loop_count, maxnumber, nuc_file_name, pro_file_name)
            elif cds_seq[0:3] in start_codon_list:
                start_flag = True
                print('\n'+cds_seq)
                # 20220808 以下自动返回正确位置，也就是开头往后挪3bp
                if posstr.split(':')[-1] == '+':
                    new_pos_str = posstr.replace(posstr.split(
                        '-')[0], str(int(posstr.split('-')[0])+3))
                elif posstr.split(':')[-1] == '-':
                    new_pos_str = posstr.replace(re.findall(
                        r'\d+', posstr)[-1], str(int(re.findall(
                            r'\d+', posstr)[-1])-3))
                tmp_flag, inter_number, acid = trans2acid(cds_seq)

                if maxnumber > 0:
                    if tmp_flag != 1:
                        print('The correct starting codon was found after {} searches / Total times: {}'.format(
                            len(loop_count_flag1), len(loop_count_flag)))
                if tmp_flag == 0:  # 20221020 其他地方还会出错,所以要再次检查
                    print('Correct Position: [{}]'.format(new_pos_str))
                    if pro_file_name != 'NULL':
                        with open(os.path.join(current_abs_path, pro_file_name+'.acid'), 'w') as f_handle:
                            f_handle.write(str(acid)+'\n')
                else:  # 还没完全找对
                    if loop_count <= maxnumber:
                        loop_look(infasta, new_pos_str, trans_flag, loop_count,
                                  maxnumber, nuc_file_name, pro_file_name)
        # ################################################################################################################################
        # 第一层if else
            """考虑细分情况 20221020考虑终止子错误的查找"""
        elif tmp_flag == 3:  # 3代表未终止 序列长度为3种情况,3n 3n+1 3n+2
            # posstr 能传到这里  形如 1-7:+;14020-14078:+
            # tmp_pos_list 也能传到这里  形如['1-7:1', '14020-14078:1']
            loop_count_flag.append(000)
            loop_count_flag3.append(3)
            loop_count = len(loop_count_flag)
            if maxnumber > 0:
                print(
                    '\n[STOP CONDON]Start search......Times:{} / Total times:{}'.format(len(loop_count_flag3), len(loop_count_flag)))
                print('old pos:{}'.format(tmp_pos_list))
                print('old pos:{}'.format(posstr))
            # if (abs(int(re.findall(r'\d+', posstr)[0])-int(re.findall(r'\d+', posstr)[1]))+1) % 3 == 0:
            # print(0)
            if len(cds_seq) % 3 == 0:
                if maxnumber > 0:
                    print(
                        'old len(sequence) is a multiple of three! {}=3n'.format(len(cds_seq)))
                if posstr.split(':')[-1] == '+':
                    new_pos_str = posstr.replace(posstr.split(
                        '-')[-1].split(':')[0], str(int(posstr.split('-')[-1].split(':')[0])+3))
                elif posstr.split(':')[-1] == '-':
                    new_pos_str = posstr.replace(re.findall(
                        r'\d+', posstr)[0], str(int(re.findall(
                            r'\d+', posstr)[0])-3))
            else:
                if len(cds_seq) % 3 == 1:
                    if maxnumber > 0:
                        print(
                            'old len(sequence) not a multiple of three! {}=3n+1'.format(len(cds_seq)))
                elif len(cds_seq) % 3 == 2:
                    if maxnumber > 0:
                        print(
                            'old len(sequence) not a multiple of three! {}=3n+2'.format(len(cds_seq)))
                if posstr.split(':')[-1] == '+':
                    new_pos_str = posstr.replace(posstr.split(
                        '-')[-1].split(':')[0], str(int(posstr.split('-')[-1].split(':')[0])+1))
                elif posstr.split(':')[-1] == '-':
                    new_pos_str = posstr.replace(re.findall(
                        r'\d+', posstr)[0], str(int(re.findall(
                            r'\d+', posstr)[0])-1))

            if loop_count <= maxnumber:
                loop_look(infasta, new_pos_str, trans_flag,
                          loop_count, maxnumber, nuc_file_name, pro_file_name)

        else:
            maxnumber == 0  # 赋值为0  相当于一个假查找
            # 一般是在外面赋值,这里因为修改tmp_falg==1的查找,不能在参数设置时默认为0,否则会影响tmp_falg==1的情况
            loop_count += 1
            print('第{}次查找中'.format(loop_count))
            new_pos_str = '124353-124892:-;126001-126552:-'

            if loop_count <= maxnumber:
                loop_look(infasta, new_pos_str, trans_flag,
                          loop_count, maxnumber, nuc_file_name, pro_file_name)
            else:
                print('{}次查找未有结果,取消第{}次查找'.format(loop_count-1, loop_count))
    #print('Already Run {}s'.format(time.time()-begin_time))
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
    loop_count_flag1 = []
    loop_count_flag3 = []
    loop_count_flag = []
    tmp_pos_list, inter_number = loop_look(
        args.infasta, args.posstr, args.trans_flag, loop_count, args.maxnumber, args.nuc_file_name, args.pro_file_name)
    if type(inter_number) == type(1):
        get_current_first_end_pos(tmp_pos_list, inter_number)
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
