#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_extract2mafft.py
#         Author:   yujie
#    Description:   mt_extract2mafft.py
#        Version:   1.0
#           Time:   2022/03/28 17:00:58
#  Last Modified:   2022/03/28 17:00:58
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
# from icecream import ic
from Bio import SearchIO
import argparse
from Bio import SeqIO
import os
import re
import time

parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   mt_extract2mafft.py\n\
提取同名基因序列\n\
mafft比对\n\
V1.0串行运行')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[dir]', help='输入cds所在目录', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\out\\cds', required=False)
optional.add_argument('-o1', '--outdir1',
                      metavar='[dir]', help='序列提取后存放位置', type=str, default='F:\\ref_tre\\gene\\blast\\fasta', required=False)
optional.add_argument('-o2', '--outdir2',
                      metavar='[dir]', help='比对好的序列', type=str, default='F:\\ref_tre\\gene\\mafft', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='是否运行mafft?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

#################################################################
# 格式化成2016-03-20 11:45:39形式
begin_time = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('Start Time : {}'.format(start_time))
#################################################################
"""模板函数"""


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        # format_seq += "\n"
    return note.strip() + "\n" + format_seq + "\n"


def read_fasta_to_dic1(infasta):  # 最简单,针对普通fasta文件 >物种名
    with open(infasta, 'r') as f:
        seq_id = ''
        id_index = []
        dict_seq = {}
        for line in f:
            # 如果是">"开头的，就创建一个key键
            if line.startswith('>'):
                seq_id = line.strip('\n')  # ID为键
                id_index.append(line.replace(
                    "\n", "").replace(">", ""))  # 顺便创建索引的列表
                dict_seq[seq_id] = ''  # 有key无value
            # 如果不是">"开头的，在这个key键下添加value值
            else:
                dict_seq[seq_id] += line.strip('\n')
    print(len(dict_seq))  # 包含说明行序列行的字典结构
    return dict_seq


# 格式化基因名字,可重复使用,首次出现于mt_from_gbk_get_cds.py
def gene_name_standardization(gene_name):
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 13:
            if all_gene_list_lower[i] == gene_name:
                gene_name = all_gene_list_upper[i]
                break
            else:
                i += 1
        if i >= 13:
            print(gene_name)
            print('WARNING!Please check!')
    return gene_name


if __name__ == '__main__':
    if args.flag1:
        all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                               'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
        all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                               'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
        dict_gene_id_seq = {}
        for i in all_gene_list_upper:
            dict_gene_id_seq[i] = ''

        file_list = os.listdir(args.input)
        file_list.sort()
        for file in file_list:
            species_id = file.replace('_cds.fasta', '')
            # ic(species_id)
            infasta = os.path.join(args.input, file)
            with open(infasta, 'r') as fi:
                dict_seq = read_fasta_to_dic1(infasta)
                for j in dict_seq.keys():
                    gene = j.split()[-1].split('=')[-1].rstrip(']')
                    gene = gene_name_standardization(gene)
                    dict_gene_id_seq[gene] += format_fasta(j, dict_seq[j], 70)
                    # ic(j)
                    """#利用python本身的index函数,不通用,舍去
                    if gene in all_gene_list_upper:
                        # ic()
                        n = all_gene_list_upper.index(gene)
                        ic(n, gene)
                        dict_gene_id_seq[gene].append([j, dict_seq[j]])
                    elif gene in all_gene_list_lower:
                        n = all_gene_list_lower.index(gene)
                        gene = all_gene_list_upper[n]
                        ic(n, gene)
                        dict_gene_id_seq[gene].append([j, dict_seq[j]])
                    else:
                        print('warning {}'.format(gene))
                        list_warning.append(gene)
                    """
        all_gene_list_upper2 = ['COX1', 'ND1', 'ND2', 'ND4L', 'COX2',
                                'CYTB', 'ATP8', 'ND4', 'ATP6', 'ND3', 'ND5', 'ND6', 'COX3']
        all_gene_list_lower2 = ['cox1', 'nad1', 'nad2', 'nad4L', 'cox2',
                                'cob', 'atp8', 'nad4', 'atp6', 'nad3', 'nad5', 'nad6', 'cox3']
        n = 0
        for i in all_gene_list_upper2:
            n += 1
            filename = 'gene{0}.{1}.fasta'.format(n, i)
            with open(os.path.join(args.outdir1,  filename), 'wb') as f:
                f.write(dict_gene_id_seq[i].encode())

    if args.flag2:
        #########################################
        file_list1 = os.listdir(args.outdir1)
        file_list1.sort()
        for file1 in file_list1:
            infasta1 = os.path.join(args.outdir1, file1)
            cmd1 = "mafft --auto {0} > {1}/{2}.aln".format(
                infasta1, args.outdir2, file1)
            # print(cmd1)
            os.system(cmd1)
    #########################################
        file_list2 = os.listdir(args.outdir2)
        file_list2.sort()
        for file2 in file_list2:
            inaln2 = os.path.join(args.outdir2, file2)
            cmd2 = "perl /share/nas6/xul/program/mt2/phytree/gene_tree/src/fasta2line.pl -i {0} -o {1}/{2}".format(
                inaln2, args.outdir2, file2.replace('.aln', ''))
            # print(cmd2)
            os.system(cmd2)

        cmd3 = "rm {0}/*.aln".format(args.outdir2)
        # print(cmd3)
        os.system(cmd3)


###############################################################
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
print('Done')
###############################################################
