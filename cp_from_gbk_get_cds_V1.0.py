#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_from_gbk_get_cds.py
# Original Author:
#    Description:   cp_from_gbk_get_cds.py
#        Version:   1.0
#           Time:   2022/03/09 15:21:51
#  Last Modified:   2022/03/09 15:21:51
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import argparse
from Bio import SeqIO
import os
import re

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   cp_from_gbk_get_cds.py')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[dir]', help='输入gbk所在目录', type=str, default="F:/Hibiscus_sabdariffa/111", required=False)
optional.add_argument('-o', '--output',
                      metavar='[dir]', help='输出的路径', type=str, default="F:/Hibiscus_sabdariffa/out", required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note + format_seq + "\n"


def get_cds(gbk_file, f_cds):
    seq_record = SeqIO.read(gbk_file, "genbank")
    # print(seq_record.description.split('mitochondrion')[0])#调试处
    complete_seq = str(seq_record.seq)
    complete_note = (">" + seq_record.id + '_' +
                     (seq_record.description.split('mitochondrion')
                      [0]).replace(' ', '_')).rstrip('_') + "\n"
    complete_fasta = format_fasta(complete_note, complete_seq, 70)  # 70换行本例不采用
    count = 0
    cds_fasta = ""
    gene_name_count_list = []
    for ele in seq_record.features:
        if ele.type == "CDS":
            count += 1
            # print(ele.qualifiers)
            cds_seq = ""
            tmp_list = []
            # print(ele.location.parts)
            for ele1 in ele.location.parts:
                tmp_list.append(re.findall(
                    r'\d+', str(ele1.start))[0])  # 取位置出来
                tmp_list.append(re.findall(r'\d+', str(ele1.end))[0])  # 取位置出来
                cds_seq += complete_seq[ele1.start:ele1.end]

            cds_note = ">" + seq_record.id + \
                " [" + str(int(tmp_list[0])+1)+".." + tmp_list[-1]+"]" + \
                " [gene=" + ele.qualifiers['gene'][0] + "]" + \
                "\n"  # '>'后的格式和已有脚本兼容
            gene_name_count_list.append(ele.qualifiers['gene'][0])
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
            # print(cds_note.strip())

            if (f_cds):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
    s = '文件{0}有{1}个CDS'.format(os.path.basename(gbk_file), count)
    print(s)
    print(gene_name_count_list)
    return cds_fasta, complete_fasta, count, os.path.basename(gbk_file), s, gene_name_count_list


if __name__ == '__main__':
    #all_gene_list = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2','COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    # 文件输出路径
    out_cds_file_path = os.path.join(args.output, 'cds.fasta')
    out_complete_file = os.path.join(args.output, 'complete.fasta')
    out_log_file = os.path.join(args.output, 'log')
    # genbank 文件路径
    genbank_dir_path = args.input
    out_cds_file_path_obj = open(out_cds_file_path, "w")
    out_complete_file_obj = open(out_complete_file, "w")
    out_log_file_obj = open(out_log_file, 'w')
    count_dict = {}
    missing_gene_dict = {}
    out_log_file_obj.write(
        'gene{0}ATP6{0}ATP8{0}CYTB{0}COX1{0}COX2{0}COX3{0}ND1{0}ND2{0}ND3{0}ND4{0}ND4L{0}ND5{0}ND6\n\n'.format('\t'))

    path_list = os.listdir(genbank_dir_path)
    path_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字
    for file in path_list:
        # cds_fasta, complete_fasta = get_cds(genbank_dir_path + os.sep + file, False)#另一种写法
        (cds_fasta, complete_fasta, count, file_name, s, gene_name_count_list) = get_cds(
            os.path.join(genbank_dir_path, file), False)
        count_dict[file_name] = count

        out_cds_file_path_obj.write(cds_fasta)
        out_complete_file_obj.write(complete_fasta)
        # out_log_file_obj.write(s+'\n')
        out_log_file_obj.write('>'+file.rstrip('.gbk')+'\t')
        missing_gene_list = []
        """
        for ele in all_gene_list:
            if ele in gene_name_count_list:
                out_log_file_obj.write(ele+'\t')
            else:
                out_log_file_obj.write('NULL'+'\t')
                missing_gene_list.append(ele)
                missing_gene_dict['>'+file.rstrip('.gbk')] = missing_gene_list
        """
        out_log_file_obj.write('\n')
        [out_log_file_obj.write(tmp+'\t') for tmp in missing_gene_list]
        out_log_file_obj.write('\n')

    # print(count_dict)
    # out_log_file_obj.write(str(count_dict))
    print(missing_gene_dict)
    out_log_file_obj.write(str(missing_gene_dict))
    out_cds_file_path_obj.close()
    out_complete_file_obj.close()
    out_log_file_obj.close()
    print('Done')
