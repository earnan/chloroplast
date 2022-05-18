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
import time

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   cp_from_gbk_get_cds.py')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[dir]', help='输入gbk所在目录', type=str, default='E:\\Examples\\cp_from_gbk_get_cds\\gbk', required=False)
optional.add_argument('-o', '--output',
                      metavar='[dir]', help='输出的路径', type=str, default='E:\\Examples\\cp_from_gbk_get_cds\\out', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note.strip() + "\n" + format_seq + "\n"


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
    cds_seq = ""
    tmp_list = []  # 位置列表
    for ele1 in ele.location.parts:
        if ele1.strand == (-1):
            # print('minus')
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际起点,从end中取不用+1
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际终点,从start取+1
            cds_seq += ir(complete_seq[ele1.start:ele1.end])
        elif ele1.strand == (1):
            # print('plus')
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际起点,要+1
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际终点,不用+1
            # 切片没问题,索引从start到end-1,也就是对应start+1到end的序列
            cds_seq += complete_seq[ele1.start:ele1.end]
    # print(tmp_list)
    return tmp_list, cds_seq


def get_complete_note(seq_record):  # 获取整个完整基因组ID
    seq_id = ''
    # if seq_record.description.find('chloroplast'):#有bug,用str格式化后就没问题了
    if str(seq_record.description).find('chloroplast') or seq_record.description.split(',')[-2].split()[-1] == 'chloroplast' or seq_record.description.split(',')[-2].split()[-1] == 'plastid':
        seq_id = seq_record.description.split(
            'chloroplast')[0].replace(' ', '_').rstrip('_')
        name = seq_record.name
        if seq_id == name:
            seq_id = seq_id
        elif seq_id != name:
            seq_id = seq_id+'_'+name
        complete_note = ">" + seq_id + "\n"  # chloroplast--叶绿体
    elif seq_record.description.split(',')[-2].split()[-1] == 'mitochondrion':
        seq_id = seq_record.description.split(
            'mitochondrion')[0].replace(' ', '_').rstrip('_')
        name = seq_record.name
        if seq_id == name:
            seq_id = seq_id
        elif seq_id != name:
            seq_id = seq_id+'_'+name
        complete_note = ">" + seq_id + "\n"  # mitochondrion--线粒体
    else:
        print('WARNING! {}!'.format(
            seq_record.description.split(', ')[-2].split()[-1]))
        complete_note = ">" + (seq_record.description.split('chloroplast')
                               [0]).replace(' ', '_').rstrip('_') + "\n"
    return complete_note, seq_id


def get_cds_note(ele, complete_seq, seq_id):  # 获取cds的id及序列
    if len(ele.location.parts) == 3:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+';' + tmp_list[2]+".." + tmp_list[3]+';' + \
            tmp_list[4]+".." + tmp_list[5]+"]" + " [gene=" + \
            ele.qualifiers['gene'][0] + "]" + "\n"  # '>'后的格式和已有脚本兼容
    elif len(ele.location.parts) == 2:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+';' + tmp_list[2]+".." + \
            tmp_list[3]+"]" + " [gene=" + ele.qualifiers['gene'][0] + \
            "]" + "\n"               # '>'后的格式和已有脚本兼容
    elif len(ele.location.parts) == 1:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+"]" + \
            " [gene=" + ele.qualifiers['gene'][0] + "]" + \
            "\n"    # '>'后的格式和已有脚本兼容
    return cds_note, cds_seq


def get_cds(gbk_file, flag):  # 解析gbk文件获取cds
    """完整基因组"""
    seq_record = SeqIO.read(gbk_file, "genbank")
    complete_seq = str(seq_record.seq)
    complete_note, seq_id = get_complete_note(seq_record)
    complete_fasta = format_fasta(complete_note, complete_seq, 80)  # 80换行本例不采用
    """cds序列"""
    count = 0  # 对cds数量计数
    cds_fasta = ""
    list_gene_name = []  # 统计cds
    for ele in seq_record.features:
        if ele.type == "CDS":
            count += 1
            cds_note, cds_seq = get_cds_note(ele, complete_seq, seq_id)
            list_gene_name.append(ele.qualifiers['gene'][0])
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
            if (flag):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
    s = '{0}有{1}个CDS'.format(os.path.basename(gbk_file), count)
    if count == 0:
        s = '-----------------------Warning!!! {0}有{1}个CDS-----------------------\n\
            -----------------------There may be no comments!!!-----------------------'.format(
            os.path.basename(gbk_file), count)
    print(s)
    print(list_gene_name)
    return cds_fasta, complete_fasta, count, os.path.basename(gbk_file), list_gene_name, s


if __name__ == '__main__':
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    begin_time = time.time()
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('Start Time : {}'.format(start_time))
    #################################################################
    dict_file_cds_count = {}  # 每个文件中cds计数
    file_list = os.listdir(args.input)
    file_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字
    if os.path.exists(args.output) == False:
        os.mkdir(args.output)
    for file in file_list:
        (cds_fasta, complete_fasta, count, file_name,  list_gene_name, s) = get_cds(
            os.path.join(args.input, file), False)
        dict_file_cds_count[file_name] = count  # 每个文件中cds计数
        with open((args.output+os.sep+file_name.rstrip('.gbk')+'_complete.fasta'), 'w') as f_complete, \
                open((args.output+os.sep+file_name.rstrip('.gbk')+'_cds.fasta'), 'w') as f_cds:
            f_cds.write(cds_fasta)
            f_complete.write(complete_fasta)
    print('\n')
    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
