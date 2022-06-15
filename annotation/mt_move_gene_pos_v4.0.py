#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_move_gene_pos_v3.0.py
#         Author:   yujie
#    Description:   mt_move_gene_pos_v3.0.py
#        Version:   4.0
#           Time:   2022/06/07 09:38:21
#  Last Modified:   2022/06/07 18:38:21
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
\npython3   线粒体平移基因修改位置V4.0\n\
不分开操作,则n1=0,n2有效\n\
即输入-n2 x 即可\n\
E:\OneDrive\jshy信息部\Script\mitochondrion\annotation\mt_move_gene_pos_v4.0.py\n\
V4.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[file]', type=str, help='要修改的文件,需输入', required=False)
optional.add_argument(
    '-fa', '--infasta', metavar='[file]', type=str, help='要修改的fsa,使用时需要输入', required=False)
optional.add_argument(
    '-o', '--output', metavar='[file]', type=str, help='输出文件,需输入', required=False)
optional.add_argument(
    '-ln', '--line_number', metavar='[int]', type=int, help='从第几行开始分开操作,默认不分开', default=1,  required=False)
optional.add_argument(
    '-n1', '--number1', metavar='[int]', type=int, help='平移距离1,可能不用输入',  required=False)
optional.add_argument(
    '-n2', '--number2', metavar='[int]', type=int, help='平移距离2,一定要输入',  required=False)
optional.add_argument(
    '-h1', '--info', metavar='[完整的帮助信息]', type=bool, help="使用时'-h1 1'即可", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()
if args.info:
    print(os.path.abspath(__file__))
    print("\n以下是帮助信息: ")
    print("     适用于内含子外显子")
    print("     #20220210第三版考虑了分段操作,仍需进一步排序\n      分段操作即修改环状的起点\n      此外第三版都改为了子函数")
    print("     增加了排序函数")
    print("     #20220607考虑跨首尾基因排序到开头,其自己的计数")
    print("\n")


# 替换多个指定位置上的字符(原字符串,位置列表,新的字符列表),使用时位置列表按自然语言习惯

def multi_replace(origin_string, pos_list, new_str_list):  # 把 1-10:- 变为 1-10:minus
    origin_string_list = []
    for s in origin_string:
        origin_string_list.append(s)
    new_string_list = origin_string_list
    # enumerate(sequence, [start=0]) start定义角标从哪开始,默认为0,同索引计数一致
    for i, pos_info in enumerate(pos_list):
        new_string_list[pos_info-1] = new_str_list[i]
    new_string = ''.join(new_string_list)
    return new_string


def edit_pos(pos_info, n):  # 1-10:+ 修改为 101-110:+
    if pos_info.endswith('-'):  # 仅针对负链方向,会被修改为minus
        pos_info = multi_replace(pos_info, [0], ['minus'])  # 索引-1对应自然的0
    pos_content = pos_info.split('-')
    new_start = int(pos_content[0])+n  # 例如1
    end_pos_all = pos_content[1]  # 例如10:+
    end_pos_content = end_pos_all.split(':')
    new_end = int(end_pos_content[0])+n  # 例如10
    new_pos_info = (str(new_start)+'-'+str(new_end)+':' +
                    end_pos_content[1]).replace('minus', '-')
    return new_pos_info  # 101-110:+


def get_new_line(line, n, max_len):
    # 20220607
    # max_len =
    line_content = line.split()
    pos_info = line_content[1]  # 例如67574-67687:-;135718-135975:+
    s_content = pos_info.split(';')  # 分号隔开
    if len(s_content) == 1:
        new_pos_info = edit_pos(pos_info, n)
        if int(new_pos_info.split('-')[0]) > max_len:  # 20220610考虑没有内含子的跨首尾的基因
            start = '1'
            end = new_pos_info.split('-')[1]
            if new_pos_info.endswith('-'):  # split会把末尾的 负链也当成分隔符
                new_pos_info = start+'-'+end+'-'
            else:
                new_pos_info = start+'-'+end
    elif len(s_content) == 2:
        new_pos_info1 = edit_pos(s_content[0], n)
        new_pos_info2 = edit_pos(s_content[1], n)
        if int(new_pos_info1.split('-')[0]) > max_len:  # 20220607考虑跨首尾的基因
            if new_pos_info2.endswith('-'):
                new_pos_info = '1'+'-'+new_pos_info2.split('-')[1]+'-'
            else:
                new_pos_info = '1'+'-'+new_pos_info2.split('-')[1]
        else:
            new_pos_info = new_pos_info1+';'+new_pos_info2
    elif len(s_content) == 3:
        new_pos_info1 = edit_pos(s_content[0], n)
        new_pos_info2 = edit_pos(s_content[1], n)
        new_pos_info3 = edit_pos(s_content[2], n)
        new_pos_info = new_pos_info1+';'+new_pos_info2+';'+new_pos_info3
    new_line = line.replace(pos_info, new_pos_info)
    print(new_line.strip('\n'))
    return new_line


# #############################################################################################挪动碱基
if args.line_number == 1 and args.number2 > 0 and args.infasta:  # 仅考虑把末尾n2 bp碱基挪到开头
    abs_path = os.path.abspath(args.infasta)
    indir_path = os.path.dirname(abs_path)
    file_prefix = os.path.basename(args.infasta).split('.')[0]
    with open(args.infasta, 'r') as fi_handle:
        tmp_dict = {}
        seq_id = fi_handle.readline()
        seq = fi_handle.readline().strip()  # 注意末尾有换行
        max_len = len(seq)
        tmp_dict[seq_id] = seq
    """末尾挪到开头"""
    n2 = args.number2
    last = seq[-n2:]
    # print(last+seq.rstrip(last)) #可能会产生bug
    s = ''
    tmp_list = list(seq)
    for i in range(n2):
        tmp_list.pop()  # pop函数默认返回被删除的值  直接用就好
    for i in tmp_list:
        s += i
    new_seq = last+s
    # print(s)
    with open(os.path.join(indir_path, file_prefix+'.fsa2'), 'w') as fo_handle:
        fo_handle.write(seq_id)
        fo_handle.write(new_seq+'\n')

# #############################################################################处理注释信息

"""平移的主函数"""
fi = open(args.input, 'r')
abs_path = os.path.abspath(args.input)
indir_path = os.path.dirname(abs_path)
tmp_path = os.path.join(indir_path, 'tmp')
tmp = open(tmp_path, 'w')
ln = args.line_number
n1 = args.number1
n2 = args.number2
for i in range(0, ln-1):  # 前(ln-1)行平移n1距离,从第ln行开始平移n2距离
    line = fi.readline()
    new_line = get_new_line(line, n1, max_len)
    tmp.write(new_line)
for line in fi:
    if line.strip() != '':  # 20220609 考虑输入的注释信息 下面有几行空行
        new_line = get_new_line(line, n2, max_len)
        tmp.write(new_line)
print('\n')
fi.close()
tmp.close()


# 排序的子函数
def pos_sort(input_file, output_file):
    fi = open(input_file, 'r')  # 打开tmp
    fo = open(output_file, 'w')
    dic = {}
    for line in fi:
        id = int(line.split()[1].split('-')[0])  # 位置的第一个起点
        dic[id] = line

    """计数"""
    count, cds_n, trn_n, rrn_n, dloop_n, ol_n = 0, 0, 0, 0, 0, 0  # 20220607考虑跨首尾基因排序到开头,其自己的计数
    """判断"""
    for i in sorted(dic):  # 排序
        line = dic[i].strip('\n')  # 去除了两端换行
        if line.startswith('rRNA'):  # rrna
            rrn_n += 1
            n = 'rRNA'+str(rrn_n)
        elif line.startswith('tRNA'):  # trna
            trn_n += 1
            n = 'tRNA'+str(trn_n)
        elif line.startswith('D-loop'):  # dloop
            dloop_n += 1
            n = 'D-loop'
        elif line.startswith('OL'):  # ol区  复制起始区域
            ol_n += 1
            n = 'OL'+str(ol_n)
        elif line.startswith('CDS'):  # cds
            cds_n += 1
            n = 'CDS'+str(cds_n)

        old_str = line.split('\t')[0]
        new_str = n
        print(line.replace(old_str, new_str))
        fo.write(line.replace(old_str, new_str)+'\n')
    return 0


pos_sort(tmp_path, args.output)
