#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   phytree_trans_nwk_name_V3.0.py
#         Author:   yujie
#    Description:   phytree_trans_nwk_name_V3.0.py
#        Version:   3.0
#           Time:   2022/05/18 15:33:22
#  Last Modified:   2022/08/30 18:00:07
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
#import pretty_errors
import re
import sys
import time
#import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   phytree_trans_nwk_name_V3.0.py\n\
#         Author:   yujie\n\
#    Description:   phytree_trans_nwk_name_V3.0.py\n\
#        Version:   3.0\n\
#           Time:   2022/05/18 15:33:22\n\
#  Last Modified:   2022/08/30 18:00:46\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
#\n\
##########################################################\n\
\n\
\npython3   phytree_trans_nwk_name_V3.0.py\n\
Function:\n\
1.方法1 -f1 -1  [id.list]   -2  [sample.genome.nwk] -3  [out nwk]\n\
2.方法2 -f2 -1  [id.list]   -2  [sample.genome.nwk] -3  [out nwk]   -c1/-c2\n\
    2.1原始树只有登录号(无版本) -c1\n\
    2.2NC后为空格的情况,适用于bayes -c2\n\
3.方法3 -f3 -id1    [初始ID文件]    -id2    [修改后的id文件]    -2  [sample.genome.nwk] -3  [out nwk]\n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\chloroplast\phytree\phytree_trans_nwk_name_V3.0.py\n\
Path: /share/nas1/yuj/script/chloroplast/phytree/phytree_trans_nwk_name_V3.0.py\n\
Version: 3.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')

optional.add_argument(
    '-f1', '--function1', action='store_true', required=False)
optional.add_argument(
    '-1', '--idlist', metavar='[id.list]', type=str,  required=False)  # id.list文件中登录号一列 \t 物种名占两列
optional.add_argument(
    '-2', '--treenwk', metavar='[sample.genome.nwk]', type=str,  default='sample.genome.nwk',   required=False)
optional.add_argument(
    '-3', '--outfile', metavar='[out nwk]', type=str,   default='phytree.nwk', required=False)

optional.add_argument(
    '-f2', '--function2', action='store_true', required=False)
optional.add_argument(
    '-c1', '--check1', help="原始树只有登录号(无版本),使用则-c1", action='store_true', required=False)
optional.add_argument(
    '-c2', '--check2', help="NC后为空格的情况,适用于bayes,使用时-c2", action='store_true', required=False)

optional.add_argument(
    '-f3', '--function3', action='store_true', required=False)
optional.add_argument(
    '-id1', '--idfile1', metavar='[初始ID文件]', type=str,  required=False)
optional.add_argument(
    '-id2', '--idfile2', metavar='[修改后的id文件]', type=str,   required=False)

optional.add_argument(
    '-i', '--info', help="完整的帮助信息,运行则-i", action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()


# 第一二部分的子函数,判断args.idlist文件版本号以及重排获得新名字
def get_id_list_from_id_list(id_list_lines):
    # id.list # MT157619.1	Camellia petelotii var. microcarpa
    id_list = []  # 原来树文件的id列表
    id_new_list = []  # 新的树文件的id列表
    """对每一行id进行前处理,带上版本号"""
    for id in id_list_lines:
        content = id.strip('\n').split('\t')  # 制表符切开
        """ #登录号,20220107修改,需要加一个判断,保证登录号都带有版本信息".1" """
        if content[0].find('.') < 0:  # 如果不带版本信息
            print('\n请注意{}没有版本信息! '.format(id.strip('\n')))
            accession_0 = content[0]+'.1'
            accession = '_'+content[0]+'.1'
            # MT157619.1	Camellia petelotii var. microcarpa
            id = id.replace(content[0], accession_0)
        else:
            accession = '_'+content[0]

        """构建新的id"""
        id_new = ''
        for i in range(len(content)):
            if i > 0:
                id_new_tmp = id_new+content[i].replace(" ", "_")+'_'
                id_new = id_new_tmp  # 仅物种名
        id_new = id_new+accession  # 物种名+登录号(带v)
        id_new_list.append(id_new)

        """获取或者构建原始树文件的id"""
        if args.check1:  # 如果原始树仅有登录号(还不带v)
            accession = accession.lstrip('_').rstrip(
                '.1')  # 这一步将idlist里的id变为原id形式 NC_028725
            id_list.append(accession)

        elif args.check2:  # 这一部分主要是针对贝叶斯原始树中"NC"开头,其后直接为空格的情况
            # 上面accession为_NC_028725.1,id为MT157619.1	Camellia petelotii var. microcarpa
            id = id.strip('\n').replace("\t", "_").replace(
                " ", "_")  # 之前id已经带有版本,这一步变为 NC_028725.1_Tyrophagus_longior
            # 构建出 NC 028725 1 Tyrophagus longior
            id = id.replace('.', ' ').replace('_', ' ')
            id_list.append(id)

        else:  # 正常情况,原始树文件没啥幺蛾子的情况
            id = id.strip('\n').replace("\t", "_").replace(
                " ", "_")  # 之前id已经带有版本,这一步变为 NC_028725.1_Tyrophagus_longior
            id_list.append(id.rstrip('_'))

    return id_list, id_new_list

# 第一二部分的子函数,判断args.idlist文件版本号以及重排获得新名字


def get_id_list_from_id_list2(id_list_lines):
    # id.list # MT157619.1	Camellia petelotii var. microcarpa

    id_dict = {}  # 原来
    id_new_dict = {}  # 新的
    """对每一行id进行前处理,带上版本号"""
    for id in id_list_lines:
        content = id.strip('\n').split('\t')  # 制表符切开
        tmp_key = content[0].rstrip('.1').split('_')[-1]
        id_dict[tmp_key] = ''  # 原来
        id_new_dict[tmp_key] = ''   # 新的
        """ #登录号,20220107修改,需要加一个判断,保证登录号都带有版本信息".1" """
        if content[0].find('.') < 0:  # 如果不带版本信息
            print('\n请注意{}没有版本信息! '.format(id.strip('\n')))
            accession_0 = content[0]+'.1'
            accession = '_'+content[0]+'.1'
            # MT157619.1	Camellia petelotii var. microcarpa
            id = id.replace(content[0], accession_0)
        else:
            accession = '_'+content[0]

        """构建新的id"""
        id_new = ''
        for i in range(len(content)):
            if i > 0:
                id_new_tmp = id_new+content[i].replace(" ", "_")+'_'
                id_new = id_new_tmp  # 仅物种名
        id_new = id_new+accession  # 物种名+登录号(带v)
        id_new_dict[tmp_key] = id_new

        """获取或者构建原始树文件的id"""
        if args.check1:  # 如果原始树仅有登录号(还不带v)
            accession = accession.lstrip('_').rstrip(
                '.1')  # 这一步将idlist里的id变为原id形式 NC_028725
            id_dict[tmp_key] = accession

        elif args.check2:  # 这一部分主要是针对贝叶斯原始树中"NC"开头,其后直接为空格的情况
            # 上面accession为_NC_028725.1,id为MT157619.1	Camellia petelotii var. microcarpa
            id = id.strip('\n').replace("\t", "_").replace(
                " ", "_")  # 之前id已经带有版本,这一步变为 NC_028725.1_Tyrophagus_longior
            # 构建出 NC 028725 1 Tyrophagus longior
            id = id.replace('.', ' ').replace('_', ' ')
            id_dict[tmp_key] = id

        else:  # 正常情况,原始树文件没啥幺蛾子的情况
            id = id.strip('\n').replace("\t", "_").replace(
                " ", "_")  # 之前id已经带有版本,这一步变为 NC_028725.1_Tyrophagus_longior
            id_dict[tmp_key] = id.rstrip('_')

    return id_dict, id_new_dict


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        int(s)
        return True
    except ValueError:
        pass

    return False


# 20220706 新增子函数
def get_ori_id_list_from_tree2(tree_nwk_line):
    ori_id_dict = {}
    partten = r'([a-zA-Z#0-9_.]+)'
    list_all = re.findall(partten, tree_nwk_line)

    for i in list_all:
        if not is_number(i):
            tmp_key = i.rstrip('.1').split('_')[-1]
            ori_id_dict[tmp_key] = i
    return ori_id_dict


def replace_with_str(tree_nwk_line, output_file, id_list, id_new_list):  # 3种情况下最后一步通用替换函数
    for i in range(len(id_new_list)):
        print('第{0}次替换'.format(i+1))
        tree_nwk_tmp = tree_nwk_line.replace(id_list[i], id_new_list[i])
        tree_nwk_line = tree_nwk_tmp
    output_file.write(tree_nwk_line)
    return 0


def replace_with_str2(tree_nwk_line, output_file, id_dict, id_new_dict):  # 3种情况下最后一步通用替换函数
    tmp_count = 0
    for i in id_new_dict.keys():
        tmp_count += 1
        if i in id_dict.keys():
            tree_nwk_tmp = tree_nwk_line.replace(id_dict[i], id_new_dict[i])
            tree_nwk_line = tree_nwk_tmp
    print('{0}次替换'.format(tmp_count))
    output_file.write(tree_nwk_line)
    return 0

# 功能3,将文件直接读取为要当参数传递的列表,等同于get_id_list_from_id_list(id_list_lines)的构建


def get_id_list_from_file(file):
    id_list = []
    with open(file, 'r') as f:
        for line in f:
            id_list.append(line.strip())
    print('\n'+str(id_list)+'\n')
    print(len(id_list))
    return id_list


if __name__ == '__main__':

    if args.info:
        print("\n以下是帮助信息: ")
        print("适用于最大似然法及贝叶斯法的结果修改")
        print("\t登录号(带版本)+物种 -> 物种+登录号(带版本)")
        print("\t登录号(不带版本) -> 物种+登录号(带版本)")
        print("#20220107修改,需要加一个判断,保证登录号都带有版本信息'.1'")
        print("#20220124增加一个读取文件为名称列表的子函数,程序有修改\n\t新增参数id1/id2\n\t简写名字 -> 属名+种名+品系")
        print("#20220518 V3.0 版本,大幅修改原代码逻辑")
        print("#20220830 V3.1 版本,更新交互信息")
        sys.exit(0)  # 正常退出
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    begin_time = time.time()
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('Start Time : {}'.format(start_time))
    #################################################################

    # 功能1
    if args.function1:
        with open(args.idlist, 'r') as f_id, open(args.treenwk, 'r') as f_tree1, open(args.outfile, 'w') as f_out:
            tree_nwk_line = f_tree1.read()  # 直接读成一个长字符串
            id_list_lines = f_id.readlines()  # id_list_lines类型为列表

            count = []  # 计数
            (id_dict, id_new_dict) = get_id_list_from_id_list2(id_list_lines)
            # 20220706 被替换的id列表是从id.list文件获取的,可能与树文件中要被替换的并不一致
            # 因此  需要一个从树文件读取列表的函数
            ori_id_dict = get_ori_id_list_from_tree2(tree_nwk_line)
            """
            if tree_nwk_line.startswith('#NEXUS'):
                print("\nBayes! Try use '_' to replace '.' ")
                for i in range(len(id_list)):
                    id_list[i] = id_list[i].replace('.', '_')
            """
            print('\n', ori_id_dict, '\n', '\n', id_new_dict, '\n')
            replace_with_str2(tree_nwk_line, f_out, ori_id_dict, id_new_dict)

    # 功能2
    if args.function2:
        with open(args.idlist, 'r') as f_id, open(args.treenwk, 'r') as f_tree1, open(args.outfile, 'w') as f_out:
            tree_nwk_line = f_tree1.read()  # 直接读成一个长字符串
            id_list_lines = f_id.readlines()  # id_list_lines类型为列表
            id_list = []  # 原来树文件的id列表
            id_new_list = []  # 新的树文件的id列表
            count = []  # 计数
            (id_list, id_new_list) = get_id_list_from_id_list(id_list_lines)
            print('\n', id_list, '\n', '\n', id_new_list, '\n')
            replace_with_str(tree_nwk_line, f_out, id_list, id_new_list)

    # 功能3
    if args.function3:  # 原名字,新名字,旧可编辑树文件,新的树文件
        with open(args.treenwk, 'r') as f_tree1, open(args.outfile, 'w') as f_out:
            tree_nwk_line = f_tree1.read()  # 直接读成一个长字符串
            id_list = get_id_list_from_file(args.idfile1)
            id_new_list = get_id_list_from_file(args.idfile2)
            print('\n', id_list, '\n', '\n', id_new_list, '\n')
            replace_with_str(tree_nwk_line, f_out, id_list, id_new_list)

    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
