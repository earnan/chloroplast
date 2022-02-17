#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename: 编程处理模板.py
#         Author: yuj@genepioneer.cn
#    Description: sample
#  Last Modified: 2021-xx-xx 16:29:29
#
# Copyright (C) 2021xxxx genepioneer Corporation
##########################################################
# id.list文件中登录号一列 \t 物种名占两列
# 导入包
import os
import argparse
# 用户交互信息
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3 修改nwk树文件的物种名称')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-1', '--idlist', metavar='[file]', type=str, help='输入id.list',  required=False)
optional.add_argument(
    '-2', '--treenwk', metavar='[file]', type=str, help='默认sample.genome.nwk',  default='sample.genome.nwk',   required=False)
optional.add_argument(
    '-o', '--output', metavar='[file]', type=str, help='默认new.nwk',  default='new.nwk', required=False)


optional.add_argument(
    '-c', '--check', metavar='[检查,见注释]', type=bool, help="直接替换登录号,使用时'-c 1'即可", default='', required=False)
optional.add_argument(
    '-a', '--check2', metavar='[NC后为空格的情况]', type=bool, help="适用于bayes,使用时'-a 1'即可", default='', required=False)


optional.add_argument(
    '-id1', '--idfile1', metavar='[初始ID]', type=str, help="输入原id",   required=False)
optional.add_argument(
    '-id2', '--idfile2', metavar='[修改后的id]', type=str, help="输入新id",   required=False)


optional.add_argument(
    '-i', '--info', metavar='[完整的帮助信息]', type=bool, help="使用时'-i 1'即可", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()
if args.info:
    print("\n以下是帮助信息: ")
    print("适用于最大似然法及贝叶斯法的结果修改")
    print("登录号(带版本)+物种  改为  物种+登录号(带版本)")
    print("登录号(不带版本)  改为  物种+登录号(带版本)")
    print("#20220107修改,需要加一个判断,保证登录号都带有版本信息'.1'")
    print("#20220124增加一个读取文件为名称列表的子函数,程序有修改\n新增参数id1/id2")
    print("         简写名字 替换为 属名+种名+品系 \n")


# 以下为核心程序

# 打开文件,初始化
if args.treenwk and args.output:  # 第一部分功能的初始化
    if args.idlist:
        id_list_file = open(args.idlist, 'r')
        id_list_lines = id_list_file.readlines()
        id_list = []
        id_new_list = []
        count = []
    tree_nwk_file = open(args.treenwk, 'r')
    tree_nwk_line = tree_nwk_file.read()  # 直接读成一个长字符串
    output_file = open(args.output, 'w')


# 定义子函数
def get_id_list(id_list_lines):  # 第一二部分的子函数,判断args.idlist文件版本号以及重排获得新名字
    for id in id_list_lines:
        content = id.split()  # 空格切开
        # 登录号,20220107修改,需要加一个判断,保证登录号都带有版本信息".1"
        if content[0].find('.') < 0:  # 如果不带版本信息
            print('请注意{} '.format(id.strip('\n')))
            accession_0 = content[0]+'.1'
            accession = '_'+content[0]+'.1'
            id = id.replace(content[0], accession_0)
        else:
            accession = '_'+content[0]
        id_new = ''
        for i in range(len(content)):
            if i > 0:
                id_new_tmp = id_new+content[i]+'_'
                id_new = id_new_tmp  # 仅物种名
        id_new = id_new+accession  # 物种名+登录号(带v)
        id_new_list.append(id_new)

        if args.check:  # 如果仅有登录号(不带v)
            accession = accession.lstrip('_').rstrip('.1')  # 这一步变为 NC_028725
            id_list.append(accession)

        elif args.check2:  # 上面accession为_NC_028725.1
            id = id.strip('\n').replace("\t", "_").replace(
                " ", "_")  # 之前id已经带有版本,这一步变为 NC_028725.1_Tyrophagus_longior

            # 构建出 NC 028725 1 Tyrophagus longior
            id = id.replace('.', ' ').replace('_', ' ')
            id_list.append(id)

        else:
            id = id.strip('\n').replace("\t", "_").replace(
                " ", "_")  # 之前id已经带有版本,这一步变为 NC_028725.1_Tyrophagus_longior
            id_list.append(id)

    print('\n', id_list, '\n', '\n', id_new_list, '\n')
    return id_list, id_new_list


def replace_with_str(tree_nwk_line, output_file, id_list, id_new_list):  # 3种情况下最后一步通用替换函数
    for i in range(len(id_new_list)):
        print('第{0}次替换'.format(i+1))
        tree_nwk_tmp = tree_nwk_line.replace(id_list[i], id_new_list[i])
        tree_nwk_line = tree_nwk_tmp
    output_file.write(tree_nwk_line)
    return 0


def get_id_list_from_file(file):  # 第三种功能,将文件读取为要当参数传递的列表
    id_list = []
    f = open(file, 'r')
    for line in f:
        # print(line.strip())
        id_list.append(line.strip())
    print('\n'+str(id_list))
    print(len(id_list))
    return id_list


# 第三种功能
if args.idfile1 and args.idfile2 and args.treenwk and args.output:  # 原名字,新名字,旧可编辑树文件,新的树文件
    id_list = get_id_list_from_file(args.idfile1)
    id_new_list = get_id_list_from_file(args.idfile2)
    replace_with_str(tree_nwk_line, output_file, id_list, id_new_list)


# 第一种功能
if args.idlist and args.treenwk and args.output:
    (id_list, id_new_list) = get_id_list(id_list_lines)
    replace_with_str(tree_nwk_line, output_file, id_list, id_new_list)

# 关闭文件
if args.treenwk and args.output:
    if args.idlist:
        id_list_file.close()
    tree_nwk_file.close()
    output_file.close()
print('Done')

# 测试用
# if args.idfile1:
# get_id_list_from_file(args.idfile1)
