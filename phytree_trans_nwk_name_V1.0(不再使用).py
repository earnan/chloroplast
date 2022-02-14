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
import os
import argparse
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3 修改nwk树文件的物种名称')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i1', '--idlist', metavar='[file]', type=str, help='默认id.list', default='id.list', required=False)
optional.add_argument(
    '-c', '--check', metavar='[检查]', type=bool, help="直接替换登录号,使用时'-c 1'即可", default='', required=False)
optional.add_argument('-i2', '--treenwk',
                      metavar='[file]', type=str, help='默认sample.genome.nwk',  default='sample.genome.nwk', required=False)
optional.add_argument(
    '-o', '--output', metavar='[file]', type=str, help='默认new.nwk', default='new.nwk', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()


id_list_file = open(args.idlist, 'r')
id_list_lines = id_list_file.readlines()
id_list = []
id_new_list = []
count = []

tree_nwk_file = open(args.treenwk, 'r')
tree_nwk_line = tree_nwk_file.read()  # 直接读成一个长字符串

output_file = open(args.output, 'w')

for id in id_list_lines:
    content = id.split()  # 空格切开
    accession = '_'+content[0]
    id_new = ''
    for i in range(len(content)):
        if i > 0:
            id_new_tmp = id_new+content[i]+'_'
            id_new = id_new_tmp  # 只有物种名
    id_new = id_new+accession  # 物种名 登录号
    id_new_list.append(id_new)
    if args.check:
        accession = content[0]
        id_list.append(accession)
    else:
        id = id.strip('\n').replace("\t", "_").replace(" ", "_")  # 原有id用_连接
        id_list.append(id)  # 原有id用_连接

    """
    content = id.split('_')
    if len(content) == 3:
        id_new = content[-2]+'_'+content[-1]+'_'+content[0]
    else:
        id_new = content[-2]+'_'+content[-1] + '_'+content[0]+'_'+content[1]
    """


for i in range(len(id_list)):
    print('第{0}次替换'.format(i+1))
    tree_nwk_tmp = tree_nwk_line.replace(id_list[i], id_new_list[i])
    tree_nwk_line = tree_nwk_tmp
# print(tree_nwk_line)
output_file.write(tree_nwk_line)

id_list_file.close()
tree_nwk_file.close()
output_file.close()
print('Done')
