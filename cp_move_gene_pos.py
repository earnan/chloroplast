#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_move_gene_pos.py
#         Author:   yujie
#    Description:   cp_move_gene_pos.py
#        Version:   1.0
#           Time:   2022/02/17 10:30:23
#  Last Modified:   2022/02/17 10:30:23
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import argparse
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   叶绿体ssc方向及基因修改V3.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[file]', type=str, help='要修改的文件,默认', default='gene.annotation.info', required=False)
optional.add_argument(
    '-o', '--output', metavar='[file]', type=str, help='输出文件,默认', default='final_gene.annotation.info', required=False)
optional.add_argument(
    '-ln', '--line_number', metavar='[int]', type=int, help='从第几行开始分开操作,默认不分开', default=1,  required=False)
optional.add_argument(
    '-n1', '--number1', metavar='[int]', type=int, help='平移距离1',  required=False)
optional.add_argument(
    '-n2', '--number2', metavar='[int]', type=int, help='平移距离2',  required=False)
optional.add_argument(
    '-h1', '--info', metavar='[完整的帮助信息]', type=bool, help="使用时'-h1 1'即可", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()
if args.info:
    print("\n以下是帮助信息: ")
    print("     适用于内含子外显子")
    print("     #20220210第三版考虑了分段操作,仍需进一步排序\n      分段操作即修改环状的起点\n      此外第三版都改为了子函数")
    print("     增加了排序函数")
    print("\n")
