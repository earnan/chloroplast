#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_pi_plot.py
#         Author:   yujie
#    Description:   cp_pi_plot.py
#        Version:   1.0
#           Time:   2022/08/25 09:44:25
#  Last Modified:   2022/08/25 09:44:25
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import pretty_errors
import re
import sys
import time
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
# \n\
# \n\
#       Filename:   cp_pi_plot.py\n\
#         Author:   yujie\n\
#    Description:   cp_pi_plot.py\n\
#        Version:   1.0\n\
#           Time:   2022/08/25 09:47:27\n\
#  Last Modified:   2022/08/25 09:47:27\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
# \n\
# \n\
\n\
\npython3   cp_pi_plot.py\n\
功能：\n\
1.常规使用\n\
1.1 -i [ ] -o [ ] \n\
2.其他使用\n\
2.1 -i [ ] -o [ ] \n\
\n\
# \n\
Path: E:\OneDrive\jshy信息部\Script\chloroplast\advance\cp_pi_plot.py\n\
Path: /share/nas1/yuj/script/chloroplast/advance/cp_pi_plot.py\n\
Version: 1.0\n\
# \n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

infile = 'F:\QA4159\data.txt'
with open(infile, 'r') as infile_handle:
    window_list = []
    midpoint_list = []
    pi_list = []
    n = 0

    for line in infile_handle:
        line = line.strip()
        if not line.startswith('Win') and line != '':
            content = line.split()
            window_list.append(content[0])
            midpoint_list.append(content[1])
            pi_list.append(round(float(content[2]), 2))
            n += 1
        if n > 100:
            break

# ic(pi_list)
# 折线图
x = window_list
# x = range(len(pi_list))   # 横坐标

# x3 = list_all_keys
y = pi_list  # 纵坐标

plt.plot(x, y,  color='r', label="pi")  # s-:方形
plt.xticks(rotation=90)
# plt.plot(x1, k1, 's-', color='b', label="high coverage")  # s-:方形
# plt.plot(x2, k2, 'o-', color='g', label="low coverage")  # o-:圆形

plt.xlabel("Nucleotide position (bp)")  # 横坐标名字
plt.ylabel("pi")  # 纵坐标名字
plt.title("title")  # 标题
plt.legend(loc="best")  # 图例
# plt.savefig("F:\QA4159\saved.png")  # 保存图片名称
plt.show()
