#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   get_advanced_config_yaml.py
#         Author:   yujie
#    Description:   get_advanced_config_yaml.py
#        Version:   1.0
#           Time:   2022/08/19 09:55:22
#  Last Modified:   2022/08/19 09:55:22
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
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   get_advanced_config_yaml.py\n\
#         Author:   yujie\n\
#    Description:   get_advanced_config_yaml.py\n\
#        Version:   1.0\n\
#           Time:   2022/08/19 09:55:36\n\
#  Last Modified:   2022/08/19 09:55:36\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
#\n\
##########################################################\n\
\n\
\npython3   get_advanced_config_yaml.py\n\
功能：\n\
1.常规使用\n\
1.1 -i [ ] -o [ ] \n\
2.其他使用\n\
2.1 -i [ ] -o [ ] \n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\chloroplast\get_advanced_config_yaml.py\n\
Path: /share/nas1/yuj/script/chloroplast/get_advanced_config_yaml.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--instr', metavar='[instr]', help='instr', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()
