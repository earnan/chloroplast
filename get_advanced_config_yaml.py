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
    '-i1', '--indir1', metavar='[ass dir]', help='analysis/ or sample/', type=str, default='E:/', required=False)
optional.add_argument(
    '-i2', '--indir2', metavar='[ref dir]', help='ref_adv/', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[out dir]', help='default indir', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

# 传入参数
ass_dir_path = args.indir1  # 组装
ref_dir_path = args.indir2  # 参考
out_dir_path = args.outdir  # 输出
ref_gbk_path = os.path.join(ref_dir_path, 'gbk')
ref_genome_path = os.path.join(ref_dir_path, 'fasta')
# 初始化
ass_gbk_list = []
ass_genome_list = []
ref_gbk_list = []
ref_genome_list = []
# 列表赋值
ref_gbk_list = os.listdir(ref_gbk_path)
ref_gbk_list.sort()

ref_genome_list = os.listdir(ref_genome_path)
ref_genome_list.sort()
'''
ass_gbk_list = [1, 2]
ass_genome_list = [1, 2]
ref_gbk_list = [1, 2]
ref_genome_list = [1, 2]
'''

outfile_path = args.outdir+'test.yaml'
with open(outfile_path, 'a+') as outfile_handle:
    outfile_handle.write('Assembly:\n\tGbk:\n')
    for i in ass_gbk_list:
        outfile_handle.write('\t\t- {}\n'.format(i))

    outfile_handle.write('\tGenome:\n')
    for i in ass_genome_list:
        outfile_handle.write('\t\t- {}\n'.format(i))

    outfile_handle.write('Ref:\n\tGbk:\n')
    for i in ref_gbk_list:
        outfile_handle.write('\t\t- {}\n'.format(i))

    outfile_handle.write('\tGenome:\n')
    for i in ref_genome_list:
        outfile_handle.write('\t\t- {}\n'.format(i))

    outfile_handle.write('Library:\n\tIR  : {}\n'.format(i))
    outfile_handle.write('\tGene_anno:\n\t\t- {}\n'.format(i))
    outfile_handle.write('\n')
    outfile_handle.write(
        "#Fill manully. Root, if not defined ,fill 0. if defined,fill id. eq 'id1, id2'.")
    outfile_handle.write(
        '\tPhytree:\n\t\tFa  : {}\n\t\tRoot  : 0\n'.format('null'))

os.system("cat test.yaml")
print('Done')
