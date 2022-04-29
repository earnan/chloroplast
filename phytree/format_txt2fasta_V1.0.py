#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   format_txt2fasta_V1.0.py
#         Author:   yujie
#    Description:   format_txt2fasta_V1.0.py
#        Version:   1.0
#           Time:   2022/04/29 14:31:07
#  Last Modified:   2022/04/29 14:31:07
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import argparse
import os
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   格式化为fa标准文件,追加修改内容(替换)')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[file/dir]', help='目录的话末尾不要加\\,linux不加/', type=str, required=False)
optional.add_argument(
    '-c', '--check', metavar='[适用于linux]', type=bool, help="linux下使用'-c 1'", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def formatting1(file, abs_path):  # 绝对路径
    with open(file, 'r+') as f:
        content = f.read()
        f.seek(0, 0)  # 开头偏移0,即开头
        if args.check:  # linux运行
            seq_name = str(file).replace(abs_path+'/', '').replace('.seq', '')
            print('linux')
        else:  # 否则说明在windows下运行
            seq_name = str(file).replace(abs_path+'\\', '').replace('.seq', '')
            print('windows')
        if seq_name.find("("):
            seq_name = seq_name.split(')')[0].split('(')[1]  # PA3-10
            seq_name = seq_name.split('-')[1]  # 10
        f.write('>'+seq_name+'\n'+content.strip()+'\n')
        print(seq_name)
    return 0


def get_judge_from_str(abs_path):  # 判断是目录还是单文件,对每个文件进行同样操作
    item_list = []

    if os.path.isfile(abs_path):
        file = abs_path
        fa = open(file, 'r')
        line = fa.readlines()[0]
        fa.close()
        if line.startswith('>'):
            print(file+' ok')
        else:
            formatting1(file, abs_path)

    elif os.path.isdir(abs_path):
        print('目录')
        for item in os.listdir(abs_path):
            item = os.path.join(abs_path, item)  # 确保是绝对路径
            item_list.append(item)
            get_judge_from_str(item)

    return len(item_list)  # 文件个数


if __name__ == '__main__':  # 该块程序用于安全测试,当整个文件作为模块导入使用时,这块不会运行
    path = args.input
    abs_path = os.path.abspath(path)  # 确保是绝对路径
    number = get_judge_from_str(abs_path)
    print(number)
    print('Done')
