#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename: format_txt2fasta_V1.0.py
#         Author: yuj@genepioneer.cn
#    Description: sample
#  Last Modified: 2021-xx-xx 16:29:29
#
# Copyright (C) 2021xxxx genepioneer Corporation
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
    '-c', '--check', metavar='[判断系统]', type=bool, help="linux使用时'-c 1'即可", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def formatting(file):
    with open(file, 'r+') as f:
        content = f.read()
        f.seek(0, 0)  # 开头偏移0,即开头
        if args.check:
            seq_name = str(file).replace(path+'/', '').replace('.seq', '')
            print('linux')
        else:
            seq_name = str(file).replace(path+'\\', '').replace('.seq', '')
            print('windows')
        if seq_name.find("("):
            seq_name = seq_name.split(')')[0].split('(')[1]  # PA3-10
            seq_name = seq_name.split('-')[1]  # 10
        f.write('>'+seq_name+'\n'+content.strip()+'\n')
        print(seq_name)
    return 0


def get_judge_from_str(path):
    item_list = []

    if os.path.isfile(path):
        fa = open(path, 'r')
        line = fa.readlines()[0]
        fa.close()
        if line.startswith('>'):
            print(path+' ok')
        else:
            formatting(path)

    elif os.path.isdir(path):
        print('目录')
        for item in os.listdir(path):
            item = os.path.join(path, item)
            item_list.append(item)
            get_judge_from_str(item)

    return len(item_list)


if __name__ == '__main__':  # 该块程序用于安全测试,当整个文件作为模块导入使用时,这块不会运行
    path = args.input
    path = os.path.abspath(path)
    i = get_judge_from_str(path)
    print(i)
    print('Done')
