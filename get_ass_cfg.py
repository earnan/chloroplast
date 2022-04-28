#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   get_ass_cfg.py
#         Author:   yujie
#    Description:   get_ass_cfg.py
#        Version:   1.0
#           Time:   2022/04/28 16:20:22
#  Last Modified:   2022/04/28 16:20:22
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import os
import argparse
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   get ass.cfg')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[input-dir]', type=str, help='例如/share/nas1/seqloader/yelvti/GP-20*_yelvti/2.cleandata', required=False)  # 默认clean_data
optional.add_argument('-o', '--output',
                      metavar='[output-file]', type=str, help='默认输出为当前目录ass.cfg文件', default='ass.cfg', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()


def find(path):
    dir_list = os.listdir(path)
    print('手动输入路径包含{0}'.format(str(dir_list)))
    for i in dir_list:

        out_file.write('{0}{1} :{2}'.format(
            '\t', i.replace('-', '_'), '\n'))  # 输出文件的第二行格式,物种名
        abs_path = os.path.join(path, i)  # abspath
        print('物种目录为{0}'.format(abs_path))

        file_list = os.listdir(abs_path)  # 列出子目录下的文件

        for j in range(len(file_list)):
            if file_list[j].endswith('.gz'):
                file_name = os.path.join(abs_path, file_list[j])
                print('即将写入{0}'.format(file_name))
                out_file.write(
                    '        read{0} : {1}{2}'.format(j+1, file_name, '\n'))  # 输出文件的第三四行格式


print('手动输入的路径是{0}'.format(args.input))
file_list = []
dir_list = []
out_file = open(args.output, 'a+')
out_file.write('Samples:' + '\n')  # 输出文件的第一行格式
find(args.input)

out_file.close()
print('done')
