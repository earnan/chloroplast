#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   ogdraw.py
#         Author:   yujie
#    Description:   ogdraw.py
#        Version:   1.0
#           Time:   2022/05/05 16:07:53
#  Last Modified:   2022/05/05 16:07:53
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
import re
import time
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   ogdraw.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='gbk文件', type=str,  required=False)
optional.add_argument(
    '-ir', '--irregion', metavar='[irregion]', help="形如'80005-101556,114348-135899'", type=str,  required=False)
optional.add_argument(
    '-p', '--prefix', metavar='[prefix]', help='prefix', type=str,  required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[outdir]', help='outdir', type=str, required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

drawgenemap = "/share/nas6/pub/pipline/genome-assembly-seq/chloroplast-genome-seq/v1.2/annotation/etc/plot_ogdraw_cpgenome_genbank/etc/GeneMap-1.1.1/bin/drawgenemap"
plastid_conf = "/share/nas6/pub/pipline/genome-assembly-seq/chloroplast-genome-seq/v1.2/annotation/etc/plot_ogdraw_cpgenome_genbank/conf/OGDraw_plastid_set.xml"
# .format(drawgenemap,args.infile,args.outdir,args.prefix,plastid_conf)
# .format(drawgenemap,args.infile,args.outdir,args.prefix,plastid_conf,args.irregion)
if args.irregion:
    cmd1 = " perl {0} --infile {1} --format png --density 600 --outfile {2}/{3}.dpi600.circular.png --useconfig {4} --ir_region {5}".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf, args.irregion)
    cmd2 = " perl {0} --infile {1} --format tif --density 600 --outfile {2}/{3}.dpi600.circular.tif --useconfig {4} --ir_region {5}".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf, args.irregion)
    cmd3 = " perl {0} --infile {1} --format png --density 300 --outfile {2}/{3}.dpi300.circular.png --useconfig {4} --ir_region {5}".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf, args.irregion)
    cmd4 = " perl {0} --infile {1} --format tif --density 300 --outfile {2}/{3}.dpi300.circular.tif --useconfig {4} --ir_region {5}".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf, args.irregion)
    cmd9 = " perl {0} --infile {1} --format ps --outfile {2}/{3}.circular --useconfig {4} --ir_region {5}".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf, args.irregion)
else:
    cmd1 = " perl {0} --infile {1} --format png --density 600 --outfile {2}/{3}.dpi600.circular.png --useconfig {4} ".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
    cmd2 = " perl {0} --infile {1} --format tif --density 600 --outfile {2}/{3}.dpi600.circular.tif --useconfig {4} ".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
    cmd3 = " perl {0} --infile {1} --format png --density 300 --outfile {2}/{3}.dpi300.circular.png --useconfig {4} ".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
    cmd4 = " perl {0} --infile {1} --format tif --density 300 --outfile {2}/{3}.dpi300.circular.tif --useconfig {4} ".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
    cmd9 = " perl {0} --infile {1} --format ps --outfile {2}/{3}.circular --useconfig {4}".format(
        drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)

cmd5 = " perl {0} --infile {1} --format png --density 600 --outfile {2}/{3}.dpi600.linear --useconfig {4} --force_linear ".format(
    drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
cmd6 = " perl {0} --infile {1} --format tif --density 600 --outfile {2}/{3}.dpi600.linear --useconfig {4} --force_linear ".format(
    drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
cmd7 = " perl {0} --infile {1} --format png --density 300 --outfile {2}/{3}.dpi300.linear --useconfig {4} --force_linear ".format(
    drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
cmd8 = " perl {0} --infile {1} --format tif --density 300 --outfile {2}/{3}.dpi300.linear --useconfig {4} --force_linear ".format(
    drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)
cmd10 = " perl {0} --infile {1} --format ps --outfile {2}/{3}.linear --useconfig {4} --force_linear ".format(
    drawgenemap, args.infile, args.outdir, args.prefix, plastid_conf)

print('dpi600.circular.png')
os.system(cmd1)
print('dpi600.circular.tif')
os.system(cmd2)
print('dpi300.circular.png')
os.system(cmd3)
print('dpi300.circular.tif')
os.system(cmd4)
print('circular.ps')
os.system(cmd9)

print('dpi600.linear.png')
os.system(cmd5)
print('dpi600.linear.tif')
os.system(cmd6)
print('dpi300.linear.png')
os.system(cmd7)
print('dpi300.linear.tif')
os.system(cmd8)
print('linear.ps')
os.system(cmd10)
