#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   cp_from_gbk_get_cds.py
#         Author:   yujie
#    Description:   cp_from_gbk_get_cds.py
#        Version:   3.0
#           Time:   2022/03/09 15:21:51
#  Last Modified:   2022/12/17 11:20:02
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from icecream import ic
import argparse
# import linecache
import os
# import pretty_errors
import re
import sys
import time
# import copy
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
python3 cp_from_gbk_get_cds_V3.0.py -i [gbk dir] -o [out dir]\n\
Each species generates three files of cds,trna and complete sequence.')
optional = parser.add_argument_group('optional')
required = parser.add_argument_group('required')
required.add_argument(
    '-i', '--input',  help='gbk dir path', type=str)
required.add_argument(
    '-o', '--output',  help='output dir path', type=str)
optional.add_argument('-c', '--check', help='create genes via gap,not running by default,input "-c" when using',
                      action='store_true', required=False)
optional.add_argument('-info', '--info', help='show update log and exit',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help',
                      help='show this help message and exit')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t20221217 feat: ✨ 对gbk文件进行去重')
    print('\t20221219 🐞fix(get_gene_note): 修改没有gene标签的cds的默认ID')
    print('\t20221220 ✨feat(main): 写入去重后的登录号')
    sys.exit(0)

# ############################################################################################


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        # format_seq += "\n"
    return note.strip() + "\n" + format_seq + "\n"


def ir(s):  # 反向互补
    re = s[::-1]  # 字符串反向
    c = ""  # 定义字符串c接收互补序列
    for i in re:
        if i == 'A':
            c = c + 'T'
        elif i == 'G':
            c = c + 'C'
        elif i == 'T':
            c = c + 'A'
        elif i == 'C':
            c = c + 'G'
    return c


def merge_sequence(ele, complete_seq):  # 合并获取到的序列
    gene_seq = ""
    tmp_list = []  # 位置列表
    for ele1 in ele.location.parts:
        if ele1.strand == (-1):
            # print('minus')
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际起点,从end中取不用+1
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际终点,从start取+1
            gene_seq += ir(complete_seq[ele1.start:ele1.end])
        elif ele1.strand == (1):
            # print('plus')
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际起点,要+1
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际终点,不用+1
            # 切片没问题,索引从start到end-1,也就是对应start+1到end的序列
            gene_seq += complete_seq[ele1.start:ele1.end]
    return tmp_list, gene_seq


def get_complete_note(seq_record):  # 获取整个完整基因组ID
    try:
        seq_id = ''
        # 20220819 NC_044756.1.gbk voucher Liu HM/CP02 chloroplast  有特殊符号，需要处理
        seq_record.description = seq_record.description.replace('/', '_')
        # if seq_record.description.find('chloroplast'):#有bug,用str格式化后就没问题了
        # 20220627 if str(seq_record.description).find('chloroplast') -1也成立,判断时一定要以False True为准
        # or seq_record.description.split(',')[-2].split()[-1] == 'chloroplast' or seq_record.description.split(',')[-2].split()[-1] == 'plastid':
        if str(seq_record.description).find('chloroplast') > 0 \
            or str(seq_record.description).find('plastid') > 0 \
                or seq_record.description.split(',')[-2].split()[-1] == 'chloroplast' \
                or seq_record.description.split(',')[-2].split()[-1] == 'plastid':  # gbk 属于 叶绿体
            if str(seq_record.description).find('chloroplast') > 0:  # 20220819 根据不同关键词分割
                seq_id = seq_record.description.split(
                    'chloroplast')[0].replace(' ', '_').rstrip('_')
            elif str(seq_record.description).find('plastid') > 0:
                seq_id = seq_record.description.split(
                    'plastid')[0].replace(' ', '_').rstrip('_')
            name = seq_record.name
            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"  # chloroplast--叶绿体
        elif seq_record.description.split(',')[-2].split()[-1] == 'mitochondrion':
            seq_id = seq_record.description.split(
                'mitochondrion')[0].replace(' ', '_').rstrip('_')
            name = seq_record.name
            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"  # mitochondrion--线粒体
        else:
            print('Genome Type WARNING! {}!'.format(
                seq_record.description.split(', ')[-2].split()[-1]))
            complete_note = ">" + (seq_record.description.split('chloroplast')
                                   [0]).replace(' ', '_').rstrip('_') + "\n"
    except:  # 如果遇到任何出错
        print('try/except')
        complete_note = ''
        # gbk_type = input('genome type(1:chloroplast;2:mitochondrion): ')
        gbk_type = 1
        if gbk_type == 1:
            seq_id = seq_record.description.split(
                'chloroplast')[0].replace(' ', '_').rstrip('_')  # 物种或样品名

            if seq_id.startswith('UNVERIFIED:_'):  # 去掉 UNVERIFIED:_
                seq_id = seq_id.lstrip('UNVERIFIED:_')
            # 去掉Cerion_watlingense_voucher_USNM:1514170_MN904501 中 冒号后的内容
            if len(seq_id.split(':')) > 1:
                seq_id = seq_id.split(':')[0]

            name = seq_record.name  # 要么是登录号  要么是样本

            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"

    return complete_note, seq_id


def get_gene_note(file_no, file_name, ele, complete_seq, seq_id, tmp_gene_name):  # 获取gene的id及序列
    '''
    传入上一个基因信息tmp_gene_name,返回当前基因信息tmp_gene_name
    OrderedDict(
        [
        ('codon_start', ['1']),
        ('gene', ['rps19']),
        ('product', ['ribosomal protein S19']),
        ('protein_id', ['UKP82335.1']),
        ('transl_table', ['11']),
        ('translation', [
         'MTRSLKKNPFVANHLLRKINKLNTKAEKEIIITWSRASTIIPTMIGHTIAIHNGKEHLPIYITDRMVGHKLGEFSPTLNFRGHAKNDNRSRR'])
        ]
        )
    '''
    if 'gene' not in ele.qualifiers.keys():
        # 返回上一个基因,好从其他参考找这个没名字的
        try:
            tmp_gene_name = tmp_gene_name+'_next'
        except:
            tmp_gene_name = input(
                "\n{0}: {1} Previous: {2}. Current: {3}.\nPlease input current gene name:".format(file_no, file_name, tmp_gene_name, ele.location.parts))
    else:
        tmp_gene_name = ele.qualifiers['gene'][0]

    tmp_list, gene_seq = merge_sequence(ele, complete_seq)
    pos_info_in_gene_note = " ["
    for i in range(2*len(ele.location.parts)):
        if i % 2 == 0:
            pos_info_in_gene_note += tmp_list[i]+'..'
        elif i % 2 == 1:
            pos_info_in_gene_note += tmp_list[i]+';'
    pos_info_in_gene_note = pos_info_in_gene_note.rstrip(';')+"]"
    gene_note = ">" + seq_id + pos_info_in_gene_note + " [gene=" + \
        tmp_gene_name + "]" + "\n"  # '>'后的格式和已有脚本兼容
    return gene_note, gene_seq, tmp_gene_name


def gene_name_standardization(gene_name):  # 格式化基因名字,可重复使用
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 13:
            if all_gene_list_lower[i] == gene_name:
                gene_name = all_gene_list_upper[i]
                break
            else:
                i += 1
        if i >= 13:
            print(gene_name)
            print('WARNING!Please check!')
    return gene_name


def get_gene(gbk_file_path, flag, dict_gene_len, file_no):  # 解析gbk文件获取cds
    """完整基因组"""
    file_name = os.path.basename(gbk_file_path)
    seq_record = SeqIO.read(gbk_file_path, "genbank")
    complete_seq = str(seq_record.seq)
    complete_note, seq_id = get_complete_note(seq_record)
    complete_fasta = format_fasta(complete_note, complete_seq, 80)  # 80换行本例不采用
    """gene序列"""
    cds_count = 0  # 对cds数量计数
    trna_count = 0
    cds_fasta = ""
    trna_fasta = ""
    list_cds_name = []  # 统计cds种类，记录重复
    list_trna_name = []  # 统计trna种类,记录重复
    tmp_gene_name = ''  # 上一个基因名字,为子函数get_gene_note()准备的
    for ele in seq_record.features:
        '''
        CDS
        tRNA
        rRNA
        repeat_region
        '''
        if ele.type == "CDS":
            cds_count += 1
            cds_note, cds_seq, tmp_gene_name = get_gene_note(
                file_no, file_name, ele, complete_seq, seq_id, tmp_gene_name)
            # list_cds_name.append(tmp_gene_name)  # 本次的基因名字 复用
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
            gene_name = tmp_gene_name
            gene_name = gene_name  # gene_name_standardization(gene_name)
            list_cds_name.append(gene_name)  # 存入列表
            # dict_gene_len[gene_name].append(3*(len(ele.qualifiers['translation'][0])+1))  # cds序列长度
            if (flag):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
        elif ele.type == 'tRNA':
            trna_count += 1
            trna_note, trna_seq, tmp_gene_name = get_gene_note(
                file_no, file_name, ele, complete_seq, seq_id, tmp_gene_name)
            trna_fasta += format_fasta(trna_note, trna_seq, 70)
            gene_name = tmp_gene_name
            gene_name = gene_name  # gene_name_standardization(gene_name)
            list_trna_name.append(gene_name)  # 存入列表

    s = '{2}: {0} has {1} CDS'.format(file_name, cds_count, file_no)
    if cds_count == 0:
        # --------There may be no comments--------'.format(
        s = '{2}: {0} has {1} CDS'.format(file_name, cds_count, file_no)
        print(s.ljust(50), '----------There may be no comments----------')
    elif cds_count != 0:
        # 20220811 输出左对齐  str.ljust(50)  达到50个字符宽度
        print(s.ljust(50), '+', 'list_cds_name')

    return file_name, seq_id, complete_fasta, cds_fasta, cds_count, list_cds_name,  trna_fasta, trna_count, list_trna_name, dict_gene_len, s


def create_gene_by_gap(dict_missing_gene, dict_gene_len, cds_file_path):  # 用gap构造没有的基因
    for i in dict_missing_gene.keys():
        cds_fasta = ''
        for j in dict_missing_gene[i]:
            ave = round(sum(dict_gene_len[j]) /
                        len(dict_gene_len[j]))  # 该基因平均长度
            gene_note = (i+' [0..0]'+' [gene={}]').format(j)
            gene_seq = ave*'-'
            cds_fasta += format_fasta(gene_note, gene_seq, 70)
        print(cds_fasta)
        file_name = (i.split('_')[-2]+'_' +
                     i.split('_')[-1]+'.1').lstrip('>')
        with open(cds_file_path, 'ab+') as f_cds:
            f_cds.write(cds_fasta.encode())


if __name__ == '__main__':
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    begin_time = time.time()
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('Start Time : {}'.format(start_time))
    #################################################################
    print('\n')
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if (not os.path.exists(os.path.join(args.output, 'trna'))):
        os.makedirs(os.path.join(args.output, 'trna'))
    if (not os.path.exists(os.path.join(args.output, 'cds'))):
        os.makedirs(os.path.join(args.output, 'cds'))
    if (not os.path.exists(os.path.join(args.output, 'complete'))):
        os.makedirs(os.path.join(args.output, 'complete'))

    """统计初始化"""
    dict_missing_gene = {}  # 每个文件中缺失的基因统计,总 字典
    dict_gene_len = {}  # 统计每个基因在不同物种中的长度,取平均
    list_seq_id_de_duplication = []  # 去重的物种
    list_unique_accession = []  # 去重的登录号
    """初始化"""
    dict_file_cds_count = {}  # 每个文件中cds计数
    file_list = [x for x in os.listdir(
        args.input) if os.path.isfile(os.path.join(args.input, x))]
    file_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字
    """主程序"""
    file_no = 0
    for file in file_list:
        file_no += 1
        gbk_file_path = os.path.join(args.input, file)
        file_name, seq_id, complete_fasta, cds_fasta, cds_count, list_cds_name,  \
            trna_fasta, trna_count, list_trna_name, dict_gene_len, s = get_gene(
                gbk_file_path, False, dict_gene_len, file_no)
        dict_file_cds_count[file_name] = cds_count  # 每个文件中cds计数
        '''
        20221217对gbk文件去重
        '''
        seq_id_content = seq_id.split('_')
        if seq_id.find('NC_') > 0:
            species = '_'.join(seq_id_content[:-2])
            accession = '_'.join(seq_id_content[-2:])
        else:
            species = '_'.join(seq_id_content[:-1])
            accession = '_'.join(seq_id_content[-1:])
        if species not in list_seq_id_de_duplication:
            list_seq_id_de_duplication.append(species)
            list_unique_accession.append(accession)
            """写入文件"""
            with open((os.path.join(args.output, 'complete', seq_id+'.fasta')), 'wb') as f_complete, \
                open((os.path.join(args.output, 'cds', 'cds_'+file_name.rstrip('.gbk')+'.fasta')), 'wb') as f_cds,\
                    open((os.path.join(args.output, 'trna', 'trna_'+file_name.rstrip('.gbk')+'.fasta')), 'wb') as f_trna:
                f_complete.write(complete_fasta.encode())
                f_cds.write(cds_fasta.encode())
                f_trna.write(trna_fasta.encode())
    print('{} left after removing duplicates'.format(len(list_unique_accession)))
    with open((os.path.join(args.output, 'list_tre')), 'w') as list_handle:
        for i in list_unique_accession:
            list_handle.write(i+'\n')
    print('\n')
    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
