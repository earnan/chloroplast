from Bio import SeqIO
from Bio import Entrez
import os
import sys
import argparse
parser = argparse.ArgumentParser(
    description='This script was used to download gb or fasta file of cp genome from NCBI nucleotides database')
parser.add_argument(
    '-f', '--format', help='Please input file format you want to download', required=True)
parser.add_argument('-a', '--accession',
                    help='file name contain accession number of cp genome you want to download', required=True)
args = parser.parse_args()


fr = open(args.accession, 'r')
acc = []
for line in fr:
    acc.append(line.strip())
print("The number of sequence will be downloaded is: ", len(acc))

Entrez.email = ""  # 在这里添加邮箱

fold_name = "downladed_seq_"+args.format
os.mkdir(fold_name)
os.chdir(fold_name)

for line in acc:
    print(line + " will be downloaded!")
    fw = open(line+"."+args.format, 'w')
    cp = Entrez.efetch(db='nucleotide', id=[line], rettype=args.format)
    seq = SeqIO.read(cp, args.format)
    SeqIO.write(seq, fw, args.format)
    fw.close()
    print(line + " is downloaded!")

print(len(os.listdir()), "sequence have been downloaded")
