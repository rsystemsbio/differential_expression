#!/usr/bin/env python3
#Rachel Calder 2022


#[1] Import Desired Packages
import argparse
import os
from os import walk

#[2] Establish commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, required=True, help='input directory with files')
parser.add_argument('-o', type=str, required=True, help='output directory')
args = parser.parse_args()



#[3] Acquire a list of all desired files
f = []
for (dirpath, dirnames, filenames) in walk(args.i):
    f.extend(filenames)
    break
f.sort()
z = 0

#index desired pairs
p = f[:]
print('Your files:', p)

#[4] Implementing Pre-FastQC, Trimmomatic, And Post FastQC of all files

while z != len(p):
    a = p[z]
    z += 1
    b = p[z]
    z += 1
   # print('------Initializing FastQC-------')
    os.system('fastqc ' + args.i + a + ' -o ' + args.o)
    os.system('fastqc ' + args.i + b + ' -o ' + args.o)
    print('-------Trimming Files--------')
    os.system('trimmomatic PE -phred33 -threads 20 ' + args.i + a + ' ' + args.i + b + ' ' +args.o+ a.split('.')[0] + '.trimmed.fastq ' + args.o + a.split('.')[0] + 'un.trimmed.fastq ' + 
    args.o+ b.split('.')[0] + '.trimmed.fastq ' + args.o+ b.split('.')[0] + 'un.trimmed.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36')
   # print('----------Initializing FastQC: Trimmed Files-------')
    os.system('fastqc ' + a.split('.')[0] + '.trimmed.fastq -o ' + args.o)
    os.system('fastqc ' + b.split('.')[0] + '.trimmed.fastq -o ' + args.o)

