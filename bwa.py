#!/usr/bin/env python3
#Rachel Calder 2022


#[1] Import Desired Packages
import argparse
import os
from os import walk

#[2] Establish commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, required=True, help='input directory with files')
parser.add_argument('-r', type=str, required=True, help='index file')
parser.add_argument('-o', type=str, required=False, help='output directory')
args = parser.parse_args()

#[3] Acquire a list of all desired files
f = []
for (dirpath, dirnames, filenames) in walk(args.i):
    f.extend(filenames)
    break
f.sort()
z = 0

p = f[:]
print('Your files:', p)

#[4] Implementing Mapping (BAM OUT)
while z != len(p):
    a = p[z]
    z += 1
    b = p[z]
    z += 1
    print('------Initializing Mapping for '+ a + '  -------')
    name = a.split('.')[0]
    name1 = b.split('.')[0]
    # aligns | converts to bam
    os.system('bwa mem -t 20 ' + args.r + ' ' + args.i + a +  ' ' + args.i + b + ' |/usr/local/pace-apps/manual/packages/samtools/1.14/intel-19.0.5/bin/samtools view -b > ' + name + '.bam') 