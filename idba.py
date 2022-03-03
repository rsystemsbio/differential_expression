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

p = f[:]
print('Your files:', p)

#[4] IDBA Assembly
while z != len(p):
    a = p[z]
    z += 1
    print('------Initializing Assembly or for '+ a+ ' and '+b+ '  -------')
    name = a.split('.')[0]
    #Make a directory for the out
    os.system('mkdir ' + args.o + '/' + name)
    #Assemble the reads 
    os.system('idba_ud -r ' +name+'.fa --num_threads 20 -o ' + args.o + '/' + name)
