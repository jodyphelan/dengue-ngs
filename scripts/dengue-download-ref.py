#! /usr/bin/env python3
import os
import subprocess as sp
from tqdm import tqdm
import sys
from dengue_wgs import run_cmd


def stream_fasta(f):
    seq = ""
    header = ""
    for l in open(f):
        if l.startswith(">"):
            if seq!="":
                yield(header,seq)
            header = l.strip().split()[0][1:]
            seq = ""
        else:
            seq+=l.strip()
    yield(header,seq)

data_dir = os.path.expanduser('~')+"/.dengue-ngs/"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
os.chdir(data_dir)
ref_dir = data_dir + "/ref/"
if not os.path.exists(ref_dir):
    os.makedirs(ref_dir)

run_cmd("datasets download virus genome taxon 12637 --complete-only")
run_cmd("unzip ncbi_dataset.zip")

sys.stderr.write("Processing reference files\n")
for name,seq in tqdm(stream_fasta('ncbi_dataset/data/genomic.fna')):
    with open(ref_dir + name + ".fasta",'w') as O:
        O.write(">%s\n%s\n" % (name,seq))

sys.stderr.write("Creating sourmash signature\n")
run_cmd("sourmash sketch dna --singleton ncbi_dataset/data/genomic.fna -o dengue.sig")