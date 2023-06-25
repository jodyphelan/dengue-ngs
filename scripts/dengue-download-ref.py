#! /usr/bin/env python3
import os
import subprocess as sp
from tqdm import tqdm
import sys
from dengue_ngs import run_cmd

patterns = {
    "dengue virus 1":"DENV1",
    "dengue virus 2":"DENV2",
    "dengue virus 3":"DENV3",
    "dengue virus 4":"DENV4",
    "dengue virus type 1":"DENV1",
    "dengue virus type 2":"DENV2",
    "dengue virus type 3":"DENV3",
    "dengue virus type 4":"DENV4",
    "dengue virus i":"DENV1",
    "dengue virus type i":"DENV1",

}

taxid = {
    "DENV1":11053,
    "DENV2":11060,
    "DENV3":11069,
    "DENV4":11070,
}

def stream_fasta(f):
    seq = ""
    header = ""
    for l in open(f):
        if l.startswith(">"):
            if seq!="":
                yield(header,seq,serotype)
            header = l.strip().split()[0][1:]
            serotype = None
            for pattern in patterns:
                if pattern in l.lower():
                    serotype = patterns[pattern]
            seq = ""
        else:
            seq+=l.strip()
    yield(header,seq,serotype)

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
id2tax = {}
for name,seq,serotype in tqdm(stream_fasta('ncbi_dataset/data/genomic.fna')):
    if len(seq)<9000:
        continue
    if serotype is not None:
        id2tax[name] = taxid[serotype]
        with open(ref_dir + name + ".fasta",'w') as O:
            O.write(">%s\n%s\n" % (name,seq))

with open("taxid.map",'w') as O:
    for name in id2tax:
        O.write("%s\t%s\n" % (name,id2tax[name]))

sys.stderr.write("Creating sourmash signature\n")
run_cmd("sourmash sketch dna --singleton ncbi_dataset/data/genomic.fna -o dengue.sig")

sys.stderr.write("Creating kmcp database\n")
run_cmd("kmcp compute --circular -k 31 -n 10 -l 150 -I ref -O refs.tmp --force")
run_cmd("kmcp index -f 0.001 -I refs.tmp/ -O refs.kmcp --force")
run_cmd("mv taxid.map refs.kmcp/")

run_cmd("wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
run_cmd("mkdir -p taxdump")
run_cmd("tar -zxvf taxdump.tar.gz -C taxdump/")