import os
import re
import subprocess as sp
import sys
import pathogenprofiler as pp
from uuid import uuid4
import json
import csv
import statistics as stats

__version__ = "0.0.6"

class Report:
    def __init__(self,report_file):
        self.report = {}
        self.report_file = report_file
    def set(self,key,value):
        self.report[key] = value
        self.dump()
    def set_dict(self,d):
        for key,value in d.items():
            self.set(key,value)
        self.dump()
    def dump(self):
        with open(self.report_file,"w") as O:
            json.dump(self.report,O,indent=4)


def get_fastq_stats(fastq_files):
    tmpfile = "%s.txt" % uuid4()
    run_cmd(f"seqkit stats -T {' '.join(fastq_files)} > {tmpfile}")
    numreads = 0
    lengths = []
    for row in csv.DictReader(open(tmpfile),delimiter="\t"):
        numreads  += int(row['num_seqs'])
        lengths.append(float(row['avg_len']))
    
    os.remove(tmpfile)

    return {
        "Number of reads":numreads,
        "Average read length": stats.mean(lengths),
    }

def kreport_extract_human(kreport_file):
    """
    Extract the human reads from a kraken report file.
    """
    human_reads = 0
    for line in open(kreport_file):
        if line.strip().split("\t")[4] == "9606":
            human_reads = float(line.strip().split("\t")[0])
    return {"Read percent human":human_reads}

def kreport_extract_dengue(kreport_file):
    """
    Extract the Dengue reads from a kraken report file.
    """
    dengue_reads = 0
    d1_reads = 0
    d2_reads = 0
    d3_reads = 0
    d4_reads = 0
    for line in open(kreport_file):
        if line.strip().split("\t")[4] == "12637":
            dengue_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11053":
            d1_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11060":
            d2_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11069":
            d3_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11070":
            d4_reads = float(line.strip().split("\t")[0])
    return {
        "Read percent dengue":dengue_reads, 
        "Read percent dengue 1":d1_reads, 
        "Read percent dengue 2":d2_reads, 
        "Read percent dengue 3":d3_reads, 
        "Read percent dengue 4":d4_reads
    }

def run_cmd(cmd):
    sys.stderr.write(f"Running: {cmd}\n")
    sp.call(cmd,shell=True)

class Sample:
    def __init__(self,prefix,r1,r2):
        self.prefix = prefix
        self.r1 = r1
        self.r2 = r2

    def __repr__(self):
        return f"Sample(prefix={self.prefix},r1={self.r1},r2={self.r2})"


def sort_out_paried_files(files,r1_suffix="_L001_R1_001.fastq.gz",r2_suffix="_L001_R2_001.fastq.gz"):
    prefixes = set()
    for f in files:
        tmp1 = re.search("(.+)%s" % r1_suffix,f)
        tmp2 = re.search("(.+)%s" % r2_suffix,f)
        if tmp1:
            prefixes.add(tmp1.group(1))
        if tmp2:
            prefixes.add(tmp2.group(1))
    runs = []
    for p in prefixes:
        r1 = p + r1_suffix
        r2 = p + r2_suffix
        if r1 not in files:
            raise Exception("%s is present in data file but not %s. Please check." % (r2,r1))
        if r2 not in files:
            raise Exception("%s is present in data file but not %s. Please check." % (r1,r2))
        runs.append(
            Sample(p.split("/")[-1],r1,r2))
    return runs

def find_fastq_files(directory):
    """
    Find fastq files in a directory and return a 
    list of tuples with the sample name and the 
    path to the fastq files from both pairs.
    """
    files = [f"{os.path.abspath(directory)}/{f}" for f in os.listdir(directory)]
    fastq_files = sort_out_paried_files(files)
    
    return fastq_files


def get_fasta_missing_content(fasta_file):
    """
    Get the missing content of a fasta file.
    """
    fasta = pp.fasta(fasta_file)
    seq = list(fasta.fa_dict.values())[0]
    missing = round(seq.count("N")/len(seq)*100)
    return missing

def mask_fasta(input,output,positions,newchrom=None):
    """
    Mask the fasta file with Ns.
    """
    fasta = pp.fasta(input)
    seq = list(list(fasta.fa_dict.values())[0])

    for chrom,pos in positions:
        seq[pos-1] = "N"
    
    with open(output,"w") as O:
        O.write(">%s\n%s\n" % (newchrom,''.join(seq)))

def get_missing_positions(bam,depth_cutoff=50):
    tmpfile = "%s.bed" % uuid4()
    run_cmd(f"bedtools genomecov -ibam {bam} -d > {tmpfile}")
    positions = []
    for line in open(tmpfile):
        chrom,pos,depth = line.strip().split("\t")
        if int(depth) < depth_cutoff:
            positions.append((chrom,int(pos)))
    os.remove(tmpfile)
    return positions

def fasta_depth_mask(input,output,bam_file,depth_cutoff=50,newchrom=None):
    """
    Mask the fasta file with Ns based on the depth of the bam file.
    """
    positions = get_missing_positions(bam_file,depth_cutoff=depth_cutoff)
    mask_fasta(input,output,positions,newchrom=newchrom)


def return_seqs_by_size(fasta_file,seq_size_cutoff=1000):
    """
    Return a list of sequences from a fasta file that are 
    above a certain size cutoff.
    """
    fasta = pp.fasta(fasta_file)
    seqs = {}
    for name,seq in fasta.fa_dict:
        if len(seq) > seq_size_cutoff:
            seqs[name] = seq
    return seqs

def filter_seqs_by_size(fasta_file,output, seq_size_cutoff=1000):
    """
    Filter a fasta file by a size cutoff.
    """
    seqs = return_seqs_by_size(fasta_file,seq_size_cutoff=seq_size_cutoff)
    if len(seqs) == 1:    
        with open(output,"w") as O:
            for name,seq in seqs.items():
                O.write(">%s\n%s\n" % (name,seq))
        return True
    else:
        return False
