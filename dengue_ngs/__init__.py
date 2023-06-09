import os
import re
import subprocess as sp
import sys
__version__ = "0.0.3"

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