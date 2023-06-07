#! /usr/bin/env python3
import argparse
import subprocess as sp
import sys
import csv
import os
import dengue_ngs as dwgs



def main(args):
    ref_dir = os.path.expanduser('~')+"/.dengue-ngs/ref/"
    ref_sourfile = os.path.expanduser('~')+"/.dengue-ngs/dengue.sig"
    print(ref_sourfile)
    if not os.path.exists(ref_dir) or not os.path.exists(ref_sourfile):
        sys.stderr.write("Reference files not found. Please run dengue-download-ref.py first.\n")
        sys.exit(1)

    runs = dwgs.find_fastq_files(args.folder)
    for run in runs:
        dwgs.run_cmd(f"dengue-pipeline.py --db {ref_sourfile} --refdir {ref_dir} --read1 {run.r1} --read2 {run.r2} --prefix {run.prefix}")
    

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--folder',type=str,help='File with samples',required = True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)