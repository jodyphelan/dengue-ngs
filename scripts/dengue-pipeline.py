#! /usr/bin/env python3
import argparse
import subprocess as sp
import sys
from uuid import uuid4
import csv
import os
from dengue_ngs import run_cmd

def run_cmd(cmd):
    sys.stderr.write(f"Running: {cmd}\n")
    sp.call(cmd,shell=True)

def main(args):
    run_cmd("sourmash sketch dna %(read1)s  %(read2)s --merge  %(prefix)s  -o %(prefix)s.sig" % vars(args))
    run_cmd("sourmash gather --threshold-bp 1000 %(prefix)s.sig %(db)s -o %(prefix)s.gather.csv" % vars(args))
    if not os.path.isfile(f"{args.prefix}.gather.csv"):
        quit()
    rows = [row for row in csv.DictReader(open(f"{args.prefix}.gather.csv"))]
    args.ref = rows[0]["name"].split(" ")[0]+".fasta"

    if not os.path.isfile(f"{args.refdir}/{args.ref}.bwt"):
        run_cmd("bwa index %(refdir)s/%(ref)s" % vars(args))
    run_cmd("bwa mem -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina' %(refdir)s/%(ref)s %(read1)s %(read2)s | samtools sort -o %(prefix)s.bam" % vars(args))
    run_cmd("samtools index %(prefix)s.bam" % vars(args))
    run_cmd("pilon --genome %(refdir)s/%(ref)s --frags %(prefix)s.bam --output %(prefix)s.temp" % vars(args))
    run_cmd("sed 's/>/>%(prefix)s /' %(prefix)s.temp.fasta > %(prefix)s.fasta" % vars(args))
    os.remove(f"{args.prefix}.temp.fasta")
    
    run_cmd("bwa index %(prefix)s.fasta" % vars(args))
    run_cmd("bwa mem -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina' %(prefix)s.fasta %(read1)s %(read2)s | samtools sort -o %(prefix)s.consensus.bam" % vars(args))
    run_cmd("samtools index %(prefix)s.consensus.bam" % vars(args))
    run_cmd("lofreq call -f %(prefix)s.fasta -o %(prefix)s.lofreq.vcf %(prefix)s.consensus.bam" % vars(args))

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d','--db',type=str,help='File with samples',required = True)
parser.add_argument('-r','--refdir',type=str,help='File with samples',required = True)
parser.add_argument('-1','--read1',type=str,help='File with samples',required = True)
parser.add_argument('-2','--read2',type=str,help='File with samples',required = True)
parser.add_argument('-p','--prefix',type=str,help='File with samples',required = True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)