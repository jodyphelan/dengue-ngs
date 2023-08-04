#! /usr/bin/env python
import argparse
import dengue_ngs as dn
import os

"""
This script reads in a alignment file and for each sequence in the alignment, it will
align the reads to the sequence and run lofreq to call variants
"""

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

parser = argparse.ArgumentParser(description='Run lofreq on aligned consensus',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a','--alignment',type=str,required=True,help='prefix for output files')
parser.add_argument('-t','--threads',type=int,default=1,help='Number of threads')

args = parser.parse_args()


for name,seq in stream_fasta(args.alignment):
    ref = "%s.%s.fa" % (args.alignment,name)
    with open(ref,"w") as f:
        f.write(">%s\n%s\n" % (name,seq))
    dn.run_cmd(f"bwa index {ref}")
    dn.run_cmd(f"bwa mem -t {args.threads} {ref} {name}.kraken_filtered.1.fq.gz {name}.kraken_filtered.2.fq.gz | samtools sort -@ {args.threads} -o {args.alignment}.{name}.bam -")
    dn.run_cmd(f"samtools index {args.alignment}.{name}.bam")
    dn.run_cmd(f"samtools faidx {ref}")
    if os.path.exists(f"{name}.lofreq.vcf"):
        os.remove(f"{name}.lofreq.vcf")
    dn.run_cmd(f"lofreq call-parallel --pp-threads {args.threads} --force-overwrite --ref {ref} -o {name}.lofreq.vcf {args.alignment}.{name}.bam")
    dn.run_cmd(r"bcftools query -f '%POS\t%REF\t%ALT\t%AF\t%DP\t%QUAL\n' " + f"{name}.lofreq.vcf > {name}.lofreq.tsv")

