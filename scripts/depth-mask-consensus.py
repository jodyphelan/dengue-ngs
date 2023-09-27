#! /usr/bin/env python
import argparse
import dengue_ngs as dn

parser = argparse.ArgumentParser(description='Mask low depth regions on consensus using bam file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',type=str,required=True,help='Input fasta file')
parser.add_argument('--bam',type=str,required=True,help='Input bam file')
parser.add_argument('--output',type=str,required=True,help='Output fasta file')
parser.add_argument('--depth-cutoff',type=int,default=50,help='Depth cutoff')

args = parser.parse_args()


dn.fasta_depth_mask(
    input=args.fasta,
    output=args.output,
    bam_file=args.bam,
    depth_cutoff=args.depth_cutoff
)