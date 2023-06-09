#! /usr/bin/env python3
import argparse
import subprocess as sp
import sys
import csv
import os
import dengue_ngs as dngs
import pathogenprofiler as pp


def main(args):
    ref_dir = os.path.expanduser('~')+"/.dengue-ngs/ref/"
    ref_sourfile = os.path.expanduser('~')+"/.dengue-ngs/dengue.sig"
    print(ref_sourfile)
    if not os.path.exists(ref_dir) or not os.path.exists(ref_sourfile):
        sys.stderr.write("Reference files not found. Please run dengue-download-ref.py first.\n")
        sys.exit(1)

    runs = dngs.find_fastq_files(args.folder)
    results = []
    for run in runs:
        try:
            dngs.run_cmd(f"dengue-pipeline.py --threads {args.threads} --db {ref_sourfile} --refdir {ref_dir} --read1 {run.r1} --read2 {run.r2} --prefix {run.prefix}")
            bam = pp.bam(
                bam_file=run.prefix+".consensus.bam",
                prefix=run.prefix,
                platform="Illumina"
            )
            results.append({
                "Sample ID":run.prefix,
                "median_dp":bam.get_median_coverage(f"{run.prefix}.fasta")
            })
        except:
            pass
    with open("run_results.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--folder',type=str,help='File with samples',required = True)
parser.add_argument('-t','--threads',type=int,help='File with samples',default=4)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)