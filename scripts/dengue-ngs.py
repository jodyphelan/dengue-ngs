#! /usr/bin/env python3
import argparse
from uuid import uuid4
import sys
import csv
import os
import dengue_ngs as dngs
import json


def main(args):
    ref_dir = os.path.expanduser('~')+"/.dengue-ngs/ref/"
    ref_sourfile = os.path.expanduser('~')+"/.dengue-ngs/dengue.sig"
    if not os.path.exists(ref_dir) or not os.path.exists(ref_sourfile):
        sys.stderr.write("Reference files not found. Please run dengue-download-ref.py first.\n")
        sys.exit(1)

    runs = dngs.find_fastq_files(args.folder)
    run_script = "%s.sh" % uuid4()
    with open(run_script,"w") as O:
        for run in runs:
            O.write(f"dengue-pipeline.py --threads {args.threads_per_job} --read1 {run.r1} --read2 {run.r2} --prefix {run.prefix} > {run.prefix}.log 2>&1\n")
    
    dngs.run_cmd(f"cat {run_script} | parallel -j {args.jobs} --bar")
    
    results = []
    for run in runs:
        if os.path.isfile(f"{run.prefix}.json"):
            results.append(json.load(open(f"{run.prefix}.json")))
        else:
            results.append({"Sample ID":run.prefix,"Median depth":"NA","Reference_coverage":"NA"})
    
    with open("run_results.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--folder',type=str,help='File with samples',required = True)
parser.add_argument('-t','--threads-per-job',type=int,help='File with samples',default=10)
parser.add_argument('-j','--jobs',type=int,help='File with samples',default=10)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)