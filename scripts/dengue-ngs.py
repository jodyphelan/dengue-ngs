#! /usr/bin/env python3
import argparse
from uuid import uuid4
import sys
import csv
import os
import dengue_ngs as dngs
import json
import multiprocessing as mp


threads_per_job = mp.cpu_count()//4

def main(args):
    data_dir = os.path.expanduser('~')+"/.dengue-ngs/"
    ref_dir = os.path.expanduser('~')+"/.dengue-ngs/ref/"
    if not args.kraken_db:
        args.kraken_db = f"{data_dir}/kraken2" 
    ref_sourfile = os.path.expanduser('~')+"/.dengue-ngs/dengue.sig"
    if not os.path.exists(ref_dir) or not os.path.exists(ref_sourfile):
        sys.stderr.write("Reference files not found. Please run dengue-download-ref.py first.\n")
        sys.exit(1)

    runs = dngs.find_fastq_files(args.folder)
    run_script = "%s.sh" % uuid4()
    with open(run_script,"w") as O:
        for run in runs:
            tmp_args = "--fix-ref" if args.fix_ref else ""
            O.write(f"dengue-pipeline.py {tmp_args} --kraken-db {args.kraken_db} --threads {args.threads_per_job} --read1 {run.r1} --read2 {run.r2} --prefix {run.prefix} > {run.prefix}.log 2>&1\n")
    
    if not args.collate:
        dngs.run_cmd(f"cat {run_script} | parallel -j {args.jobs}")
    
    results = []

    fieldnames = [
        "Sample ID","Number of reads","Average read length","Read percent human",
        "Read percent dengue","Read percent dengue 1","Read percent dengue 2",
        "Read percent dengue 3", "Read percent dengue 4","Consensus type",
        "Median depth","Reference_coverage"
    ]

    for run in runs:
        res = {}
        if os.path.isfile(f"{run.prefix}.json"):
            d = json.load(open(f"{run.prefix}.json"))
            for f in fieldnames:
                res[f] = d.get(f,"NA")
        else:
            res = {f:"NA" for f in fieldnames}
        results.append(res)
    
    with open("run_results.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--folder',type=str,help='File with samples',required = True)
parser.add_argument('-t','--threads-per-job',type=int,help='File with samples',default=threads_per_job)
parser.add_argument('-j','--jobs',type=int,help='File with samples',default=4)
parser.add_argument('-k','--kraken-db',type=str,help='Kraken2 database directory')
parser.add_argument('-c','--collate',action="store_true",help='Only collate existing results')
parser.add_argument('--fix-ref',action="store_true",help='Use serotype reference instead of building one')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)