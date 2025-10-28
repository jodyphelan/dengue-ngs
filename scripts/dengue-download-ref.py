#! /usr/bin/env python3
import os
import subprocess as sp
from tqdm import tqdm
import sys
from dengue_ngs import run_cmd
import multiprocessing as mp
import argparse
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(sys.stderr)])

parser = argparse.ArgumentParser()
parser.add_argument("--threads",default=mp.cpu_count()//4,type=int)
parser.add_argument("--no-kraken",action="store_true")
parser.add_argument("--add-human",action="store_true")
parser.add_argument("--kmcp",action="store_true")
parser.add_argument("--db-dir",default=os.path.expanduser("~")+"/.dengue-ngs",type=str,help="Directory to store the database files")

args = parser.parse_args()



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
    "chikungunya virus":"CHIKV",
    "zika virus":"Zika virus",
    "west nile virus":"West Nile virus",

}

taxid = {
    "DENV1":11053,
    "DENV2":11060,
    "DENV3":11069,
    "DENV4":11070,
    "CHIKV":37124,
    "Zika virus":64320,
    "West Nile virus":11082
}

def get_exclude_list():
    exclude_list = []
    for l in open("%s/share/dengue-ngs/sample_exclusion.txt" % sys.base_prefix):
        exclude_list.append(l.strip())
    return(exclude_list)

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

data_dir = args.db_dir #os.path.expanduser('~')+"/.dengue-ngs/"
logging.info(f"Using data directory: {data_dir}")
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
os.chdir(data_dir)
ref_dir = data_dir + "/ref/"
kraken_ref_dir = data_dir + "/kraken_ref/"
if not os.path.exists(ref_dir):
    os.makedirs(ref_dir)
if not os.path.exists(kraken_ref_dir):
    os.makedirs(kraken_ref_dir)

exclude_list = get_exclude_list()


def process_fasta_file(filename,taxdict):
    count = 0
    id2tax = {}
    for name,seq,serotype in stream_fasta(filename):
        if name in exclude_list:
            continue
        if len(seq)<9000:
            continue
        if serotype is not None:
            id2tax[name] = taxdict[serotype]
            with open(kraken_ref_dir + name + ".fasta",'w') as O:
                O.write(f">{name}|kraken:taxid|{taxid[serotype]}\n{seq}\n")
            with open(ref_dir + name + ".fasta",'w') as O:
                O.write(f">{name}\n{seq}\n")
            count += 1
    logging.info(f"Found {count} sequences in {filename}")
    return id2tax

sys.stderr.write("Downloading reference files\n")

id2tax = {}

## downloading the dengue reference genomes from NCBI
logging.info("Downloading dengue virus genomes from NCBI")
run_cmd("datasets download virus genome taxon 12637 --complete-only")
run_cmd("unzip -o ncbi_dataset.zip")
sys.stderr.write("Processing reference files\n")
id2tax.update(process_fasta_file('ncbi_dataset/data/genomic.fna', taxid))

## adding an additional reference genome
logging.info("Adding additional dengue virus genome JQ045626.1")
run_cmd("datasets download virus genome accession JQ045626.1")
run_cmd("unzip -o ncbi_dataset.zip")
id2tax.update(process_fasta_file('ncbi_dataset/data/genomic.fna', taxid))

## downloading the chikungunya genomes
logging.info("Downloading chikungunya virus genomes from NCBI")
run_cmd("datasets download virus genome taxon 37124 --complete-only")
run_cmd("unzip -o ncbi_dataset.zip")
id2tax.update(process_fasta_file('ncbi_dataset/data/genomic.fna', taxid))

## downloading the zika virus genomes
logging.info("Downloading zika virus genomes from NCBI")
run_cmd("datasets download virus genome taxon 64320 --complete-only")
run_cmd("unzip -o ncbi_dataset.zip")
id2tax.update(process_fasta_file('ncbi_dataset/data/genomic.fna', taxid))

## downloadding WNV
logging.info("Downloading West Nile virus genomes from NCBI")
run_cmd("datasets download virus genome taxon 11082 --complete-only")
run_cmd("unzip -o ncbi_dataset.zip")
id2tax.update(process_fasta_file('ncbi_dataset/data/genomic.fna', taxid))

with open("taxid.map",'w') as O:
    for name in id2tax:
        O.write("%s\t%s\n" % (name,id2tax[name]))

sys.stderr.write("Creating sourmash signature\n")
run_cmd(f"sourmash sketch dna -p scaled=10 {data_dir}/ref/*.fasta --name-from-first -o sourmash.sig")

if args.kmcp:
    sys.stderr.write("Creating kmcp database\n")
    run_cmd(f"kmcp compute -j {args.threads} --circular -k 31 -n 10 -l 150 -I ref -O refs.tmp --force")
    run_cmd(f"kmcp index -j {args.threads} -f 0.01 -I refs.tmp/ -O refs.kmcp --force")
    run_cmd("cp taxid.map refs.kmcp/")

run_cmd("wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
run_cmd("mkdir -p taxdump")
run_cmd("tar -zxvf taxdump.tar.gz -C taxdump/")

if not args.no_kraken:
    sys.stderr.write("Creating kraken2 database\n")
    run_cmd(f"kraken2-build --threads {args.threads}  --download-taxonomy --skip-maps --db kraken2 --use-ftp")
    if args.add_human:
        logging.info("Adding human genome to kraken2 database")
        run_cmd(f"kraken2-build --threads {args.threads} --download-library human --db kraken2 --use-ftp")
    run_cmd("ls %s/ | parallel -j %s --bar kraken2-build --add-to-library %s/{} --db kraken2" % (kraken_ref_dir,args.threads,kraken_ref_dir))
    run_cmd(f"kraken2-build --build --db kraken2 --threads {args.threads}")
    run_cmd("kraken2-build --clean --db kraken2")

print("\nDone!\n")