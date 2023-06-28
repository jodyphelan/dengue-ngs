#! /usr/bin/env python3
import argparse
from uuid import uuid4
import csv
import os
from dengue_ngs import *
import pathogenprofiler as pp
import shutil

def main(args):
    tmp = str(uuid4())
    args.data_dir = os.path.expanduser('~')+"/.dengue-ngs/"
    args.refdir = "%(data_dir)s/ref/" % vars(args)
    if not args.kraken_db:
        args.kraken_db = "%(refdir)s/kraken2" % vars(args)

    report = Report(args.prefix+".json")
    report.set("Sample ID",args.prefix)

    report.set_dict(get_fastq_stats([args.read1,args.read2]))

    run_cmd("kraken2 --db %(kraken_db)s --memory-mapping --report %(prefix)s.kreport.txt --output %(prefix)s.koutput.txt --threads 20 %(read1)s %(read2)s" % vars(args))

    report.set_dict(kreport_extract_human(f"{args.prefix}.kreport.txt"))
    report.set_dict(kreport_extract_dengue(f"{args.prefix}.kreport.txt"))

    run_cmd("extract_kraken_reads.py --include-children --fastq-output -t 12637  -k %(prefix)s.koutput.txt -r %(prefix)s.kreport.txt -s %(read1)s -s2 %(read2)s -o %(prefix)s.kraken_filtered.1.fq -o2 %(prefix)s.kraken_filtered.2.fq > /dev/null" % vars(args))
    run_cmd("pigz --force -p %(threads)s %(prefix)s.kraken_filtered.1.fq %(prefix)s.kraken_filtered.2.fq" % vars(args))

    args.read1 = f"{args.prefix}.kraken_filtered.1.fq.gz"
    args.read2 = f"{args.prefix}.kraken_filtered.2.fq.gz"
    
    run_cmd(f"megahit -t {args.threads} -1 {args.read1} -2 {args.read2} -o {args.prefix}.megahit -t {args.threads}")
    run_cmd(f"mv {args.prefix}.megahit/final.contigs.fa {args.prefix}.assembly.fasta")
    consensus_by_assembly = filter_seqs_by_size(
        fasta_file = f"{args.prefix}.assembly.fasta",
        output = f"{args.prefix}.temp.fasta",
        seqname = args.prefix,
        seq_size_cutoff=10500
    )
    shutil.rmtree(f"{args.prefix}.megahit")


    if not consensus_by_assembly:

        if args.reference_assignment_method=="sourmash":
            args.db = "%(data_dir)s/dengue.sig" % vars(args)
            run_cmd("sourmash sketch dna %(read1)s  %(read2)s --merge  %(prefix)s  -o %(prefix)s.sig" % vars(args))
            run_cmd("sourmash gather --threshold-bp 1000 %(prefix)s.sig %(db)s -o %(prefix)s.gather.csv" % vars(args))
            if not os.path.isfile(f"{args.prefix}.gather.csv"):
                quit()
            rows = [row for row in csv.DictReader(open(f"{args.prefix}.gather.csv"))]
            if len(rows)==0:
                quit("Can't find reference\n")
            args.ref = rows[0]["name"].split(" ")[0]+".fasta"
        else:
            args.db = os.path.expanduser('~')+"/.dengue-ngs/refs.kmcp/"
            run_cmd("kmcp search -j %(threads)s -d %(db)s %(read1)s  %(read2)s -o %(prefix)s.kmcp.tsv.gz" % vars(args))
            run_cmd("kmcp profile -X %(data_dir)s/taxdump/ -T %(db)s/taxid.map -m 1 %(prefix)s.kmcp.tsv.gz -o %(prefix)s.k.profile" % vars(args))
            if not os.path.isfile(f"{args.prefix}.k.profile"):
                quit()
            rows = [row for row in csv.DictReader(open(f"{args.prefix}.k.profile"),delimiter="\t")]
            if len(rows)==0:
                quit("Can't find reference\n")
            args.ref = rows[0]["ref"]+".fasta"

        pilon_correct(
            ref=f"{args.refdir}/{args.ref}",
            r1=args.read1,
            r2=args.read2,
            prefix=args.prefix+".temp",
            threads=args.threads,
        )
        
        fasta_depth_mask(
            input=f"{args.prefix}.temp.fasta",
            output=f"{args.prefix}.consensus.fasta",
            bam_file=f"{args.prefix}.bam",
            depth_cutoff=args.min_dp,
            newchrom=args.prefix
        )
        os.remove(f"{args.prefix}.temp.fasta")

        report.set('Consensus type','mapping')
        ## end of consensus by reference
    else:
        pilon_correct(
            ref=f"{args.prefix}.temp.fasta",
            r1=args.read1,
            r2=args.read2,
            prefix=args.prefix+".consensus",
            threads=args.threads,
        )
        report.set('Consensus type','assembly')


    run_cmd("bwa index %(prefix)s.consensus.fasta" % vars(args))
    run_cmd("bwa mem -t %(threads)s -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina' %(prefix)s.consensus.fasta %(read1)s %(read2)s | samtools sort -@ %(threads)s -o %(prefix)s.consensus.bam" % vars(args))
    run_cmd("samtools index %(prefix)s.consensus.bam" % vars(args))
    run_cmd("bedtools genomecov -d -ibam %(prefix)s.consensus.bam > %(prefix)s.consensus.depth.txt" % vars(args))
    run_cmd("lofreq call -f %(prefix)s.consensus.fasta -o %(prefix)s.lofreq.vcf %(prefix)s.consensus.bam" % vars(args))
    run_cmd(r"bcftools query -f '%POS\t%REF\t%ALT\t%AF\t%DP\t%QUAL\n' " + f"{args.prefix}.lofreq.vcf > {args.prefix}.lofreq.tsv")
    remove_bwa_index(f"{args.prefix}.consensus.fasta")
    bam = pp.bam(
        bam_file=args.prefix+".consensus.bam",
        prefix=args.prefix,
        platform="Illumina"
    )
    
    plot_lofreq_results(args.prefix)

    report.set("Median depth",bam.get_median_coverage(f"{args.prefix}.fasta"))
    report.set("Reference_coverage", 100 - get_fasta_missing_content(f"{args.prefix}.consensus.fasta"))   



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-1','--read1',type=str,help='Forward read',required = True)
parser.add_argument('-2','--read2',type=str,help='Reverse read',required = True)
parser.add_argument('-p','--prefix',type=str,help='Prefix for output files',required = True)
parser.add_argument('-t','--threads',type=int,help='Number of threads',default=4)
parser.add_argument('--min-dp',type=int,default=50,help='Minimum depth for consensus')
parser.add_argument('--reference-assignment-method',type=str,choices=['kmcp','sourmash'],default='kmcp',help='Minimum depth for consensus')
parser.add_argument('--kraken-db',type=str,help='Kraken2 database directory')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)