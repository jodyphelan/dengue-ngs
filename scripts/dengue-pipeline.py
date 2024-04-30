#! /usr/bin/env python3
import argparse
from uuid import uuid4
import csv
import os
from dengue_ngs import *
import pathogenprofiler as pp
import shutil
import yaml
import logging
from rich.logging import RichHandler

def main(args):
    tmp = str(uuid4())
    args.data_dir = os.path.expanduser('~')+"/.dengue-ngs/"
    args.refdir = "%(data_dir)s/ref/" % vars(args)
    if not args.kraken_db:
        args.kraken_db = f"{args.data_dir}/kraken2" 

    report = Report(args.prefix+".json")
    report.set("Sample ID",args.prefix)
    report.set("Analysis completed","no")

    report.set_dict(get_fastq_stats(args.read1,args.read2))
    
    if args.read2:
        run_cmd("kraken2 --db %(kraken_db)s --report %(prefix)s.kreport.txt --output %(prefix)s.koutput.txt --paired --threads %(threads)s %(read1)s %(read2)s" % vars(args))
    else:
        run_cmd("kraken2 --db %(kraken_db)s --report %(prefix)s.kreport.txt --output %(prefix)s.koutput.txt --threads %(threads)s %(read1)s" % vars(args))
    
    report.set_dict(kreport_extract_human(f"{args.prefix}.kreport.txt"))
    report.set_dict(kreport_extract_dengue(f"{args.prefix}.kreport.txt"))

    major_serotype = sorted([(x,report.get('Read percent dengue %s' % x)) for x in range(1,5)],key=lambda x:x[1],reverse=True)[0][0]
    
    
    kraken_output = f"{args.prefix}.koutput.txt"
    filter_fastq_by_taxon(kraken_output=kraken_output, serotype = major_serotype, reads=args.read1, output=f"{args.prefix}.kraken_filtered.1.fq")
    if args.read2:
        filter_fastq_by_taxon(kraken_output=kraken_output, serotype = major_serotype, reads=args.read2, output=f"{args.prefix}.kraken_filtered.2.fq")



    args.read1 = f"{args.prefix}.kraken_filtered.1.fq"
    if args.read2:
        args.read2 = f"{args.prefix}.kraken_filtered.2.fq"

    print("adjaiodjsa")
    tmp_stats = (get_fastq_stats(args.read1,args.read2))
    if tmp_stats['Number of reads']==0:
        logging.critical("No reads left after filtering, exiting")
        quit()
    if args.platform=="nanopore":
        args.assemble = False
        args.serotype_ref = True

    if args.assemble==True:
    
        run_cmd(f"megahit -t {args.threads} -1 {args.read1} -2 {args.read2} -o {args.prefix}.megahit -t {args.threads}")
        run_cmd(f"mv {args.prefix}.megahit/final.contigs.fa {args.prefix}.assembly.fasta")
        consensus_by_assembly = filter_seqs_by_size(
            fasta_file = f"{args.prefix}.assembly.fasta",
            output = f"{args.prefix}.temp.fasta",
            seqname = args.prefix,
            min_seq_size_cutoff=10000,
            max_seq_size_cutoff=12000
        )
        shutil.rmtree(f"{args.prefix}.megahit")
    else:
        consensus_by_assembly = False

    
    report.set("Serotype",major_serotype)
    if not consensus_by_assembly:
        if args.serotype_ref:
            serotype_references =  yaml.safe_load(open("%s/share/dengue-ngs/variables.yaml" % sys.base_prefix))['references']
            args.ref = serotype_references[major_serotype]+".fasta"
            args.ref = f"{args.refdir}/{args.ref}"
            report.set("Consensus type","Serotype reference mapping")
        elif args.fix_ref:
            args.ref = args.fix_ref
            report.set("Consensus type","Fixed reference mapping")
        else:

            if args.reference_assignment_method=="sourmash":
                args.db = "%(data_dir)s/dengue.sig" % vars(args)
                run_cmd("sourmash sketch dna -p abund,scaled=10 %(read1)s  %(read2)s --merge  %(prefix)s  -o %(prefix)s.sig" % vars(args))
                run_cmd("sourmash gather --threshold-bp 1000 %(prefix)s.sig %(db)s -o %(prefix)s.gather.csv" % vars(args))
                if not os.path.isfile(f"{args.prefix}.gather.csv"):
                    quit()
                rows = [row for row in csv.DictReader(open(f"{args.prefix}.gather.csv"))]
                if len(rows)==0:
                    quit("Can't find reference\n")
                args.ref = rows[0]["name"].split(" ")[0]+".fasta"
                args.ref = f"{args.refdir}/{args.ref}"
                report.set("Consensus type","Closest reference mapping")
                report.set("Reference used",args.ref)
                
            else:
                args.db = os.path.expanduser('~')+"/.dengue-ngs/refs.kmcp/"
                run_cmd("kmcp search -j %(threads)s -d %(db)s %(read1)s  %(read2)s -o %(prefix)s.kmcp.tsv.gz" % vars(args))
                run_cmd("kmcp profile -j %(threads)s -X %(data_dir)s/taxdump/ -T %(db)s/taxid.map -m 1 %(prefix)s.kmcp.tsv.gz -o %(prefix)s.k.profile" % vars(args))
                if not os.path.isfile(f"{args.prefix}.k.profile"):
                    quit()
                rows = [row for row in csv.DictReader(open(f"{args.prefix}.k.profile"),delimiter="\t")]
                if len(rows)==0:
                    quit("Can't find reference\n")
                args.ref = rows[0]["ref"]+".fasta"
                args.ref = f"{args.refdir}/{args.ref}"
                report.set("Consensus type","Closest reference mapping")
                report.set("Reference used",args.ref)

        pilon_correct(
            ref=args.ref,
            r1=args.read1,
            r2=args.read2,
            platform=args.platform,
            consensus_name=args.prefix+".temp.fasta",
            bam_file=f"{args.prefix}.ref.bam",
            threads=args.threads,
            min_depth=args.min_dp
        )
        
        fasta_depth_mask(
            input=f"{args.prefix}.temp.fasta",
            output=f"{args.prefix}.temp.consensus.fasta",
            bam_file=f"{args.prefix}.ref.bam",
            depth_cutoff=args.min_dp,
            newchrom=args.prefix
        )
        os.remove(f"{args.prefix}.temp.fasta")

        ## end of consensus by reference
    else:
        pilon_correct(
            ref=f"{args.prefix}.temp.fasta",
            r1=args.read1,
            r2=args.read2,
            consensus_name=args.prefix+".consensus.fasta",
            threads=args.threads,
        )
        report.set('Consensus type','Assembly')

        strand_direction = get_strand_direction(f"{args.prefix}.consensus.fasta",f"{args.refdir}/EU677137.1.fasta")
        if strand_direction=="-":
            tmpfile = str(uuid4())
            run_cmd(f"seqkit seq --reverse --complement -t DNA {args.prefix}.consensus.fasta > {tmpfile}")
            shutil.move(tmpfile,f"{args.prefix}.consensus.fasta")

    sys.stderr.write("Consensus sequence generated\n")
    if args.platform=="illumina":
        run_cmd("bwa index %(prefix)s.temp.consensus.fasta" % vars(args))
        run_cmd("bwa mem -t %(threads)s -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' %(prefix)s.temp.consensus.fasta %(read1)s %(read2)s | samtools sort -@ %(threads)s -o %(prefix)s.consensus.bam" % vars(args))
        remove_bwa_index(f"{args.prefix}.temp.consensus.fasta")
    else:
        run_cmd("minimap2 -ax map-ont -t %(threads)s %(prefix)s.temp.consensus.fasta %(read1)s | samtools sort -@ %(threads)s -o %(prefix)s.consensus.bam" % vars(args))
    run_cmd("samtools index %(prefix)s.consensus.bam" % vars(args))
    
    bam = pp.bam(
        bam_file=args.prefix+".consensus.bam",
        prefix=args.prefix,
        platform=args.platform.title(),
    )

    freebayes_correct(
        ref=f"{args.prefix}.temp.consensus.fasta",
        bam=args.prefix+".consensus.bam",
        platform=args.platform,
        prefix=args.prefix,
        output=f"{args.prefix}.final.consensus.fasta",
        threads=args.threads,
        min_depth=args.min_dp,
        min_freq=args.consensus_variant_frequency
    )

    report.set("Median depth",bam.get_median_coverage(f"{args.prefix}.fasta"))
    report.set("Reference_coverage", 100 - get_fasta_missing_content(f"{args.prefix}.final.consensus.fasta"))   

    report.set("Analysis completed","yes")


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-1','--read1',type=str,help='Forward read',required = True)
parser.add_argument('-2','--read2',type=str,help='Reverse read')
parser.add_argument('-p','--prefix',type=str,help='Prefix for output files',required = True)
parser.add_argument('--platform',type=str,choices=['illumina','nanopore'],help='Sequencing platform',required=True)
parser.add_argument('-t','--threads',type=int,help='Number of threads',default=4)
parser.add_argument('--min-dp',type=int,default=50,help='Minimum depth for consensus')
parser.add_argument('--reference-assignment-method',type=str,choices=['kmcp','sourmash'],default='sourmash',help='Minimum depth for consensus')
parser.add_argument('--kraken-db',type=str,help='Kraken2 database directory')
parser.add_argument('--serotype-ref',action="store_true",help='Use serotype reference instead of building one')
parser.add_argument('--fix-ref',help='Force a reference instead of building one')
parser.add_argument('--consensus-variant-frequency',default=0.01,type=float,help='Minimum frequency of variant to be considered in final reference')
parser.add_argument('--assemble',action="store_true",help='Try assembly')
parser.add_argument('--logging',default="INFO",choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],help='Logging level')
parser.set_defaults(func=main)

args = parser.parse_args()

logging.basicConfig(
    level=args.logging, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]
)

args.func(args)