import os
import re
import subprocess as sp
import sys
import pathogenprofiler as pp
from uuid import uuid4
import json
import csv
import statistics as stats
from collections import defaultdict
from glob import glob
import pandas as pd
import re
import logging
import pysam
from os.path import expanduser
import yaml

__version__ = "0.0.13"

conf = yaml.safe_load(open("%s/share/dengue-ngs/variables.yaml" % sys.base_prefix))
taxid2name = conf['taxid2name']
sero2tax = {1:11053, 2:11060,3:11069,4:11070}

def get_conf():
    """
    Get the configuration dictionary.
    """
    return conf

def get_strand_direction(consensus,ref):
    """
    Get the strand direction of the consensus sequence relative to a reference sequence
    """
    tmpfile = str(uuid4())
    run_cmd(f"minimap2 {ref} {consensus} > {tmpfile}.paf")
    strand = open(f"{tmpfile}.paf").readline().strip().split()[4]
    return strand 

def file_line_count(filename):
    """
    Count the number of lines in a file
    """
    return int(sp.check_output(f"wc -l {filename}",shell=True).decode().split()[0])

def plot_lofreq_results(prefix,lofreq_tsv,depth_file):
    # Read lofreq results
    if file_line_count(lofreq_tsv)==0:
        return
    df = pd.read_csv(lofreq_tsv,delimiter="\t",header = None)
    df.rename(columns={0: 'position', 1:'ref',2:'alt',3:'frequency',4:'depth',5:'qual'}, inplace=True)
    
    # Read depth
    dp = pd.read_csv(depth_file, delimiter="\t",header = None)
    dp.rename(columns={0: 'chrom', 1:'pos',2:'depth'}, inplace=True)

    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # Add traces colouring by depth setting the max and min of the colour scale
    fig.add_trace(
        go.Scatter(
            x=df.position, 
            y=df.frequency, 
            name="Lofreq SNP",
            mode='markers', 
            marker=dict(color=df.qual, colorscale='Viridis',showscale=True,cmin=65,cmax=200),
        ),
        secondary_y=False, 
    )

    # Add depth trace
    fig.add_trace(
        go.Scatter(x=dp.pos, y=dp.depth, name="Depth",mode='lines'),
        secondary_y=True
    )
    fig.update_layout(
        yaxis2=dict(type='log')
    ) # update

    # add legend on top left
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )

    )

    # use simple white template
    fig.update_layout(template='simple_white')


    # write to pdf
    fig.write_html(f"{prefix}.lofreq.html")


def remove_bwa_index(ref):
    for e in ["amb","ann","bwt","pac","sa"]:
        os.remove(f"{ref}.{e}")

def pilon_correct(ref,r1,r2,consensus_name,platform,sample_name,bam_file=None,threads=1,min_depth=50):
    tmp = str(uuid4())
    run_cmd(f"cp {ref} {tmp}.ref.fasta")
    if platform.lower()=="illumina":
        run_cmd(f"bwa index {tmp}.ref.fasta")
    
        run_cmd(f"bwa mem -t {threads} {tmp}.ref.fasta {r1} {r2} | samtools sort -@ {threads} -o {tmp}.bam")
        run_cmd(f"samtools index {tmp}.bam")
        run_cmd(f"pilon -Xmx10g --genome {tmp}.ref.fasta --frags {tmp}.bam --output {tmp}.consensus --mindepth {min_depth} --vcf")
    elif platform.lower()=="nanopore":
        run_cmd(f"minimap2 -ax map-ont {ref} {r1} | samtools sort -@ {threads} -o {tmp}.bam")
        run_cmd(f"samtools index {tmp}.bam")
        run_cmd(f"pilon -Xmx10g --genome {ref} --nanopore {tmp}.bam --output {tmp}.consensus --mindepth {min_depth} --vcf")
    
    run_cmd(f"bcftools view -v snps -i 'FILTER=\"PASS\" & QUAL>0'  -c 2 {tmp}.consensus.vcf -Oz -o {tmp}.variants.vcf.gz")
    run_cmd(f"tabix {tmp}.variants.vcf.gz")
    run_cmd(f"bcftools consensus -I -f {tmp}.ref.fasta {tmp}.variants.vcf.gz > {tmp}.consensus.fasta")
    
    # run_cmd(f"bcftools query -f '%POS\\t%FILTER\\n' {tmp}.consensus.vcf > {tmp}.consensus.info")
    # mask_positions = []
    # for line in open(f"{tmp}.consensus.info"):
    #     pos,filter = line.strip().split("\t")
    #     if filter!="PASS":
    #         mask_positions.append(('chromosome',int(pos)))

    # mask_fasta(f"{tmp}.consensus.fasta",consensus_name,mask_positions,newchrom=consensus_name)
    run_cmd(f"sed -r 's/>.+/>{sample_name}/' {tmp}.consensus.fasta > {consensus_name}")

    run_cmd(f"mv {tmp}.consensus.vcf {consensus_name}.vcf")
    if bam_file:
        run_cmd(f"mv {tmp}.bam {bam_file}")
        run_cmd(f"mv {tmp}.bam.bai {bam_file}.bai")
    for f in glob(f"{tmp}.*"):
        os.remove(f)

def medaka_correct(ref,reads,output, prefix=None,min_depth=50):
    logging.debug("Correcting consensus using medaka %s" % ref)
    tmpdir = str(uuid4())
    run_cmd(f"medaka_consensus -i {reads} -d {ref} -o {tmpdir}")
    run_cmd(f"minimap2 -ax map-ont {tmpdir}/consensus.fasta {reads} | samtools sort -o {tmpdir}/reads_to_consensus.bam")
    run_cmd(f"samtools index {tmpdir}/reads_to_consensus.bam")
    
    fasta_depth_mask(
        input=f"{tmpdir}/consensus.fasta",
        output=output,
        bam_file=f"{tmpdir}/reads_to_consensus.bam",
        depth_cutoff=min_depth,
        newchrom=prefix
    )

    # run_cmd(f"mv {tmp}.medaka/round_1.fasta {output}")
    run_cmd(f"rm -rf {tmpdir}")
    
def freebayes_correct(ref,output,platform,bam=None,r1=None,r2=None,prefix=None,threads=1,min_depth=50, min_freq=0.01, min_ad = 50):
    logging.debug("Correcting consensus using freebayes %s" % ref)
    tmp = str(uuid4())
    run_cmd(f"cp {ref} {tmp}.ref.fasta")
    if bam is None:
        if platform.lower()=="illumina":
            run_cmd(f"bwa index {tmp}.ref.fasta")
            run_cmd(f"bwa mem -t {threads} {tmp}.ref.fasta {r1} {r2} | samtools sort -@ {threads} -o {tmp}.bam")
        elif platform.lower()=="nanopore":
            run_cmd(f"minimap2 -t {threads} -ax map-ont {ref} {r1} | samtools sort -@ {threads} -o {tmp}.bam")
        bam = f"{tmp}.bam"
    run_cmd(f"samtools index {bam}")
    run_cmd(f"freebayes -f {tmp}.ref.fasta {bam} -F {min_freq} -C {min_ad} | bcftools norm -a | bcftools view -v snps -Oz -o {tmp}.variants.vcf.gz")
    
    run_cmd(f"tabix {tmp}.variants.vcf.gz")
    run_cmd(f"bcftools consensus -I -f {tmp}.ref.fasta {tmp}.variants.vcf.gz > {tmp}.consensus.fasta")
    
    # run_cmd(f"bcftools view -i 'GT=\"het\"' {tmp}.variants.vcf.gz | bcftools query -f '%POS\\n' > {tmp}.het.pos")
    # to_mask = []
    # for line in open(f"{tmp}.het.pos"):
        # to_mask.append(('chromosome',int(line.strip())))
    # het_masked_consensus = f"{tmp}.consensus.het_masked.fasta"
    # mask_fasta(f"{tmp}.consensus.fasta",het_masked_consensus,to_mask,newchrom=prefix)

    fasta_depth_mask(
        input=f"{tmp}.consensus.fasta",
        output=output,
        bam_file=bam,
        depth_cutoff=min_depth,
        newchrom=prefix
    )

    # run_cmd(f"mv {tmp}.consensus.fasta {output}")
    for f in glob(f"{tmp}.*"):
        os.remove(f)



class Report:
    def __init__(self,report_file):
        self.report = {}
        self.report_file = report_file
        self.tax2reads = {}
        for tax in taxid2name:
            self.tax2reads[tax] = 0
    def get(self,key):
        return self.report[key]
    def set(self,key,value):
        self.report[key] = value
        self.dump()
    def set_dict(self,d):
        for key,value in d.items():
            self.set(key,value)
        self.dump()
    def dump(self):
        with open(self.report_file,"w") as O:
            json.dump(self.report,O,indent=4)
    def load_kreport(self, kreport_file):
        """
        Load a kraken report file into the report dictionary.
        """
        for line in open(kreport_file):
            taxid = int(line.strip().split("\t")[4])
            if taxid in self.tax2reads:
                self.tax2reads[taxid] = float(line.strip().split("\t")[0])
    def get_major_virus(self):
        """
        Get the major virus from the report.
        """
        major_virus = max(self.tax2reads.items(), key=lambda x: x[1])

        if major_virus[1] > 0:
            return major_virus[0]
        else:
            logging.warning("No major virus found in the report.")

def get_fastq_stats(read1,read2=None):
    tmpfile = "%s.txt" % uuid4()
    if read2:
        run_cmd(f"seqkit stats -T {read1} {read2} > {tmpfile}")
    else:
        run_cmd(f"seqkit stats -T {read1} > {tmpfile}")
    numreads = 0
    lengths = []
    for row in csv.DictReader(open(tmpfile),delimiter="\t"):
        numreads  += int(row['num_seqs'])
        lengths.append(float(row['avg_len']))
    
    os.remove(tmpfile)

    return {
        "Number of reads":numreads,
        "Average read length": stats.mean(lengths),
    }

def kreport_extract_human(kreport_file):
    """
    Extract the human reads from a kraken report file.
    """
    human_reads = 0
    for line in open(kreport_file):
        if line.strip().split("\t")[4] == "9606":
            human_reads = float(line.strip().split("\t")[0])
    return {"Read percent human":human_reads}

def kreport_extract_dengue(kreport_file):
    """
    Extract the Dengue reads from a kraken report file.
    """
    dengue_reads = 0
    d1_reads = 0
    d2_reads = 0
    d3_reads = 0
    d4_reads = 0
    chik_reads = 0
    zika_reads = 0

    for line in open(kreport_file):
        if line.strip().split("\t")[4] == "12637":
            dengue_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11053":
            d1_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11060":
            d2_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11069":
            d3_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "11070":
            d4_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "37124":
            chik_reads = float(line.strip().split("\t")[0])
        if line.strip().split("\t")[4] == "64320":
            zika_reads = float(line.strip().split("\t")[0])
        
    return {
        "Read percent dengue":dengue_reads, 
        "Read percent dengue 1":d1_reads, 
        "Read percent dengue 2":d2_reads, 
        "Read percent dengue 3":d3_reads, 
        "Read percent dengue 4":d4_reads,
        "Read percent chikungunya":chik_reads,
        "Read percent zika":zika_reads
    }

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True
    return False

def run_cmd(cmd: str, desc=None, log: str=None) -> sp.CompletedProcess:
    cmd = "/bin/bash -c set -o pipefail; " + cmd
    logging.debug(f"Running command: {cmd}")
    output = open(log,"w") if log else sp.PIPE
    result = sp.run(cmd,shell=True,stderr=output,stdout=output)
    if result.returncode != 0:
        raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,result.stderr.decode()))
    return result



class Sample:
    def __init__(self,prefix,r1,r2):
        self.prefix = prefix
        self.r1 = r1
        self.r2 = r2

    def __repr__(self):
        return f"Sample(prefix={self.prefix},r1={self.r1},r2={self.r2})"


def sort_out_paried_files(files,r1_suffix="_S[0-9]+_L001_R1_001.fastq.gz",r2_suffix="_S[0-9]+_L001_R2_001.fastq.gz"):
    prefixes = defaultdict(lambda:{"r1":[],"r2":[]})

    for f in files:
        tmp1 = re.search("(.+)%s" % r1_suffix,f)
        tmp2 = re.search("(.+)%s" % r2_suffix,f)
        p = None
        if tmp1:
            p = tmp1.group(1).split("/")[-1]
            prefixes[p]['r1'].append(f)
        elif tmp2:
            p = tmp2.group(1).split("/")[-1]
            prefixes[p]['r2'].append(f)
            
    runs = []
    for p,vals in prefixes.items():
        if len(vals['r1'])!=1:
            raise Exception("%s has number of R1 files != 1 %s. Please check." % (p,vals['r1']))
        if len(vals['r2'])!=1:
            raise Exception("%s has number of R2 files != 1 %s. Please check." % (p,vals['r2']))
        runs.append(
            Sample(p,vals['r1'][0],vals['r2'][0]))
    return runs

def sort_out_single_files(files, r1_suffix=".fastq.gz"):
    prefixes = defaultdict(lambda:{"r1":[]})

    for f in files:
        tmp1 = re.search("(.+)%s" % r1_suffix,f)
        p = None
        if tmp1:
            p = tmp1.group(1).split("/")[-1]
            prefixes[p]['r1'].append(f)
            
    runs = []
    for p,vals in prefixes.items():
        if len(vals['r1'])!=1:
            raise Exception("%s has number of R1 files != 1 %s. Please check." % (p,vals['r1']))
        runs.append(
            Sample(p,vals['r1'][0],None))
    return runs

def find_fastq_files(directory,r1,r2):
    """
    Find fastq files in a directory and return a 
    list of tuples with the sample name and the 
    path to the fastq files from both pairs.
    """
    files = [f"{os.path.abspath(directory)}/{f}" for f in os.listdir(directory)]
    if r2:
        fastq_files = sort_out_paried_files(files,r1,r2)
    else:
        fastq_files = sort_out_single_files(files,r1)

    return fastq_files


def get_fasta_missing_content(fasta_file):
    """
    Get the missing content of a fasta file.
    """
    sys.stderr.write("Getting missing content of %s\n" % fasta_file)
    fasta = pp.Fasta(fasta_file)
    seq = list(fasta.fa_dict.values())[0]
    missing_characters = set(["n","N","-"])
    num_missing = 0
    for c in seq:
        if c in missing_characters:
            num_missing += 1
    missing_percent = round(num_missing / len(seq) * 100)
    return missing_percent

def mask_fasta(input,output,positions,newchrom=None):
    """
    Mask the fasta file with Ns.
    """
    sys.stderr.write("Masking %s\n" % input)
    fasta = pp.Fasta(input)
    seq = list(list(fasta.fa_dict.values())[0])

    if not newchrom:
        newchrom = list(fasta.fa_dict)[0]

    print(len(list(fasta.fa_dict.values())[0]))
    print(positions)
    for chrom,pos in positions:
        # print(pos)
        seq[pos-1] = "N"
    
    with open(output,"w") as O:
        O.write(">%s\n%s\n" % (newchrom,''.join(seq)))

def get_missing_positions(bam,depth_cutoff=50):
    tmpfile = "%s.bed" % uuid4()
    run_cmd(f"bedtools genomecov -ibam {bam} -d > {tmpfile}")
    positions = []
    for line in open(tmpfile):
        chrom,pos,depth = line.strip().split("\t")
        if int(depth) < depth_cutoff:
            positions.append((chrom,int(pos)))
    os.remove(tmpfile)
    return positions

def fasta_depth_mask(input,output,bam_file,depth_cutoff=50,newchrom=None):
    """
    Mask the fasta file with - based on the depth of the bam file.
    """
    logging.debug("Masking %s" % input)
    positions = get_missing_positions(bam_file,depth_cutoff=depth_cutoff)
    print(input)
    mask_fasta(input,output,positions,newchrom=newchrom)


def return_seqs_by_size(fasta_file,min_seq_size_cutoff,max_seq_size_cutoff):
    """
    Return a list of sequences from a fasta file that are 
    above a certain size cutoff.
    """
    fasta = pp.Fasta(fasta_file)
    seqs = {}
    for name,seq in fasta.fa_dict.items():
        if len(seq) > min_seq_size_cutoff and len(seq) < max_seq_size_cutoff:
            seqs[name] = seq
    return seqs

def get_megahit_contig_depth(contigs):
    """
    Get the depth of a contig from a megahit assembly.
    """
    depth = {}
    for line in open(contigs):
        if line.startswith(">"):
            r = re.search(">(\S+).+multi=(\S+)",line)
            contig = r.group(1)
            d = float(r.group(2))
            depth[contig] = d
    return depth

def filter_seqs_by_size(fasta_file,output, seqname, min_seq_size_cutoff, max_seq_size_cutoff):
    """
    Filter a fasta file by a size cutoff.
    """
    seqs = return_seqs_by_size(fasta_file,min_seq_size_cutoff=min_seq_size_cutoff,max_seq_size_cutoff=max_seq_size_cutoff)
    if len(seqs) == 1:    
        with open(output,"w") as O:
            for name,seq in seqs.items():
                O.write(">%s\n%s\n" % (seqname,seq))
        return True
    elif len(seqs) > 1:
        depths = get_megahit_contig_depth(fasta_file)
        highest_depth_contig = sorted([x for x in depths.items() if x[0] in seqs],key=lambda x: x[1],reverse=True)[0][0]
        with open(output,"w") as O:
            O.write(">%s\n%s\n" % (seqname,seqs[highest_depth_contig]))
        return True
    else:
        return False



class TaxonTree:
    """
    A data structure to store a taxonomic tree.
    
    Attributes
    ----------
    tree : dict
        A dictionary representing the tree. Keys are parent nodes and values are lists of child nodes.
    """
    def __init__(self,node_file='/Users/jody/.taxonkit/nodes.dmp'):
        """
        Parameters
        ----------
        node_file : str
            Path to the NCBI nodes.dmp file
        """
        self.tree = {}
        for l in open(node_file):
            row = [f.strip() for f in l.strip().split('|')]
            if row[0]==row[1]:
                continue
            child = int(row[0])
            parent = int(row[1])
            if parent not in self.tree:
                self.tree[parent] = []
            self.tree[parent].append(child)
    def find_descendants(self,node):
        """
        Recursively find all descendants of a given node.
        
        Parameters
        ----------
        node : int
            The node to find descendants of.
        
        Returns
        -------
        descendants : set
            A set of all descendants of the given node.
        """
        descendants = []
        children = self.tree.get(node, [])
        for child in children:
            descendants.append(child)
            # Recursively find descendants of the child
            descendants.extend(self.find_descendants(child))
        return set(descendants)
    
def filter_fastq_by_taxon(kraken_output: str,virus: int,reads: str,output: str):
    """
    Filter a fastq file by taxon.
    """
    nodedump = expanduser("~")+"/.dengue-ngs/taxdump/nodes.dmp"
    tree = TaxonTree(nodedump)

    include = [2732406]
    exclude = set(taxid2name.keys()) - set([virus])

    logging.debug(f"Filtering reads for virus {virus} with include nodes {include} and exclude nodes {exclude}")
    taxon_to_keep = set()
    for node in include:
        taxon_to_keep.add(int(node))
        taxon_to_keep.update(tree.find_descendants(int(node)))

    if exclude:
        for node in exclude:
            if int(node) in taxon_to_keep:
                taxon_to_keep.remove(int(node))
                taxon_to_keep.difference_update(tree.find_descendants(int(node)))

    read_names = set()
    for l in open(kraken_output):
        row = l.strip().split('\t')
        if int(row[2]) in taxon_to_keep:
            read_names.add(row[1])

    logging.debug(f"Found {len(read_names)} reads to keep for virus {virus} in {reads}")

    with pysam.FastxFile(reads) as fin, open(output, mode='w') as fout:
        for entry in fin:
            if entry.name in read_names:
                fout.write(str(entry) + '\n')