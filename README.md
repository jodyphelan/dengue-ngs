# dengue-ngs

## Install 

```
mamba create -n dengue samtools pilon bwa lofreq ncbi-datasets-cli sourmash tqdm bedtools bcftools parallel blast megahit
mamba activate dengue
pip3 install --force-reinstall  git+https://github.com/jodyphelan/dengue-ngs.git
```

## Set up reference database
```
dengue-download-ref.py
```

## Run
```
dengue-ngs.py -f folder_with_fastq_files
```

