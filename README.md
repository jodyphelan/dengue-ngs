# dengue-ngs

## Install 

```
mamba create -n dengue samtools pilon bwa lofreq ncbi-datasets-cli sourmash tqdm
mamba activate dengue
pip3 install git+https://github.com/jodyphelan/dengue-ngs.git
```

## Set up reference database
```
dengue-download-ref.py
```

## Run
```
dengue-ngs.py -f folder_with_fastq_files
```

