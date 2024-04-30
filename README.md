# dengue-ngs

## Install 

```
wget https://raw.githubusercontent.com/jodyphelan/dengue-ngs/main/conda/linux-env.txt
conda create --name dengue --file linux-env.txt
conda activate dengue
pip install --force-reinstall  git+https://github.com/jodyphelan/pathogen-profiler.git
pip install --force-reinstall  git+https://github.com/jodyphelan/dengue-ngs.git
```

## Set up reference database
```
dengue-download-ref.py
```

## Run
```
dengue-ngs.py -f folder_with_fastq_files
```

