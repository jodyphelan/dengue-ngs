FROM mambaorg/micromamba:latest
LABEL image.name="jodyphelan/dengue-ngs"

# This needs to be set for the PATH variable to be set correctly
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# create a directory for the internal data used by the container
USER root
RUN mkdir /input_data
RUN mkdir /output_data


# Install our software
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    pigz \
    git \
    zip \
    unzip \
    samtools \
    pilon \
    bwa \
    lofreq \
    ncbi-datasets-cli \
    sourmash \
    tqdm \
    bedtools \
    bcftools \
    parallel \
    blast \
    megahit \
    kmcp \
    kraken2 \
    seqkit \
    krakentools \
    pandas \
    plotly \
    python-kaleido \
    minimap2 \
    pyyaml \
    procps-ng \
    && \
    micromamba clean --all --yes

# Install pipeline
RUN pip install git+https://github.com/jodyphelan/dengue-ngs.git
RUN pip install git+https://github.com/jodyphelan/pathogen-profiler.git
RUN dengue-download-ref.py --threads 2 


# set `/data` as working directory so that the output is written to the
# mount point when run with `docker run -v $PWD:/data ... -o output.csv`
WORKDIR /output_data

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "dengue-ngs.py","--jobs","1", "--threads", "4", "--folder", "/input_data"]


