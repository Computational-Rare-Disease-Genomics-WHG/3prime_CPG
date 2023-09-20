# 3 Prime CPG Variants

## Background

TBD


## Pipeline

Firstly, you can install the input data using the download_init.sh script available.

```bash

# Download reference sequence
./download_init.sh --data refseq --output "./data" --yaml_file "config.yaml"

# Download refseq
./download_init.sh --data refseq_assembly --output "./data" --yaml_file "config.yaml"

# Download mane
./download_init.sh --data mane --output "./data" --yaml_file "config.yaml"

# Download loeuf
./download_init.sh --data loeuf --output "./data" --yaml_file "config.yaml"

# Download cpg_island
./download_init.sh --data cpg_island --output "./data" --yaml_file "config.yaml"

# Download dnase
./download_init.sh --data dnase --output "./data" --yaml_file "config.yaml"

# Download polya
./download_init.sh --data polya --output "./data" --yaml_file "config.yaml"

# Download ancestral sequence
./download_init.sh --data ensembl_ancestral --output "./data" --yaml_file "config.yaml"

# Download miRNA
./download_init.sh --data mirna --output "./data" --yaml_file "config.yaml"

# Download Histone Methylation
./download_init.sh --data histone_methylation --output "./data" --yaml_file "config.yaml"


```

This installs all of the data dependencies that are required. There is however additional methylation and UK Biobank / gnomAD data that is required to the pipeline.

Update `config.yaml` for the path to the respective datasets once downloaded.

The pipeline within this repository can be reproduced using Nextflow. Install Nextflow from their website. Dependencies include the version of Java recommended by Nextflow.

There is a reproducible docker image available at `docker pull elstonndsouza/3primeCpG:latest`. Although, this Docker image can be also be rebuilt if installation is difficult on your target platform of choice.

```bash

nextflow run main.nf -c config.yaml -profile local

```

## Plan

TBD

1. Filter MANE to 3' UTR regions coordinates.
2. Download gnomAD and filter using tabix to UTR regions.
3. Annotate variants with methylation level using the testis methylation dataset frome ENCODE project.
4. Annotate with CADD / PhyloP and other annotation tools.
5. Machine learning on finding those that aren't observed.
6. Profit??
