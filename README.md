# 3 Prime CPG Variants

## Background

CpG methylation plays a crucial role in gene regulation and genomic stability. It can influence gene expression by affecting the binding of transcription factors and other regulatory proteins to DNA. Generally, high levels of methylation in a specific region of DNA are associated with gene silencing, while low levels or absence of methylation allow gene expression. CpG sites can undergo methylation, leading to stable gene silencing. CpG sites are also hypermutable, given that they undergo deamination variants (C->T) at a much higher rate than other sites. 3’ CpG variants are not well understood in the literature. This project examines their methylation levels and compares methylated sites with a C-T variant to those without, assuming evolutionary pressures have constrained these sites.

## Research Questions

1. Can a machine learning model predict whether a site falls into one of two classes: Variant or Invariant?
   
2. What sequence features are necessary to train this model?
   - Local Sequence Context
   - Presence of a CpG Island
   - PolyA site
   - Conservation scores
   - [Additional factors to be researched]

3. How do certain learning models perform in comparison to other methods?

4. What are the clinical implications of the findings?

5. How to address the dataset's unbalanced size, common in CpG sites datasets?

## Steps

1. Define 3’ UTR regions using MANE (Matched Annotation from NCBI and EMBL-EBI).

2. Use methylation data from available sources to annotate site.

3. Define "methylated sites" based on a specified methylation level (e.g., 60% or higher). Extract genomic locations and relevant sequence features.

4. Create a pipeline to calculate the described features.

5. Utilize gnomAD version 3.1.2 and / or UK BB  identify 3’ UTR methylated variants (from previously defined regions) that are observed and unobserved.

## Relevant Papers

- [Mutation saturation for fitness effects at human CpG sites](https://elifesciences.org/articles/71513)
- [Intergenic, gene terminal, and intragenic CpG islands in the human genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2817693/#B24)
- [A comprehensive catalog of CpG islands methylated in human lung adenocarcinomas for the identification of tumor suppressor genes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4042824/)
- [A curated census of pathogenic and likely pathogenic UTR variants and evaluation of deep learning models for variant effect prediction](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8460052/)
  

## Pipeline

Firstly, you can install the input data using the download_init.sh script available.

```bash

data_types=("refseq" "refseq_assembly" "mane" "loeuf" "cpg_island" "dnase" "polya" "ensembl_ancestral" "mirna" "histone_methylation")

for data_type in "${data_types[@]}"; do
    ./download_init.sh --data "$data_type" --output "./data" --yaml_file "config.yaml"
done

```

This installs all of the data dependencies that are required. There is however additional methylation and UK Biobank that is required to the pipeline.

Update `config.yaml` for the path to the respective datasets once downloaded.

The pipeline within this repository can be reproduced using Nextflow. Install Nextflow from their website. Dependencies include the version of Java recommended by Nextflow.

There is a reproducible docker image available at `docker pull elstonndsouza/3primeCpG:latest`. Although, this Docker image can be also be rebuilt if installation is difficult on your target platform of choice.

```bash

nextflow run main.nf -c config.yaml -profile local

```

## Machine Learning 

TODO
