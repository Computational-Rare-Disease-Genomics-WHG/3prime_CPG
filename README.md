# 3 Prime CPG Variants 

## Background 

TBD


## Pipeline 

The following pipeline can be reproduced using Nextflow. Install Nextflow from their website. Dependencies include the version of Java recommended by Nextflow. 

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