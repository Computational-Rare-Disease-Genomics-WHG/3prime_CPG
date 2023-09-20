#!/bin/bash

set -e

# Function for displaying usage instructions
usage() {
    echo "Usage: $0 --data [data_type] --output [directory_to_output] --yaml_file [yaml_file to append path]"
    exit 1
}

# Function to download "cadd" data
download_cadd() {
    local output="$1"
    cadd_file="https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz"
    cadd_index="https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz.tbi"

    # Download CADD data and its index
    wget "$cadd_file" -O "${output}/whole_genome_SNVs_inclAnno.tsv.gz" || exit 1
    wget "$cadd_index" -O "${output}/whole_genome_SNVs_inclAnno.tsv.gz.tbi" || exit 1
}

# Function to download "refseq" data
download_refseq() {
    local output="$1"
    genomes_url="ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"

    # Download RefSeq data
    wget "$genomes_url" -O "${output}/human_genome_grch38.fna.gz" || exit 1
}

# Function to download "refseq" assembly report
download_refseq_assembly_report() {
    local output="$1"
    assembly_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt"

    # Download assembly report
    wget "$assembly_url" -O "$output/assembly_report.txt" || exit 1
}

# Function to download gnomAD vcf files
download_gnomad() {
    local output="$1"
    gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"

    # Download gnomAD data for chromosomes 1-22, X, and Y
    for i in {1..22} X Y; do
        wget "$gnomad_url/gnomad.genomes.v3.1.2.sites.chr${i}.vcf.bgz" -O "${output}/gnomad.genomes.v3.1.2.sites.chr${i}.vcf.bgz" || exit 1
    done
}

# Function to download MANE file
download_mane() {
    local output="$1"
    mane_url="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0"

    # Download MANE data
    wget "$mane_url/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz" -O "${output}/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz" || exit 1
    wget "$mane_url/MANE.GRCh38.v1.0.ensembl_rna.fna.gz" -O "${output}/MANE.GRCh38.v1.0.ensembl_rna.fna.gz" || exit 1
    wget "$mane_url/MANE.GRCh38.v1.0.summary.txt.gz" -O "${output}/MANE.GRCh38.v1.0.summary.txt.gz" || exit 1
}

# Function to download LOEUF scores
download_loeuf() {
    local output="$1"
    loeuf_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"

    # Download LOEUF scores
    wget "$loeuf_url" -O "${output}/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz" || exit 1
}

# Function to download CpG island data
download_cpg_island() {
    local output="$1"
    cpg_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz"

    # Download CpG island data
    wget "$cpg_url" -O "${output}/cpgIslandExtUnmasked.txt.gz" || exit 1
}

# Function to download DNase data
download_dnase() {
    local output="$1"
    dnase_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeRegDnaseClustered.txt.gz"

    # Download DNase data
    wget "$dnase_url" -O "${output}/wgEncodeRegDnaseClustered.txt.gz" || exit 1
}

download_polya(){
    local output="$1"

    poly_a_miner="https://raw.githubusercontent.com/YalamanchiliLab/PolyA-miner/master/Human_hg38.PolyA_DB.bed"
    poly_a_db="https://exon.apps.wistar.org/polya_db/v3/download/3.2/human_pas.zip"
    unibas_poly_asite="https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz"

    # Download PolyA data
    wget "$poly_a_miner" -O "${output}/Human_hg38.PolyA_DB.bed"
    wget "$poly_a_db" -O "${output}/human_pas.zip" 
    wget "$unibas_poly_asite" -O "${output}/atlas.clusters.2.0.GRCh38.96.bed.gz"

    #unzip the human_pas.zip
    unzip "${output}/human_pas.zip" -d "${output}/" 

    # remove the zip file
    rm "${output}/human_pas.zip"
}

download_histone_methylation(){
    local output="$1"
    histone_url="https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/E062-H3K9me3.broadPeak.gz"
    wget "$histone_url" -O "${output}/E062-H3K9me3.broadPeak.gz"
    gunzip "${output}/E062-H3K9me3.broadPeak.gz"
    mv -- "${output}/E062-H3K9me3.broadPeak" "${output}/E062-H3K9me3.bed"
}

download_ancestral_sequence(){
    local output="$1"
    ancestral_url="https://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz"
    wget "$ancestral_url" -O "${output}/homo_sapiens_ancestor_GRCh38.tar.gz"
    tar -xvf "${output}/homo_sapiens_ancestor_GRCh38.tar.gz" -C "${output}"
}


download_mirna(){
    local output="$1"
    mirna_url="https://www.targetscan.org/vert_80/vert_80_data_download/All_Target_Locations.hg19.bed.zip"
    wget "$mirna_url" -O "${output}/All_Target_Locations.hg19.bed.zip"
    unzip "${output}/All_Target_Locations.hg19.bed.zip" -d "${output}/All_Target_Locations"
}

# Function to append data to the YAML file
append_to_yaml() {
    local data_type="$1"
    local filepath="$2"
    local yaml_file="$3"
    echo "$data_type: $filepath" >> "$yaml_file"
}

# Main script
DATA=""
OUTPUT=""
YAML_FILE=""

# Parsing arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --data) DATA="$2"; shift ;;
        --output) OUTPUT="$2"; shift ;;
        --yaml_file) YAML_FILE="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Ensure all arguments are provided
if [[ -z "$DATA" || -z "$OUTPUT" || -z "$YAML_FILE" ]]; then
    usage
fi

# Determine which data type to download and call the appropriate function
case $DATA in
    cadd)
        download_cadd "$OUTPUT"
        append_to_yaml "cadd" "$(realpath "${OUTPUT}/whole_genome_SNVs_inclAnno.tsv.gz")" "$YAML_FILE"
        append_to_yaml "cadd_index" "$(realpath "${OUTPUT}/whole_genome_SNVs_inclAnno.tsv.gz.tbi")" "$YAML_FILE"
        ;;
    refseq)
        download_refseq "$OUTPUT"
        append_to_yaml "refseq" "$(realpath "${OUTPUT}/human_genome_grch38.fna.gz")" "$YAML_FILE"
        ;;
    refseq_assembly)
        download_refseq_assembly_report "$OUTPUT"
        append_to_yaml "refseq_assembly" "$(realpath "$OUTPUT/assembly_report.txt")" "$YAML_FILE"
        ;;

    ensembl_ancestral)
        download_ancestral_sequence "$OUTPUT"
        append_to_yaml "ensembl_ancestral" "$(realpath "$OUTPUT/homo_sapiens_ancestor_GRCh38")" "$YAML_FILE"
        ;;

    mane)
        download_mane "$OUTPUT"
        append_to_yaml "mane_summary" "$(realpath "$OUTPUT/MANE.GRCh38.v1.0.summary.txt.gz")" "$YAML_FILE"
        append_to_yaml "mane_transcript" "$(realpath "$OUTPUT/MANE.GRCh38.v1.0.ensembl_rna.fna.gz")" "$YAML_FILE"
        append_to_yaml "mane_gff" "$(realpath "$OUTPUT/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz")" "$YAML_FILE"
        ;;
    gnomad)
        download_gnomad "$OUTPUT"
        append_to_yaml "gnomad_dir" "$(realpath "$OUTPUT/")" "$YAML_FILE"
        ;;
    loeuf)
        download_loeuf "$OUTPUT"
        append_to_yaml "loeuf" "$(realpath "$OUTPUT/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz")" "$YAML_FILE"
        ;;
    cpg_island)
        download_cpg_island "$OUTPUT"
        append_to_yaml "cpg_island" "$(realpath "$OUTPUT/cpgIslandExtUnmasked.txt.gz")" "$YAML_FILE"
        ;;
    dnase)
        download_dnase "$OUTPUT"
        append_to_yaml "dnase_peaks" "$(realpath "$OUTPUT/wgEncodeRegDnaseClustered.txt.gz")" "$YAML_FILE"
        ;;
    polya)
        download_polya "$OUTPUT"
        append_to_yaml "polya_miner" "$(realpath "$OUTPUT/Human_hg38.PolyA_DB.bed")" "$YAML_FILE"
        append_to_yaml "polya_db" "$(realpath "$OUTPUT/human.PAS.txt")" "$YAML_FILE"
        append_to_yaml "polya_site_unibas" "$(realpath "$OUTPUT/atlas.clusters.2.0.GRCh38.96.bed.gz")" "$YAML_FILE"
        ;;
    mirna)
        download_mirna "$OUTPUT"
        append_to_yaml "mirna" "$(realpath "$OUTPUT/All_Target_Locations")" "$YAML_FILE"
        ;;
    histone_methylation)
        download_histone_methylation "$OUTPUT"
        append_to_yaml "histone_methylation" "$(realpath "$OUTPUT/E062-H3K9me3.bed")" "$YAML_FILE"
        ;;
    *)
        echo "Unknown data type: $DATA"; usage ;;
esac

echo "Download complete!"
