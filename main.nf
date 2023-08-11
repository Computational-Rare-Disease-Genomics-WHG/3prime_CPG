#!/usr/bin/env nextflow

params.config_file = "config.yaml"
config = readYaml(file: params.config_file)

// Channel for chromosome numbers
chroms = Channel.from(1..22, 'X', 'Y')

process extractColumns {
    input:
    val chrom from chroms

    output:
    file "${chrom}.tsv.gz" into extractedColumns

    script:
    """
    header=\$(zcat ${config.CADD}/whole_genome_SNVs_inclAnno.tsv.gz | head -n 2 | tail -n 1)
    desired_columns='#Chrom Pos Ref Alt Type Consequence FeatureID ConsScore GC CpG mamPhCons mamPhyloP GerpRS GerpRSpval GerpN GerpS'
    column_numbers=\$(echo "\$header" | awk -v cols="\$desired_columns" 'BEGIN {split(cols, arr); OFS="\t"} {for (i=1; i<=NF; i++) { for (j in arr) { if (\$i == arr[j]) printf "%s%s", i, (j == length(arr) ? "\\n" : OFS) } } }')
    output_file="${chrom}.tsv.gz"
    tabix -h ${config.CADD}/whole_genome_SNVs_inclAnno.tsv.gz \$chrom | cut -f"\${column_numbers}" | tail -n +2 | grep -F -f ${config.mane} | bgzip > "\${output_file}"
    """
}

process filterPossibleVariants {
    input:
    file input_file from extractedColumns.flatten()

    output:
    file "*.csv" into possibleVariants

    script:
    """
    filename=\$(basename "\$input_file")
    output_file="${config.output_dir}/possible_variants/\${filename%.tsv.gz}.csv"
    Rscript filter_possible_variants.R -i "\$input_file" -o "\$output_file"
    """
}

process gnomadProcessing {
    input:
    val chrom from chroms

    output:
    file "filtered_output.tsv" into filteredGnomad

    script:
    """
    vcf_file="${config.gnomAD}/gnomad.genomes.v3.1.1.sites.chr\${chrom}.vcf.bgz"
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/AC\\t%INFO/AF\\t%INFO/AN\\t%INFO/vep\\t%CHROM-%POS-%REF-%ALT\\n' "\$vcf_file" --include '((REF == "C" || REF == "G") && (ALT == "T" || ALT == "A"))' | grep -F -f ${config.mane} >> filtered_output.tsv
    """
}

process splitGnomad {
    input:
    file "filtered_output.tsv" from filteredGnomad

    script:
    """
    for chr in {1..22} X Y; do
        grep -P "^chr\${chr}\\t" filtered_output.tsv | cut -f1-8,10- > gnomad_chr\${chr}.txt
    done
    rm filtered_output.tsv
    """
}

process annotateVariants {
    input:
    val chrom from chroms
    file "*.csv" from possibleVariants.flatten()

    script:
    """
    Rscript annotate_variants.R -i data/possible_variants/\${chrom}.csv -g data/gnomad_chr\${chrom}.txt -m ${config.methylation} -o data/possible_variants_chr\${chrom}.txt -c \${chrom}
    """
}

process combineTSVFiles {
    input:
    file '*.txt' from annotateVariants.flatten().collect()

    output:
    file 'all_possible_cpg_variants.tsv.gz'

    script:
    """
    # Take the header from the first file
    cat \$(ls *.tsv.gz | head -n 1) | head -n 1 > all_possible_cpg_variants.tsv

    # Concatenate the content without headers
    for file in *.txt; do
        cat \$file | tail -n +2 >> all_possible_cpg_variants.tsv
    done

    # Compress the combined file
    bgzip all_possible_cpg_variants.tsv
    """
}


workflow {
    // Step 1: Extract the desired columns from CADD
    extractColumns()

    // Step 2: Filter the possible variants
    filterPossibleVariants()

    // Step 3: Process gnomAD data
    gnomadProcessing()

    // Step 4: Split the gnomAD data by chromosome
    splitGnomad()

    // Step 5: Annotate the variants with methylation level
    annotateVariants()

    // Step 6: Stack the chroms on top of each other
    combineTSVFiles()
}
