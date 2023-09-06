# Annotates the possible variants with
# gnomAD allele frequency and methylation level

# Example usage:
# Rscript annotate_variants.R \
#   -i data/possible_variants/1.csv \
#   -g data/gnomad_chr1.txt \
#   -o data/possible_variants_chr1_with_gnomad.txt \
#   -c 1

library(data.table)
library(magrittr)
library(ggplot2)
library(optparse)

# Set up the number of threads
# to maximum available
setDTthreads(0)

# Set up the command line options
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input file"),
    make_option(c("-g", "--gnomad"), type = "character", default = NULL,
                help = "gnomAD file"),
    make_option(c("-c", "--chromosome"), type = "character", default = NULL,
                help = "Chromosome"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) ||
    is.null(opt$gnomad) ||
    is.null(opt$chromosome) || # nolint
    is.null(opt$output)) {
    stop("Please provide input and output file paths using -i/--input, -g/--gnomad, and -o/--output options.") # nolint
}

# Read input data
dt <- fread(opt$input)
gnomad <- fread(opt$gnomad)
selected_chrom <- opt$chromosome

# gnomAD manipulation
names(gnomad) <- c("chrom", "pos", "ref", "alt", "qual",
"gnomad_ac", "gnomad_af", "gnomad_an", "variant_id")
gnomad %<>% .[, .(variant_id, gnomad_ac, gnomad_af, gnomad_an)]
setkey(gnomad, variant_id)
setkey(dt, variant_id)

# Merging dt with gnomad
dt <- gnomad[dt]

# Observed in gnomAD
dt[is.na(gnomad_af), observed_in_gnomad := FALSE]
dt[!is.na(gnomad_af), observed_in_gnomad := TRUE]

# Fwrite to file
fwrite(
    dt,
    opt$output,
    sep = "\t")
