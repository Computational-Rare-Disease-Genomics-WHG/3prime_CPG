# Annotates the possible variants with
# gnomAD allele frequency and methylation level

# Example usage:
# Rscript annotate_variants.R \
#   -i data/possible_variants/1.csv \
#   -g data/gnomad_chr1.txt \
#   -m data/context_prepared_methyl14_weighted_logit.methyl_level.bed.gz \
#   -o data/possible_variants_chr1.txt \
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
    make_option(c("-m", "--methylation"), type = "character", default = NULL,
                help = "Methylation file"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$input) ||
    is.null(opt$gnomad) ||
    is.null(opt$methylation) ||
    is.null(opt$output)) {
    stop("Please provide input and output file paths using -i/--input, -g/--gnomad, -m/--methylation, and -o/--output options.") # nolint
}

# Read input data
dt <- fread(opt$input)
gnomad <- fread(opt$gnomad)
meth_dt <- fread(opt$methylation)
selected_chrom <- opt$chrom

# Setting up the possible variants
dt[, variant_id := paste0("chr", variant_id)]

# Some weird character conversion needs to happen
# to fix this for whatever reason
dt[, chrom:= as.character(chrom)]
dt[, chrom := paste0("chr", selected_chrom)]
setkey(dt, variant_id)

# gnomAD manipulation
names(gnomad) <- c("chrom", "pos", "ref", "alt", "qual",
"info_ac", "info_af", "info_an", "variant_id")
gnomad %<>% .[, .(variant_id, qual, info_ac, info_af, info_an)]
setkey(gnomad, variant_id)

# Setting up methylation data
names(meth_dt) <- c("chrom", "pos", "pos_2", "methyl_level")
meth_dt %<>% .[chrom == paste0("chr", selected_chrom)]
meth_dt <- rbind(meth_dt[, .(chrom, pos, methyl_level)],
    meth_dt[, .(chrom, pos_2, methyl_level)], use.names = FALSE)
meth_dt$chrom <- paste0("chr", selected_chrom)
meth_dt <- unique(meth_dt)
setkey(meth_dt, chrom, pos)


# Merging dt with gnomad
dt <- gnomad[dt]

# Meging dt with methylation
setkey(meth_dt, chrom, pos)
setkey(dt, chrom, pos)
dt <- meth_dt[dt]

# Observed in gnomAD
dt[is.na(info_af), observed_in_gnomad := FALSE]
dt[!is.na(info_af), observed_in_gnomad := TRUE]

# Fwrite to file
fwrite(
    dt,
    opt$output,
    sep = "\t")
