# Annotates the possible variants with
# UK Biobank frequency and methylation level

# Note that UK Biobank doesn't have any frequency
# information available for the Y chromosome

# Example usage:
# Rscript annotate_biobank_variants.R \
#   -i data/possible_variants_chr1_with_gnomad.txt \
#   -u data/ukbb_chr1.afreq \
#   -o data/possible_variants_chr1_with_gnomad_and_ukbb.txt \
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
    make_option(c("-u", "--ukbb"), type = "character", default = NULL,
                help = "UKBB variant file"),
    make_option(c("-c", "--chromosome"), type = "character", default = NULL,
                help = "Chromosome"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) ||
    is.null(opt$ukbb) ||
    is.null(opt$chromosome) ||
    is.null(opt$output)) {
    stop("Please provide input and output file paths using -i/--input, -g/--ukbb, and -o/--output options.") # nolint
}

if (opt$chromosome == "Y") {
    print("UK Biobank doesn't have any frequency information available for the Y chromosome. Adding columns but doing nothing") # nolint
}

# From the summary log within PLINK2
ukbb_total_samples <- 200031
ukbb_males <- 89918
ukbb_females <- 110088
ukbb_nosex <- 25
column_names <- c("chrom", "id", "ref", "alt", "ukbb_af", "ukbb_an")

# Create an empty data.table for Y
empty_dt <- as.data.table(
    setNames(
        data.frame(matrix(
            ncol = length(column_names), nrow = 0)),
        column_names
    )
)

# Read input data
selected_chrom <- opt$chromosome
dt <- fread(opt$input)

# Read in frequency table from
# UKBB, unless if its a "Y"
# then use the empty DT
# created above
ukbb <- data.table()

if (selected_chrom == "Y"){
    ukbb <- empty_dt
} else {
    ukbb <- fread(opt$ukbb)
}

# UK Biobank manipulation
names(ukbb) <- column_names

# Filter the data to possible CpG variants only
ukbb %<>% .[(ref == "C" & alt == "T") | (ref == "G" & alt == "C")]

# Calculate the Allele Count
if (selected_chrom == "X") {
    ukbb[, ukbb_ac := round(
        ukbb_af * ((ukbb_females + ukbb_nosex) * 2 + ukbb_males))]
} else {
    ukbb[, ukbb_ac := round(ukbb_af * ukbb_total_samples * 2)]
}

# Construct variant ID

# Get the position
grep_string <- paste0("chr", selected_chrom, ":|:SG")
ukbb[, pos := as.numeric(gsub(grep_string, "", id))]
ukbb[, id := NULL]

# Create the variant ID
ukbb[, variant_id := paste0(
    "chr", selected_chrom, "-", pos, "-", ref, "-", alt)]

# Select the columns
ukbb %<>% .[, .(variant_id, ukbb_ac, ukbb_af, ukbb_an)]

# Set the key
setkey(ukbb, variant_id)

# Merge the data
dt <- ukbb[dt]

# Create a column for the observed in UKBB
dt[is.na(ukbb_af), observed_in_ukbb := FALSE]
dt[!is.na(ukbb_af), observed_in_ukbb := TRUE]

# Write the data to file
fwrite(
    dt,
    opt$output,
    sep = "\t")
