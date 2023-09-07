# Annotates the possible variants with
# gnomAD UK Biobank frequency and methylation level

# Note that UK Biobank doesn't have any frequency
# information available for the Y chromosome

# Example usage:
# Rscript annotate_frequencies.R \
#   -i data/possible_variants/possible_variants_chr${chr}_with_methylation.txt \
#   -u data/ukbb/ukbb_chr${chr}.afreq \
#   -g data/gnomad/gnomad_chr${chr}.txt \
#   -o data/possible_variants_chr${chr}_with_af.txt \
#   -c ${chr}

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
    is.null(opt$gnomad) ||
    is.null(opt$chromosome) ||
    is.null(opt$output)) {
    stop("Please provide input and output file paths using -i/--input, -g/--ukbb, and -o/--output options.") # nolint
}

if (opt$chromosome == "Y") {
    print("UK Biobank doesn't have any frequency information available for the Y chromosome. Adding columns but doing nothing") # nolint
}

# From the summary log within PLINK2
# Needed to calculate the allele count
ukbb_total_samples <- 200031
ukbb_males <- 89918
ukbb_females <- 110088
ukbb_nosex <- 25

ukbb_column_names <- c("chrom", "id", "ref", "alt", "ukbb_af", "ukbb_an")
gnomad_column_names <- c("chrom", "pos", "ref", "alt", "qual",
"gnomad_ac", "gnomad_af", "gnomad_an", "variant_id")

# Create an empty data.table for ukbb Y
empty_dt <- as.data.table(
    setNames(
        data.frame(matrix(
            ncol = length(ukbb_column_names), nrow = 0)),
        ukbb_column_names
    )
)

# Read input data
selected_chrom <- opt$chromosome
dt <- fread(opt$input)
gnomad <- fread(opt$gnomad)


# Work through gnomAD manipulation
names(gnomad) <- gnomad_column_names
gnomad %<>% .[, .(variant_id, gnomad_ac, gnomad_af, gnomad_an)]

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
names(ukbb) <- ukbb_column_names

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

# Create the variant ID
ukbb[, variant_id := paste0(
    "chr", selected_chrom, "-", pos, "-", ref, "-", alt)]

# Select the columns
ukbb %<>% .[, .(variant_id, ukbb_ac, ukbb_af, ukbb_an)]

# Merge the data
setkey(gnomad, variant_id)
setkey(dt, variant_id)
setkey(ukbb, variant_id)
dt <- ukbb[dt]
dt <- gnomad[dt]


# Create a column for the observed in UKBB
dt[is.na(ukbb_af), observed_in_ukbb := FALSE]
dt[!is.na(ukbb_af), observed_in_ukbb := TRUE]
dt[is.na(gnomad_af), observed_in_gnomad := FALSE]
dt[!is.na(gnomad_af), observed_in_gnomad := TRUE]

# Write the data to file
fwrite(
    dt,
    opt$output,
    sep = "\t")
