# Prep dataset annotations
# Usage example:
# Rscript src/prep_dataset_annotations.R \
#   -i data/possible_variants_all_chrs.txt \
#   -o data/../data/possible_variants_all_chrs_prepped.txt

# Libraries
library(data.table)
library(magrittr)
library(optparse)
library(scales)
library(tools)

# Set the number of threads to use
setDTthreads(0)

# Set up the command line options
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help =
                "Input file (all pos variants w/ gnomad and UKBB frequencies)"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


 # Read and process data
dt <- fread(opt$input)

# Fix methylation for those that are NA
dt[is.na(methyl_level), methyl_level := 0]

# Criteria for being observed_status
# If the variant is AC > 1 in either gnomad or UKBB
# If the variant is a singleton in both gnomad and UKBB

dt[, frequency_status := ifelse(
    (gnomad_ac > 1 | ukbb_ac > 1) |
    (gnomad_singleton == 1 & ukbb_singleton == 1),
    "Observed",
    "Not Observed")]

# Fix the consequence column
dt[,
    consequence := toTitleCase(
        gsub("PRIME", "'", gsub("_", " ", consequence)))
]

# Write the output
fwrite(dt, opt$output, sep = "\t")
