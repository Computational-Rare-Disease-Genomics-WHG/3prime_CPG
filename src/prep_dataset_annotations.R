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
library(logger)


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


log_info("Loading dataset")
 # Read and process data
dt <- fread(opt$input)
log_info("Finished loading dataset")


log_info("Data processing")


log_info("Fixing methylation")
# Fix methylation for those that are NA
dt[is.na(methyl_level), methyl_level := 0]



log_info("Fixing NAs in frequencies")
## Fix data_types of NA
dt[
    is.na(gnomad_ac) | is.na(gnomad_af) | gnomad_ac == 0 | gnomad_af == 0, `:=`(
        gnomad_ac = 0,
        gnomad_af = 0
    )
]

dt[
    is.na(ukbb_ac) | is.na(ukbb_ac) | ukbb_ac == 0 | ukbb_af == 0, `:=` (
        ukbb_ac = 0,
        ukbb_af = 0
    )
]

# Create a column for the observed in UKBB
dt[ukbb_ac == 0, observed_in_ukbb := FALSE]
dt[ukbb_ac != 0, observed_in_ukbb := TRUE]
dt[gnomad_ac == 0, observed_in_gnomad := FALSE]
dt[gnomad_af != 0, observed_in_gnomad := TRUE]


# Criteria for being observed_status
# If the variant is AC > 1 in either gnomad or UKBB
# If the variant is a singleton in both gnomad and UKBB

dt[, observation_status := ifelse(
    (gnomad_ac > 1 | ukbb_ac > 1) |
    (gnomad_ac == 1 & ukbb_ac == 1),
    "Observed",
    "Not Observed")]





# Fix the consequence column
dt[,
    consequence := toTitleCase(
        gsub("PRIME", "'", gsub("_", " ", consequence)))
]

log_info("Finished data processing. Writing to file")

# Write the output
fwrite(dt, opt$output, sep = "\t")

log_info("Finished writing to file")
