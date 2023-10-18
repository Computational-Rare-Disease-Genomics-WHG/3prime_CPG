library(data.table)
library(magrittr)
library(ggplot2)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input file"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$input) ||
    is.null(opt$output)) {
    stop("Please provide input and output file paths using -i/--input, -m/--methylation, and -o/--output options.") # nolint
}

# Reading in input file
dt <- fread(opt$input)

# Find multi allelic
multi_allelic_dt <- dt[alt %like% ","]
exploded_dt <- multi_allelic_dt[, .(
    chrom, pos, ref, qual, filter, info_an, info_pass_an,
    alt = unlist(strsplit(alt, ",")),
    info_ac = unlist(strsplit(info_ac, ",")),
    info_af = unlist(strsplit(info_af, ",")),
    info_pass_ac = unlist(strsplit(info_pass_ac, ","))
), by = seq_len(nrow(multi_allelic_dt))]
exploded_dt[, seq_len := NULL]

# Add back the split multi allelic
dt <- rbind(exploded_dt, dt)
dt %<>% .[!alt %like% ","]

# Remove the original (now redundant) rows with multiple alleles
dt %<>% .[!grepl(",", alt)]

# Convert columns back to their original types if needed
dt[, `:=`(
  info_ac = as.integer(info_ac),
  info_af = as.numeric(info_af),
  info_pass_ac = as.numeric(info_pass_ac)
)]

dt %<>% .[(ref == "C" & alt == "T") | (ref == "G" & alt == "A")]
fwrite(dt, opt$output, col.names = FALSE, sep = "\t")
