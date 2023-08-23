# filter_possible_variants.R
# Filters the possible variants from CADD
# to only CpG variants

library(data.table)
library(optparse)

# Set threads to all available
setDTthreads(0)

# Argument parsing
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input file path"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$input) || is.null(opt$output)) {
  stop("Please provide input and output file paths using -i/--input and -o/--output options.")
}

# Read input data
dt <- fread(opt$input)

cols <- c("chrom", "pos", "ref", "alt", "type",
  "consequence", "CpG", "cons_score",
  "gc", "transcript_id", "mam_phastcons", "mam_phylop",
  "gerp_rs", "gerp_rs_pval", "gerp_n", "gerp_s",
  "spliceai-acc-gain", "spliceai-acc-loss",
  "spliceai-don-gain", "spliceai-don-loss",
  "cadd_raw", "cadd_phred")

names(dt) <- cols

# Rscript to find CpG variants

# Forward strand CpG
# Filter rows where REF is a 'C' and the next base is a 'G'
forward_strand_cpg <- dt[ref == 'C' & shift(ref, type = 'lead', fill = '') == 'G']

# Filter rows where REF is a 'G' and the previous base is a 'C'
reverse_strand_cpg <- dt[ref == 'G' & shift(ref, type = 'lag', fill = '') == 'C']

# Filter to CpG variants
forward_strand_cpg[, mC_strand := '+']
reverse_strand_cpg[, mC_strand := '-']
cpg_dt <- rbind(forward_strand_cpg, reverse_strand_cpg)

# Filter to variants
cpg_dt <- cpg_dt[
  (mC_strand == '+' & ref == 'C' & alt == 'T') |
  (mC_strand == '-' & ref == 'G' & alt == 'A')
]

# Add variant id
cpg_dt[, variant_id := paste0(chrom, "-", pos, "-", ref, "-", alt)]

# Write output data
fwrite(cpg_dt, opt$output, sep="\t", row.names = F)
