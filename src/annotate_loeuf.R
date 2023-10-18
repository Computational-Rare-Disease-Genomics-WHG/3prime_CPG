# Merges contraint data with the variant data

library(data.table)
library(magrittr)

# Set the number of threads to use
setDTthreads(0)

# Set up the command line options
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help =
                "Input file (all pos variants w/ gnomad and UKBB frequencies)"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file"),
    make_option(c("-l", "--loeuf"), type = "character", default = NULL,
                help = "LOEUF file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read in the dataset
dt <- fread(opt$input, sep = "\t")

# Find the transcript set
transcript_list <- unique(dt$transcript_id)

# Get LOEUF scores
constraint <- fread(opt$loeuf)
constraint %<>% .[,
    .(transcript, oe_lof_lower)]
constraint <- unique(constraint)
names(constraint) <- c("transcript_id", "loeuf")
constraint %<>% .[transcript_id %in% transcript_list]

# Determine the decile breaks
deciles <- quantile(constraint$loeuf, prob = seq(
    0, 1, length.out = 10), na.rm = TRUE)

# Create the loeuf decile column
constraint[, loeuf_decile := as.numeric(
    cut(
        loeuf,
        deciles,
        include.lowest = TRUE,
        label = FALSE
        )
    )
]

# Merge the constraint data
setkey(dt, transcript_id)
setkey(constraint, transcript_id)
dt <- constraint[dt]


# Write the output
fwrite(dt, opt$output, sep = "\t")
