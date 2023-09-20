# Defines the 3' UTR regions of transcripts in the MANE GFF file

library(data.table)
library(magrittr)
library(optparse)

setDTthreads(0)


# Named command line arguments
option_list <- list(
  make_option(
    c("-m", "--mane-tsv"),
    type = "character",
    help = "Mane GFF file (as TSV)",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output-file"),
    dest = "output_file",
    type = "character",
    help = "Output file name (should end in .tsv)",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
