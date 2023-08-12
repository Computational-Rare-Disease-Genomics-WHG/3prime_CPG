#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("rtracklayer")

#BiocManager::install("Biostrings")

#install.packages("seqinr")
#install.packages("tidyverse")
#install.packages("stringr")
#install.packages("tidyr")
#install.packages("splitstackshape")
#install.packages("beepr")

library(beepr)
# beep(3)

library(dplyr)
library(tidyr)

library(splitstackshape)

library(tidyverse)
library(seqinr)
library(Biostrings)
library(rtracklayer)
library(stringr)

#reverse compliment function
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"
    if(xx[bbb]=="C") y[bbb]<-"G"
    if(xx[bbb]=="G") y[bbb]<-"C"
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE)
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)
}

setwd("C:/Users/makbo/Desktop/Whiffin")

### generating the 3 UTR regions file for tabix (it's been modified so check the output dataframe is correct format)
#read in GFF file
data <- readGFF("data/MANE.GRCh38.v1.0.ensembl_genomic.gff")

#choose 3' UTRs
three_UTRs <- data[data@listData$type == 'three_prime_UTR',]

#create dataframe with chr, start, end for each 3' UTR (for tabix)
UTR_subset <- data.frame(three_UTRs@listData$seqid, three_UTRs@listData$start,
                            three_UTRs@listData$end, three_UTRs@listData$strand)

colnames(UTR_subset) <- c("chr","UTR_start","UTR_end","strand")

UTR_subset$strand <- sub("$", "1", UTR_subset$strand)

UTR_subset$chr <- gsub("chr", "", UTR_subset$chr)


#tab delim file with 3' UTR positions (chr, start, end) (change UTR_subset code to get table in correct format)
# write.table(UTR_subset, file="data/UTR_coords.bed", sep="\t", row.names = FALSE,
#             col.names = FALSE, quote = FALSE)

###reading in and cleaning up MANE transcript data and sequences
transcripts <- readDNAStringSet(file = "data/MANE.GRCh38.v1.0.ensembl_rna.fna")

s <- strsplit(names(transcripts), "[ \t]+")

#iterate through each "row" of the "list of lists"
for (r in 1:length(s)){
  #empty string for description
  desc <- ""
  #iterate through words from 'description:' to the end of the word list
  for (w in s[[r]][8:(length(s[[r]]))]){
    #add each word after description: to the same string
    desc <- paste(desc, w)
  }
  #replace the description: with the new string
  s[[r]][8] <- desc
  #remove everything after the description string
  s[[r]] <- s[[r]][-(9:(length(s[[r]])))]
  #split the description string into the actual description and source
  desc_source <- strsplit(s[[r]][8], "\\[")
  #clean up description
  s[[r]][8] <- str_trim(desc_source[[1]][1], side = "both")
  #clean up source
  s[[r]][9] <- substr(desc_source[[1]][2],1,(nchar(desc_source[[1]][2])-1))
}
beep(3)

#dataframe with all the information of the transcripts
info <- data.frame(matrix(unlist(s), ncol=9, byrow=T))
rownames(info) <- c(1:length(names(transcripts)))
colnames(info) <- c("ID","cDNA","position","gene","gene_biotype",
                    "transcript_biotype","gene_symbol","description","HGNC")

#only look at interesting columns
info <- info[,c(1,3,4,7,8,9)]

#remove unecessary words
info$gene <- gsub("gene:", "", info$gene)
info$gene_symbol <- gsub("gene_symbol:", "", info$gene_symbol)
info$description <- gsub("description:", "", info$description)
info$HGNC <- gsub("Source:HGNC Symbol;Acc:HGNC:", "", info$HGNC)

#separate out the position information
info$position <- gsub("chromosome:GRCh38:", "", info$position)
positions <- as.data.frame(str_split_fixed(info$position, ":", 4))
colnames(positions) <- c("chr","start","end","strand")

#combine transcript and position info
info <- cbind(info, positions)[,c(-2)]
#add a + symbol in front of forward strand transcripts
info$strand[info$strand != "-1"] <-
  sub("^", "+", info$strand[info$strand != "-1"])

#make positions into integers
class(info$start) <- "numeric"
class(info$end) <- "numeric"

#add transcript sequences to the table and remove duplicates
info$sequence <- transcripts

info <- subset(info, !duplicated(info$ID))

#select important columns
info <- info[,c(1,6,9,10)]

#Nechama's way of processing the data
mane <- read.delim("data/MANE.GRCh38.v1.0.ensembl_genomic.gff", header = FALSE)


v9 <- as.data.frame(str_split_fixed(mane$V9, ";", 13))
colnames(v9) <- c("ID", "parent", "gene_id", "transcript_id", "gene_type",
                  "gene_name", "transcript_type", "transcript_name", "exon_number",
                  "exon_id", "tag", "protein_id")

mane <- cbind(mane[-9],v9)

ms <- mane[(mane$V3 == "three_prime_UTR"),]

ms <- ms[!grepl("Plus_Clinical", ms$tag),]

ms$gene_name <- gsub("gene_name=", "", ms$gene_name)
ms$transcript_id <- gsub("transcript_id=", "", ms$transcript_id)

#grouping exons, calculating exon length, and getting positions of exons

ms$exon_number <- 1
ms$exon_length <- NA
ms$total_exon_length <- NA
ms$exon_positions <- NA

pb <- txtProgressBar(min=0,max = nrow(ms), initial = 1)

for (r in 1:nrow(ms)){
  if (ms$exon_number[r] != 1){
    next
  }
  x <- TRUE
  count <- 1
  ms$exon_length[r] <- abs(ms$V5[r] - ms$V4[r])+1
  exon_lens <- c(abs(ms$V5[r] - ms$V4[r]) +1)
  if (ms$V7[r] == "+") {
    exon_positions <- c(ms$V4[r]:ms$V5[r])
  } else if (ms$V7[r] == "-") {
    exon_positions <- c(ms$V5[r]:ms$V4[r])
  }
  while (x == TRUE && r < nrow(ms)){
    if (identical(ms$transcript_id[r], ms$transcript_id[r+count])){
      ms$exon_number[r+count] <- 1+count
      ms$exon_length[r+count] <- abs(ms$V5[r+count] - ms$V4[r+count])+1
      exon_lens <- c(exon_lens, abs(ms$V5[r+count] - ms$V4[r+count]) +1)
      if (ms$V7[r] == "+") {
        exon_positions <- c(exon_positions, ms$V4[r+count]:ms$V5[r+count])
      } else if (ms$V7[r] == "-") {
        exon_positions <- c(exon_positions, ms$V5[r+count]:ms$V4[r+count])
      }
      count <- count+1
      x<-TRUE
    } else {
      ms$total_exon_length[r:(r+count)] <- sum(exon_lens)
      ms$exon_positions[r:(r+count)] <- list(exon_positions)
      x <- FALSE
    }
  }
  setTxtProgressBar(pb,r)
}
beep(3)

# df with transcript ID, exon positions, and exon lengths
UTR_length <- ms[c(12,23,24)][!duplicated(ms),]

#put the transcript data with the 3' UTR data into one df
info_ms <- merge(x = info, y = UTR_length,
                 by.x = "ID", by.y = "transcript_id")

info_ms$UTR_sequence <- NA

# get the 3' UTR sequences
pb <- txtProgressBar(min=0,max = nrow(info_ms), initial = 1)

for (r in 1:nrow(info_ms)){
  if (info_ms$strand[r] == "+1"){
    info_ms$UTR_sequence[r] <- str_sub(info_ms$sequence[r], -info_ms$total_exon_length[r], -1)
  } else if (info_ms$strand[r] == "-1") {
    info_ms$UTR_sequence[r] <- rev.comp(str_sub(info_ms$sequence[r], -info_ms$total_exon_length[r], -1), rev = FALSE)
  }
  setTxtProgressBar(pb,r)
}
beep(3)

# df with coordinates and 3' UTR sequences

info_new <- info_ms[,c(2, 6, 7)]

info_new$UTR_sequence <- strsplit(as.character(info_new$UTR_sequence), "")
# df with per base coordinates, in 3' UTRs
info_new <- unnest(info_new, c(UTR_sequence, exon_positions))

info_new <- info_new[!duplicated(info_new),]

info_new <- info_new[order(info_new$chr, info_new$exon_positions),]

colnames(info_new) <- c("chromosome","position","ref_base")

# write.table(info_new, file="data/base_positions_UTR.txt", sep="\t", row.names = FALSE,
#             col.names = TRUE, quote = FALSE)

#############################
info_new <- read.delim("data/base_positions_UTR.txt", header = TRUE)

# methylation data in 3' UTRs
UTR_meth_data <- as.data.frame(read.table(file="data/testis_cov3_38_3_UTRs.bed.gz"))

methylation <- UTR_meth_data[,c(1,3,5,6)]
colnames(methylation) <- c("chromosome", "position", "methylation", "base_strand")
methylation$chromosome <- gsub("chr", "", methylation$chromosome)
methylation <- methylation[!duplicated(methylation),]

# match methylation data to the base by position
combined_data <- merge(info_new, methylation, by = c("chromosome", "position"), all=T)

# remove extra entries (~470) that occured during merge (not sure why this happened)
combined_data <- inner_join(combined_data, info_new)

# write.table(combined_data, file="data/base_positions_methylation_UTR.txt", sep="\t", row.names = FALSE,
#             col.names = TRUE, quote = FALSE)

###########################
combined_data <- read.delim("data/base_positions_methylation_UTR.txt", header = TRUE)
