#########################################################################
# Title: Remove_known_spacers  
# Author: Sarah Cameron, University of Bath
# Email: sc3445@bath.ac.uk
# Last updated: 04/11/2024 
#########################################################################
# Description: 
# Removes any known expected flanking sequences of a repeat sequence e.g.
# leader, cloned spacer etc. These are listed in all_tags.csv and any
# sequence closely matching these (max mismatch 5) is removed also.
# This creates a .txt file with reads that contain unique flanking 
# regions, a Unique_Tags_New.fasta file with the unique sequences
# listed and a Total_Tags_New.fasta file with all flanking sequences 
# (includes duplicates of unique sequences)
#########################################################################

#Load packages
library(Biostrings)
library(tibble)

#Argument 1 is a database of all the known flanking regions of repeats 
args = commandArgs(trailingOnly = TRUE)
args[1]='all_tags.csv'
Known_Tags <- read.csv(args[1], header=FALSE, stringsAsFactors = FALSE)
Known_Tags <- DNAStringSet(Known_Tags$V1)
Known_Rev <- DNAStringSet(reverseComplement(Known_Tags))
Known_Tags <- c(Known_Tags, Known_Rev)

#Argument 2 is a database of all the flanking regions calculated from the raw reads 
args[2]='condition_4_matches.BLAST.Result.txtAll_Tags_New.csv'
Raw_Tags <- read.csv(args[2], stringsAsFactors = FALSE)


#We then need to get all the unique sequences from the raw reads and make it into a DNAStringSet 
Ups <- Raw_Tags$Upstream_Tag
Downs <- Raw_Tags$Downstream_Tag
Up_uni <- unique(Ups)
Down_uni <- unique(Downs)
up_down <- c(Up_uni, Down_uni)
total_uni <- unique(up_down)
DNA_URAW <- DNAStringSet(total_uni)


#We can use vwhichPDict to use the dictionary of patterns to search the DNAStringSet for which rows have that pattern or a inexact match of it 
what <- vwhichPDict(Known_Tags, DNA_URAW, max.mismatch = 5, min.mismatch = 0,with.indels = TRUE,fixed=FALSE, algorithm = "auto")

#Using the tibble package and function enframe we can create a dataframe of the list generated 
#install.packages("tibble")

WHAT <- enframe(what)

#Then using which we can find out the rows that got no matches i.e. are seemingly unique to the raw reads
rows_unique <- which(WHAT$value == "integer(0)")

#Then create a dataframe of these unique tags 
Unique_Tags <- as.character(DNA_URAW[rows_unique,])
Unique_Tags <- data.frame(Unique_Tags)
fasta <- DNAStringSet(Unique_Tags$Unique_Tags)

#To find which original reads have novel tags in 
read_rows <- which(Raw_Tags$Upstream_Tag %in% Unique_Tags$Unique_Tags | Raw_Tags$Downstream_Tag %in% Unique_Tags$Unique_Tags)
reads <- Raw_Tags[read_rows,]
reads_acq <- unique(reads$Read_ID)

#Write the dataframe of tags only to csv 
writeXStringSet(fasta, paste0(args[2], "Unique_Tags_New.fasta"),format = "fasta")
#Write a list of the read ids with acquisitions 
write.table(reads_acq, paste0(args[2], "read_ids.txt"), row.names=F, quote = F)


#Repeat the above but collect all the flanking sequences without removing duplicates 

raw_reads <- c(Ups, Downs)
DNA_URAW <- DNAStringSet(raw_reads)

#We can use vwhichPDict to use the dictionary of patterns to search the DNAStringSet for which rows have that pattern or a inexact match of it 
what <- vwhichPDict(Known_Tags, DNA_URAW, max.mismatch = 5, min.mismatch = 0,with.indels = TRUE,fixed=FALSE, algorithm = "auto")

#Using the tibble package and function enframe we can create a dataframe of the list generated 
#install.packages("tibble")

WHAT <- enframe(what)

#Then using which we can find out the rows that got no matches i.e. are seemingly unique to the raw reads
rows_unique <- which(WHAT$value == "integer(0)")

#Then create a dataframe of these unique tags 
Unique_Tags <- as.character(DNA_URAW[rows_unique,])
Unique_Tags <- data.frame(Unique_Tags)
fasta <- DNAStringSet(Unique_Tags$Unique_Tags)


#Write the dataframe of tags only to csv 
writeXStringSet(fasta, paste0(args[2], "Total_Tags_New.fasta"),format = "fasta")