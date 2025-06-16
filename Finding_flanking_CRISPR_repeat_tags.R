#########################################################################
# Title: Finding_flanking_CRISPR_repeat_tags  
# Author: Sarah Cameron, University of Bath
# Email: sc3445@bath.ac.uk
# Last updated: 04/11/2024 
#########################################################################
# Description: 
# Imports a BLAST result table in outfmt 6 with sstd information which
# has used the known 28bp CRISPR repeat region as query. The script then
# finds the flanking 32bp either side of the BLAST match and saves them 
# in the BLAST table. 
#########################################################################
library(dplyr)
library(Biostrings)

# Set working directory with BLAST results
setwd(getwd()) 

#First import the BLAST data file to the environment which has the repeat matches and the fasta file of reads
args = commandArgs(trailingOnly = TRUE) 
args[1]="condition_10_matches.BLAST.Result.txt"
args[2]="condition_10_trimmed_rename.fasta"

BLAST <- read.delim(args[1],stringsAsFactors=F, head=F)
#Assign the column names 
colnames(BLAST) <- c("Query_ID", "Read_ID", "%_Identical", "Align_Length", "No_Mismatches", "No_Gaps", "Start_Align_Query", "End_Align_Query", "Start_Align_Read", "End_Align_Read", "EValue", "BitScore", "Strand", "Read_Length")
reads_with_repeats <- length(unique(BLAST$Read_ID))

#Now we want to separate the BLAST results into Fwd and Rev strands 
Rev <- which(BLAST$Strand == 'minus')
Rev_Strand <- BLAST[Rev,]
Fwd_Strand <- BLAST[-Rev,]

#To speed up the tag extraction process we want to remove reads where no acquisitions are likely to have taken place
# i.e. 
repeats_in_read <- Fwd_Strand %>% count(Fwd_Strand$Read_ID)
acq <- which(repeats_in_read$n >= 3)
singles <- repeats_in_read[-acq,]
no_acq <- which(Fwd_Strand$Read_ID %in% singles$`Fwd_Strand$Read_ID`)
Fwd_Strand <- Fwd_Strand[-no_acq,]

#Filter reads out if <534 bp 
too_short <- which(Fwd_Strand$Read_Length <= 500) 
if (length(too_short) >= 1)
  Fwd_Strand <- Fwd_Strand[-too_short,]


#Add some columns to work out the upstream tag 
Fwd_Strand$Start_Up_Tag <- NA
Fwd_Strand$End_Up_Tag <- NA 
Fwd_Strand$Upstream_Tag <- NA 

#Calculate the end of the upstream tag 
Fwd_Strand$End_Up_Tag <- Fwd_Strand$Start_Align_Read - 1 

#Calculate the start of the upstream tag
Fwd_Strand$Start_Up_Tag <- Fwd_Strand$End_Up_Tag - 31

#Now we need to add a control step where the script will not try and work out any tags for alignments that have ended at the end of the read (aka have minus numbers in the position columns)
End_of_Reads <- which(Fwd_Strand$Start_Up_Tag < 1 | Fwd_Strand$End_Up_Tag < 32)
if (length(End_of_Reads) >= 1)
  Fwd_Strand <- Fwd_Strand[-End_of_Reads,]

Genome_Seq <- readDNAStringSet(args[2], format="fasta")

#Loop to extract the upstream tags 
for (line in c(1:nrow(Fwd_Strand)))
{
  Read_No <- Fwd_Strand$Read_ID[line]
  Fwd_Strand$Upstream_Tag[line] <- as.character(Genome_Seq[[Read_No]][Fwd_Strand$Start_Up_Tag[line]:Fwd_Strand$End_Up_Tag[line]])
}

#Find the downstream 32bps 
Fwd_Strand$Start_Down_Tag <- NA
Fwd_Strand$End_Down_Tag <- NA 
Fwd_Strand$Downstream_Tag <- NA 
#Calculate the start of the downstream tag 
Fwd_Strand$Start_Down_Tag <- Fwd_Strand$End_Align_Read + 1 

#Calculate the end of the downstream tag
Fwd_Strand$End_Down_Tag <- Fwd_Strand$Start_Down_Tag + 31

too_short<-which(Fwd_Strand$Read_Length < Fwd_Strand$End_Down_Tag)
if (length(too_short) >= 1)
  Fwd_Strand <- Fwd_Strand[-too_short,]

#Loop to extract the downstream tags 
for (line in c(1:nrow(Fwd_Strand)))
{
  Read_No <- Fwd_Strand$Read_ID[line]
  Fwd_Strand$Downstream_Tag[line] <- as.character(Genome_Seq[[Read_No]][Fwd_Strand$Start_Down_Tag[line]:Fwd_Strand$End_Down_Tag[line]])
}


#Do the same for the reverse strand 
repeats_in_read <- Rev_Strand %>% count(Rev_Strand$Read_ID)
acq <- which(repeats_in_read$n >= 3)
singles <- repeats_in_read[-acq,]
no_acq <- which(Rev_Strand$Read_ID %in% singles$`Rev_Strand$Read_ID`)
Rev_Strand <- Rev_Strand[-no_acq,]

#Filter reads out if <534 bp 
too_short <- which(Rev_Strand$Read_Length <= 500) 
if (length(too_short) >= 1)
  Rev_Strand <- Rev_Strand[-too_short,]

#Add some columns to work out the upstream tag 
Rev_Strand$Start_Up_Tag <- NA
Rev_Strand$End_Up_Tag <- NA 
Rev_Strand$Upstream_Tag <- NA 

#Calculate the end of the upstream tag 
Rev_Strand$End_Up_Tag <- Rev_Strand$End_Align_Read - 1 

#Calculate the start of the upstream tag
Rev_Strand$Start_Up_Tag <- Rev_Strand$End_Up_Tag - 31 

too_short <- which(Rev_Strand$Start_Up_Tag < 0)
if (length(too_short) >= 1)
  Rev_Strand <- Rev_Strand[-too_short,]


#Loop to extract the upstream tags 
for (line in c(1:nrow(Rev_Strand)))
{
  Read_No <- Rev_Strand$Read_ID[line]
  Up_Tag <- DNAStringSet(Genome_Seq[[Read_No]][Rev_Strand$End_Up_Tag[line]:Rev_Strand$Start_Up_Tag[line]])
  Rev_Strand$Upstream_Tag[line] <- reverseComplement(Up_Tag)
}

#Find the downstream 32bps 
Rev_Strand$Start_Down_Tag <- NA
Rev_Strand$End_Down_Tag <- NA 
Rev_Strand$Downstream_Tag <- NA 
#Calculate the start of the downstream tag 
Rev_Strand$Start_Down_Tag <- Rev_Strand$Start_Align_Read + 1 

#Calculate the end of the downstream tag
Rev_Strand$End_Down_Tag <- Rev_Strand$Start_Down_Tag + 31

too_short<-which(Rev_Strand$Read_Length < Rev_Strand$End_Down_Tag)
if (length(too_short) >= 1)
  Rev_Strand <- Rev_Strand[-too_short,]

#Loop to extract the downstream tags 
for (line in c(1:nrow(Rev_Strand)))
{
  Read_No <- Rev_Strand$Read_ID[line]
  Down_Tag <- DNAStringSet(Genome_Seq[[Read_No]][Rev_Strand$Start_Down_Tag[line]:Rev_Strand$End_Down_Tag[line]])
  Rev_Strand$Downstream_Tag[line] <- reverseComplement(Down_Tag)
}


#Merge both upstream tables together 
Both <- merge(Fwd_Strand, Rev_Strand, all = TRUE)

#Write to csv 
write.table(Both, paste(args[1],"All_Tags_New.csv", sep=""),row.names = F,quote = F,sep = ",")

