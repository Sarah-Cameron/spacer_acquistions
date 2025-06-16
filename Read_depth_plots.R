#########################################################################
# Title: Read_depth_plots   
# Author: Sarah Cameron, University of Bath
# Email: sc3445@bath.ac.uk
# Last updated: 04/11/2024 
#########################################################################
# Description: 
# Imports the read depth coverage file from samtools to create plot of
# coverage across the phage genome and saves as an .svg file.
#########################################################################

library(ggplot2)
library(dplyr)
library(svglite)

# Load the samtools depth output file
depth_data <- read.table("pool_1_all_read_depth.coverage", header = FALSE, col.names = c("chrom", "pos", "depth"))

# Create the read depth plot
ggplot(depth_data, aes(x = pos, y = depth)) +
  geom_bar(color='#512DA8', stat='identity') +
  scale_y_continuous(limits = c(0, 15),breaks = seq(0, 15, 5))+
  scale_x_continuous()+
  labs(x = "Phi2 Phage Genome Position", y = "Read Depth", title = "Pool 1 Total Spacers")+
  theme_minimal()

ggsave("Pool_1_Total.svg", width=10, height=4, bg='white', device="svg")


#Other colours for each genome 
#0B6625
#D21F3C