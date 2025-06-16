# spacer_acquistions
Analysis of spacer acquisitons for Elliott et al, 2025: 

1. *Finding_flanking_CRISPR_repeat_tags.R* - Takes BLAST output of reads with repeat sequence and finds flanking 32bp either side of repeat
2. *Remove_known_spacers.R* - Removes any of the flanking regions (or closely match) those that are expected e.g. leader, cloned in spacer
3. *Acquisition_mapping.sh* - Contains the commands for mapping the flanking regions to the phage genome
4. *Read_depth_plots.R* - Code for creating read depth plots of spacers that map to corresponding phage genome 
