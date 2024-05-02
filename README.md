Welcome! 

This is the GitHub page for my master project "Using metabarcoding data to elucidate spatiotemporal distribution
patterns of zooplankton in the Baltic Sea area", where I explore metabarcoding and microscopy data (years 2015-2020). 
The project was conducted at SciLifeLab in The Environmental Genomics group.

The code is organised according to following steps:

1. "read_merged_18S.R" - This file is pre-processing of the metabarcoding data and the code is a modified pipeline from
   The Environmental Genomics group. The output are processed files to use in step 2. 

2. "reads_comparisons_18S_metazoa.R" - Investigation of the distribution of annotated reads for different taxonomic levels for all 18S vs. metazoa in reference database PR^2.
   "reads_comparisons_18S_metazoa_SILVA.R" - Investigation of the distribution of annotated reads for different taxonomic levels for all 18S vs. metazoa in SILVA.
   (The codes are similar)

3. "count_table_zooplankton.R" - Summarized observations of genus/specie in each sample for microscopy data obtained from SHARKweb. 

4. "common_samples.R" - processing of data for common samples. NMDS for microscopy and metabarcoding genera/ASVs. 

5. "comparison_of_methods" - comparative barplots for microscopy and metabarcoding and scatter plots. 

