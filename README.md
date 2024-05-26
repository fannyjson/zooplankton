**Welcome! **

This is the GitHub page for my master project "Using metabarcoding data to elucidate spatiotemporal distribution
patterns of zooplankton in the Baltic Sea area", where I explore metabarcoding and microscopy data from the Swedish marine monitoring program. 

1. "pre_processing.R" - This file is pre-processing of the metabarcoding data and the code is a modified pipeline from
   The Environmental Genomics group. The output files are used in the following steps. 

2. "annotated_reads_comparison.R" - Investigation of the distribution of annotated reads for different taxonomic levels for all 18S vs. metazoa in reference database PR^2.

3. "count_table_zooplankton.R" - Summarized observations of genus/specie in each sample for microscopy data obtained from SHARKweb. 

4. "common_samples.R" : Processing of data for common samples. NMDS plots for ASVs/taxa.

5. "RF_5fold.R" : Random Forest predictions using 5-fold cross validation. The output from this was used in step 6. 

6. "comparison_of_methods.R" : Comparative barplots of different kinds for microscopy and metabarcoding + scatter plots and Pearson correlations between methods. 

7. "databases_comparison.R" : Matching-tables for SILVA, PR^2 and microscopy.
