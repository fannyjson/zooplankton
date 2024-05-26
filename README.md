**Welcome! **

This is the GitHub page for my master project "Using metabarcoding data to elucidate spatiotemporal distribution
patterns of zooplankton in the Baltic Sea area", where I explore metabarcoding and microscopy data from the Swedish marine monitoring program. 

1. "pre_processing.R" : This file is pre-processing of the metabarcoding data and the code is a modified pipeline from
   The Environmental Genomics group. The output files are used in the following steps.

2. "annotated_reads_comparison.R" : Annotated reads for non-metazoa and metazoa for different taxonomic levels. Output used in step 4.

3. "count_table_zooplankton.R" : Summarized observations of taxa in each sample for microscopy data obtained from SHARKweb. Count table for 'genus' was further used in steo 4 and 5. 

4. "common_samples.R" : Processing of data for common samples. Output was used in step 7.

5. "comparison_databases.R" : Matching-tables for SILVA, PR^2 and microscopy. 

6. "RF_5fold.R" : Random Forest predictions using 5-fold cross validation. The output was used in step 7. 

7. "comparison_of_methods.R" : Comparative barplots of different kinds for microscopy and metabarcoding + scatter plots and Pearson correlations between methods. 

