#Random Forest with 5-fold cross-validation
#This the RF is a modified version of EnvGen's code

## Set working directory
setwd("C:/Users/johan/OneDrive/R/ML/RF")

## Load libraries
library(lubridate)
library(ape)
library(phangorn)
library(seqinr)
library(vegan)
library(pheatmap)
library(randomForest)
library(ranger)
library(dplyr)
library(ggplot2)
library(caret)

## Set files
seqtab_file_18S = "C:/Users/johan/OneDrive/R/ML/RF/18S/18S/filtered_seqtab_18S.tsv" ## count matrix for each ASV & sample
taxa_file_18S = "C:/Users/johan/OneDrive/R/ML/RF/18S/18S/filtered_taxa_18S.tsv" ## taxonomic annotation for each ASV
seqtab_file_16S = "C:/Users/johan/OneDrive/R/ML/RF/16S/16S/filtered_seqtab_16S.tsv" ## count matrix for each ASV & sample
taxa_file_16S =  "C:/Users/johan/OneDrive/R/ML/RF/16S/16S/filtered_taxa_16S.tsv" ## taxonomic annotation for each ASV
phys_chem_file = "C:/Users/johan/OneDrive/R/ML/RF/physical_chemical_processed.tsv"
#bact_plan_file = "env_data/combined/bacterioplankton_processed.tsv"
#pico_plan_file = "env_data/combined/picoplankton_processed.tsv"
#phyt_plan_file = "env_data/combined/phytoplankton_processed.tsv"
#zoop_plan_file = "C:/Users/johan/OneDrive/R/ML/RF/zooplankton_processed.tsv"
zoop_plan_file = "C:/Users/johan/OneDrive/R/Master project/count_table_zooplankton_genus.tsv"
id_translation_file = "C:/Users/johan/OneDrive/R/ML/RF/physical_chemical_processed_translation.tsv"

asv_counts_18S = as.matrix(read.delim(seqtab_file_18S, row.names = 1))
asv_taxa_18S = as.matrix(read.delim(taxa_file_18S, row.names = 1))[,1:9]
asv_counts_16S = as.matrix(read.delim(seqtab_file_16S, row.names = 1))
asv_taxa_16S = as.matrix(read.delim(taxa_file_16S, row.names = 1))[,1:7]

colnames(asv_counts_18S) = gsub("^X", "", colnames(asv_counts_18S))
colnames(asv_counts_16S) = gsub("^X", "", colnames(asv_counts_16S))
rownames(asv_counts_18S) = paste(rownames(asv_counts_18S), "18S", sep = "_")
rownames(asv_taxa_18S) = paste(rownames(asv_taxa_18S), "18S", sep = "_")
rownames(asv_counts_16S) = paste(rownames(asv_counts_16S), "16S", sep = "_")
rownames(asv_taxa_16S) = paste(rownames(asv_taxa_16S), "16S", sep = "_")

norm_asv_counts_16S = t(t(asv_counts_16S)/colSums(asv_counts_16S))
norm_asv_counts_18S = t(t(asv_counts_18S)/colSums(asv_counts_18S))

identical(rownames(asv_counts_18S), rownames(asv_taxa_18S))
identical(rownames(asv_counts_16S), rownames(asv_taxa_16S))
identical(colnames(asv_counts_16S), colnames(asv_counts_18S))
colnames(asv_counts_16S) = gsub("^X", "", colnames(asv_counts_16S))
samples = colnames(asv_counts_16S)

alt_ids = as.matrix(read.delim(id_translation_file))
alt_ids = alt_ids[,c(1,ncol(alt_ids))]
samples_alt = alt_ids[match(samples, alt_ids[,1]),2]

phys_chem = read.delim(phys_chem_file)
#bact_plan = read.delim(bact_plan_file)
#pico_plan = read.delim(pico_plan_file)

#phyt_plan = read.delim(phyt_plan_file)
#colnames(phyt_plan) = gsub("^X", "", colnames(phyt_plan))
#colnames(phyt_plan) = gsub("\\.", "-", colnames(phyt_plan))
#phyt_plan = as.matrix(phyt_plan)

zoop_plan = read.delim(zoop_plan_file)
colnames(zoop_plan) = gsub("^X", "", colnames(zoop_plan))
colnames(zoop_plan) = gsub("\\.", "-", colnames(zoop_plan))
zoop_plan = as.matrix(zoop_plan)

## RANDOM FOREST WITH 5-FOLD CROSS VALIDATION 

run_randomforest_cv = function(features_matrix, responses_matrix, folds = 5) {
  n_samples = ncol(features_matrix)
  fold_size = ceiling(n_samples / folds)
  
  predicted_responses_matrix = matrix(nrow = nrow(responses_matrix), ncol = ncol(responses_matrix))
  
  for (fold in 1:folds) {
    start_idx = (fold - 1) * fold_size + 1
    end_idx = min(fold * fold_size, n_samples)
    
    features_subset = features_matrix[, start_idx:end_idx]
    responses_subset = responses_matrix[, start_idx:end_idx]
    
    for (i in 1:nrow(responses_subset)) {
      df = as.data.frame(cbind(responses_subset[i, ], t(features_subset)))
      rf = ranger(V1 ~ ., data = df, num.trees = 1000, importance = 'impurity')
      predicted_responses_matrix[i, start_idx:end_idx] = rf$predictions
    }
  }
  
  return(predicted_responses_matrix)
}


## Run randomforest predictions

## use one of the below as features_matrix_full
#features_matrix_full = asv_counts_16S 
#features_matrix_full = norm_asv_counts_16S
#features_matrix_full = asv_counts_18S
#features_matrix_full = norm_asv_counts_18S
features_matrix_full = rbind(norm_asv_counts_16S, norm_asv_counts_18S)

## use one of the below as responseses_matrix_full
responses_matrix_full = zoop_plan
#responses_matrix_full = phyt_plan

shared_samples = intersect(samples_alt, colnames(responses_matrix_full))
ix = match(shared_samples, samples_alt)
features_matrix = features_matrix_full[,ix]
binary_features_matrix = features_matrix
binary_features_matrix[which(features_matrix > 0)] = 1
ix = which(rowSums(binary_features_matrix) > 0.1*ncol(features_matrix))
features_matrix = features_matrix[ix,]
ix = match(shared_samples, colnames(responses_matrix_full))
responses_matrix = responses_matrix_full[,ix]
responses_matrix = responses_matrix[which(rowSums(responses_matrix) > 0),]

predicted_responses_matrix = run_randomforest_cv(features_matrix, responses_matrix, folds = 5)



#Correlation plot responses/predictions
cor_vec = c()
for (i in 1:nrow(responses_matrix)) {
  cor_vec[i] = cor.test(responses_matrix[i,], predicted_responses_matrix[i,])$est
}
names(cor_vec) = rownames(responses_matrix)

summary(cor_vec)

plot(log(rowSums(responses_matrix)), cor_vec)

png("C:/Users/johan/OneDrive/R/ML/RF/correlation_responses_species_240523.png", width = 10, height = 6, units = "in", res=200)
plot(log10(rowSums(responses_matrix)), cor_vec, 
     xlab = "log(responses)", 
     ylab = "Pearson correlation coefficient",
     main = "Correlation for responses")

dev.off()


cor_vec=as.matrix(cor_vec)
cor_vec=data.frame(cor_vec)
cor_vec$Taxa=rownames(cor_vec)
rownames(cor_vec) = 1:nrow(cor_vec)
cor_vec=cor_vec[,c("Taxa", "cor_vec")]
#cor_vec=cor_vec[cor_vec$cor_vec>0.6,]

#Make latex table

library(xtable)

cor_vec=cor_vec[order(cor_vec$Taxa), ]
# Generate LaTeX table
latex_table=xtable(cor_vec)

for (i in 1:nrow(cor_vec)) {
  latex_table[i, 1] = paste("\\textit{", as.character(cor_vec[i, 1]), "}", sep = "")
}

# Define a custom function to avoid escaping backslashes
sanitize_text = function(x) {
  gsub("\\", "", x, fixed = TRUE)
}

# Print the LaTeX code with the custom function
print(latex_table, include.rownames = FALSE, sanitize.text.function = sanitize_text)

#Table with actual and predicted values

actual_predicted_df = data.frame()

# Iterate over each cell in the matrices
for (i in 1:nrow(responses_matrix)) {
  for (j in 1:ncol(responses_matrix)) {
    # Get actual and predicted counts for the current cell
    actual_count = responses_matrix[i, j]
    predicted_count = predicted_responses_matrix[i, j]
    
    # Create a data frame with current sample, genus, actual count, and predicted count
    df = data.frame(Sample = colnames(responses_matrix)[j],
                     Genus = rownames(responses_matrix)[i],
                     Actual = actual_count,
                     Predicted = predicted_count)
    
    # Bind the current data frame to the result data frame
    actual_predicted_df = rbind(actual_predicted_df, df)
  }
}

write.table(actual_predicted_df, "C:/Users/johan/OneDrive/R/ML/RF/sample_predicted_240523.tsv", sep="\t")


#Plots for genera in metabar

genera_to_keep = c("Oithona", "Calanus", "Pseudocalanus", "Acartia", "Eurytemora", "Paracalanus", "Temora", "Centropages", "Penilia", "Oikopleura")

overall_correlation = actual_predicted_df %>%
  group_by(Genus) %>%
  summarize(
    Correlation = cor(Actual, Predicted, use = "pairwise.complete.obs"),
    P_Value = cor.test(Actual, Predicted)$p.value
  )

#Loop through each genus
for (genus in genera_to_keep) {
  #Subset current genus
  genus_data = actual_predicted_df[actual_predicted_df$Genus == genus, ]
  
  #Plot
  p = ggplot(genus_data, aes(x = Actual, y = Predicted)) +
    geom_point(color = "black", size = 2, shape = 16) +
    ggtitle(genus) +
    xlab("Actual microscopy count") +
    ylab("Predicted microscopy count") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 50, face = "bold")) +
    geom_smooth(method = "lm", se = TRUE, aes(fill = "Trend line")) +
    scale_fill_manual(values = "gray", guide = FALSE) +
    #Add correlation value and p-value to plot
    annotate("text", x = Inf, y = Inf, 
             label = paste0("R = ", round(overall_correlation$Correlation[overall_correlation$Genus == genus], 2), 
                            ", p = ", format.pval(overall_correlation$P_Value[overall_correlation$Genus == genus], digits = 2)), 
             color = "blue", hjust = 1, vjust = 1, size = 10)
  
  # Save each plot to a separate file
  ggsave(filename = paste0("C:/Users/johan/OneDrive/R/ML/RF/Correlation plots/", genus, "_RF_plot.png"), plot = p, bg='white', width = 10, height = 10)
}
