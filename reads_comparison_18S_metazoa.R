# Comparison of reads from 18S and only metazoan 18S 
# The files loaded are output files from read_merged_18S.R 

setwd("C:/Users/johan/OneDrive/R/Master project") 

#Taxa table
taxa_18S=read.table("C:/Users/johan/OneDrive/R/Master project/taxa_18S.tsv") 
taxa_18S = taxa_18S[!is.na(taxa_18S[, 1]), ]

#Normalized seqtab table metazoa 
norm_seqtab_18S=read.table("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_18S_240308.tsv")

#Normalized seqtab all 
norm_seqtab_18S_all=read.table("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_18S_all_240308.tsv")

#Original seqtab metazoa
seqtab_18S=read.table("C:/Users/johan/OneDrive/R/Master project/seqtab_18S_240308.tsv")

#Original seqtab all
seqtab_18S_all=read.table("C:/Users/johan/OneDrive/R/Master project/seqtab_18S_all_240308.tsv")

#Load libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

#################################################################################

#Fix sample names 

colnames(norm_seqtab_18S) = gsub(colnames(norm_seqtab_18S), pattern = '^X', replacement = '')
colnames(norm_seqtab_18S_all) = gsub(colnames(seqtab_18S_all), pattern = '^X', replacement = '')

colnames(seqtab_18S) = gsub(colnames(seqtab_18S), pattern = '^X', replacement = '')
colnames(seqtab_18S_all) = gsub(colnames(seqtab_18S_all), pattern = '^X', replacement = '')

replace_xx <- function(vector) {
  gsub("(_X|_XX)", "", vector)
}

taxa_18S[] <- lapply(taxa_18S, replace_xx)

replace <- function(vector) {
  gsub("_", " ", vector)
}


taxa_18S[] <- lapply(taxa_18S, replace)

#Number of reads for annotated 18S reads
all_reads=sum(seqtab_18S_all)
metazoa_reads=sum(seqtab_18S)

########################## #18S ALL ############################################

#Make plot of relative abundance of for different taxa levels 
#Match seqtab with taxa table
#Loop through each column of seqtab_processed and store indices>0. 

taxa_list=list()

for (i in 1:ncol(norm_seqtab_18S_all)) {
  present_indices=which(norm_seqtab_18S_all[,i] > 0)
  
  if (length(present_indices)>0) {
    present_asvs=rownames(norm_seqtab_18S_all)[present_indices]
    
    if(any(present_asvs %in% rownames(taxa_18S))){
      present_asvs=present_asvs[present_asvs %in% rownames(taxa_18S)]
      
      present_taxa=taxa_18S[rownames(taxa_18S) %in% present_asvs,]
      taxa_list[[i]]=data.frame(sample_id=colnames(norm_seqtab_18S)[i],ASV=present_asvs,
                                Taxa=present_taxa)
    }
  }
}


#Summarising counts at different taxonomic levels
#Normalized clade count

clade_counts_norm=list()

for (i in 1:ncol(taxa_18S)){
  matr_norm=NULL
  clade=unique(taxa_18S[,i])
  clade=clade[!is.na(clade)]
  for (j in 1:length(clade)){
    ix=which(clade[j]==taxa_18S[,i])
    if (length(ix)>1){
      matr_norm=rbind(matr_norm, apply(norm_seqtab_18S_all[ix,],2,sum,na.rm=TRUE))
    } else {
      matr_norm=rbind(matr_norm, norm_seqtab_18S_all[ix,])
    }
  }
  rownames(matr_norm)=clade
  
  clade_counts_norm[[i]]=matr_norm
}

matr_norm[is.na(matr_norm)]=0
matr_norm=matr_norm[rowSums(matr_norm !=0, na.rm=TRUE)>0 , , drop=FALSE]

#Non-normalized clade count

clade_counts=list()

for (i in 1:ncol(taxa_18S)){
  matr=NULL
  clade=unique(taxa_18S[,i])
  clade=clade[!is.na(clade)]
  for (j in 1:length(clade)){
    ix=which(clade[j]==taxa_18S[,i])
    if (length(ix)>1){
      matr=rbind(matr, apply(seqtab_18S_all[ix,],2,sum,na.rm=TRUE))
    } else {
      matr=rbind(matr, seqtab_18S_all[ix,])
    }
  }
  rownames(matr)=clade
  #colnames(matr)=colnames(matr)=sample_id
  
  clade_counts[[i]]=matr
}

matr[is.na(matr)]=0
matr=matr[rowSums(matr !=0, na.rm=TRUE)>0 , , drop=FALSE]


# Total reads for each taxa (18S)

domain_reads = na.omit(clade_counts[[1]])
domain_reads = sum(domain_reads)

supergroup_reads = na.omit(clade_counts[[2]])
supergroup_reads = sum(supergroup_reads)

division_reads = na.omit(clade_counts[[3]])
division_reads = sum(division_reads)

subdivision_reads = na.omit(clade_counts[[4]])
subdivision_reads = sum(division_reads)

class_reads = na.omit(clade_counts[[5]])
class_reads = sum(class_reads)

order_reads = na.omit(clade_counts[[6]])
order_reads = sum(order_reads)

family_reads = na.omit(clade_counts[[7]])
family_reads= sum(family_reads)

genus_reads = na.omit(clade_counts[[8]])
genus_reads = sum(genus_reads)

species_reads = na.omit(clade_counts[[9]])
species_reads = sum(species_reads)


# Create a list to store taxonomic levels and sums
taxonomic_levels <- list()
taxonomic_sums <- list()

# Loop through each clade_counts element
for (i in 1:9) {
  # Extract clade counts
  clade <- clade_counts[[i]]
  
  # Convert to data frame and remove NAs
  clade_df <- na.omit(data.frame(clade))
  
  # Calculate sum for each taxonomic level
  taxonomic_sums[[i]] <- sum(unlist(clade_df), na.rm = TRUE)
  
  # Store taxonomic levels
  taxonomic_levels[[i]] <- clade_df
}

# Calculate sum of all reads
sum_domain_reads <- sum(unlist(taxonomic_sums), na.rm = TRUE)

# Calculate fractions
fractions <- sapply(taxonomic_sums, function(category_sums) category_sums / taxonomic_sums[[1]] * 100)
# Create a data frame with row names
fractions_df <- data.frame(fractions)
rownames(fractions_df) <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

print(fractions_df)
                   
#################### 18S (METAZOA ONLY) #######################################

#18S metazoa only
taxa_metazoa = taxa_18S[taxa_18S[, 4] == "Metazoa", ]
taxa_metazoa = taxa_metazoa[!is.na(taxa_metazoa[, 1]), ]
write.csv(taxa_metazoa, "C:/Users/johan/OneDrive/R/Master project/taxa_metazoa.csv", row.names = FALSE)

meta_taxa_list=list()

for (i in 1:ncol(seqtab_18S)) {
  present_indices=which(norm_seqtab_18S[,i] > 0)
  
  if (length(present_indices)>0) {
    present_asvs=rownames(norm_seqtab_18S)[present_indices]
    
    if(any(present_asvs %in% rownames(taxa_metazoa))){
      present_asvs=present_asvs[present_asvs %in% rownames(taxa_metazoa)]
      
      present_taxa=taxa_metazoa[rownames(taxa_metazoa) %in% present_asvs,]
      taxa_list[[i]]=data.frame(sample_id=colnames(norm_seqtab_18S)[i],ASV=present_asvs,
                                Taxa=present_taxa)
    }
  }
}

#Summarising counts at different taxonomic levels
#Normalized

norm_meta_clade_counts=list()

for (i in 1:ncol(taxa_metazoa)){
  norm_meta_matr=NULL
  clade=unique(taxa_metazoa[,i])
  clade=clade[!is.na(clade)]
  for (j in 1:length(clade)){
    ix=which(clade[j]==taxa_metazoa[,i])
    if (length(ix)>1){
      norm_meta_matr=rbind(norm_meta_matr, apply(norm_seqtab_18S[ix,],2,sum,na.rm=TRUE))
    } else {
      norm_meta_matr=rbind(norm_meta_matr, norm_seqtab_18S[ix,])
    }
  }
  
  rownames(norm_meta_matr)=clade
  #colnames(matr)=colnames(norm_matr)=sample_id
  
  norm_meta_clade_counts[[i]]=norm_meta_matr
}

norm_meta_matr[is.na(norm_meta_matr)]=0
norm_meta_matr=norm_meta_matr[rowSums(norm_meta_matr !=0, na.rm=TRUE)>0 , , drop=FALSE]

#Non-normalized

meta_clade_counts=list()

for (i in 1:ncol(taxa_metazoa)){
  meta_matr=NULL
  clade=unique(taxa_metazoa[,i])
  clade=clade[!is.na(clade)]
  for (j in 1:length(clade)){
    ix=which(clade[j]==taxa_metazoa[,i])
    if (length(ix)>1){
      meta_matr=rbind(meta_matr, apply(seqtab_18S[ix,],2,sum,na.rm=TRUE))
    } else {
      meta_matr=rbind(meta_matr, seqtab_18S[ix,])
    }
  }
  
  rownames(meta_matr)=clade
  
  meta_clade_counts[[i]]=meta_matr
}


meta_matr[is.na(meta_matr)]=0
meta_matr=matr[rowSums(meta_matr !=0, na.rm=TRUE)>0 , , drop=FALSE]

# Total reads for each taxa (18S)

sum(seqtab_18S)

domain_reads_meta = na.omit(meta_clade_counts[[1]])
domain_reads_meta = sum(domain_reads_meta)

supergroup_reads_meta = na.omit(meta_clade_counts[[2]])
supergroup_reads_meta = sum(supergroup_reads_meta)

division_reads_meta = na.omit(meta_clade_counts[[3]])
division_reads_meta = sum(division_reads_meta)

subdivision_reads_meta = na.omit(meta_clade_counts[[4]])
subdivision_reads_meta = sum(division_reads_meta)

class_reads_meta = na.omit(meta_clade_counts[[5]])
class_reads_meta = sum(class_reads_meta)

order_reads_meta = na.omit(meta_clade_counts[[6]])
order_reads_meta = sum(order_reads_meta)

family_reads_meta = na.omit(meta_clade_counts[[7]])
family_reads_meta= sum(family_reads_meta)

genus_reads_meta = na.omit(meta_clade_counts[[8]])
genus_reads_meta = sum(genus_reads_meta)

species_reads_meta = na.omit(meta_clade_counts[[9]])
species_reads_meta = sum(species_reads_meta)


# Create a list to store taxonomic levels and sums
taxonomic_levels_meta <- list()
taxonomic_sums_meta <- list()

# Loop through each clade_counts element
for (i in 1:9) {
  # Extract clade counts
  meta_clade <- meta_clade_counts[[i]]
  
  # Convert to data frame and remove NAs
  meta_clade_df <- na.omit(data.frame(meta_clade))
  
  # Calculate sum for each taxonomic level
  taxonomic_sums_meta[[i]] <- sum(unlist(meta_clade_df), na.rm = TRUE)
  
  # Store taxonomic levels
  taxonomic_levels_meta[[i]] <- meta_clade_df
}


# Calculate sum of all reads
sum_domain_reads_meta <- sum(unlist(taxonomic_sums_meta), na.rm = TRUE)

# Calculate fractions
fractions_meta <- sapply(taxonomic_sums_meta, function(category_sums_meta) category_sums_meta/ taxonomic_sums_meta[[1]] * 100)
# Create a data frame with row names
fractions_df_meta <- data.frame(fractions_meta)
rownames(fractions_df_meta) <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

print(fractions_df_meta)

################## PLOT ########################################################

# Merge
merged_fractions = merge(fractions_df, fractions_df_meta, by = "row.names", all = TRUE)
colnames(merged_fractions)[colnames(merged_fractions)=="Row.names"]="Taxa"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions"]="All 18S"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions_met"]="Metazoa"

long_data = gather(merged_fractions, key = "Category", value = "Value", -Taxa)
taxa_order = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
long_data$Taxa <- factor(as.character(long_data$Taxa), levels = taxa_order)


p<-ggplot(long_data, aes(x = Taxa, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +
  labs(title = "Relative abundance of annotated reads per taxa level",
       x = "Taxa",
       y = "%") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = c("Metazoa" = "skyblue", "All 18S" = "orange"))  # Set colors

print(p)
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/reads_distribution.png",p)

############## NMDS & Shannon ###########################################################

# Similarity in sample composition 

matrix <- t(norm_seqtab_18S)

# Calculate Bray-Curtis dissimilarity
bray_matrix <- vegdist(matrix, method = "bray")

# Perform NMDS
nmds_result <- metaMDS(bray_matrix)

# Extract coordinates and Shannon diversity from the NMDS result
nmds_coordinates <- data.frame(nmds_result$points)
shannon_values <- diversity(matrix, index = "shannon")

# Combine NMDS coordinates and Shannon diversity values
nmds_data <- cbind(nmds_coordinates, Shannon_Diversity = shannon_values)

# Create a plot using ggplot2
library(ggplot2)
final_nmds <- ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Shannon_Diversity)) +
  geom_point() +
  scale_color_gradient(name = "Shannon Diversity", low = "black", high = "white") +
  ggtitle("NMDS plot showing shannon diversity for ASVs in each sample")+ 
  theme_bw()

# Print the plot
print(final_nmds)
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_shannon_samples_asvs.png",p)
