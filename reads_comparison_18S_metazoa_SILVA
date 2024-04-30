# Comparison of reads from 18S and only metazoan 18S 
# The files loaded are output files from read_merged_18S.R 

setwd("C:/Users/johan/OneDrive/R/Master project") 

#Taxa table

taxa_18S=read.table("C:/Users/johan/OneDrive/R/Master project/taxa_SILVA.tsv") 

#Normalized seqtab table metazoa 

norm_seqtab_18S=read.table("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_metazoa_silva_processed.tsv")

#Normalized seqtab all 

norm_seqtab_18S_all=read.table("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_18S_silva_processed.tsv")

#Original seqtab metazoa

seqtab_18S=read.table("C:/Users/johan/OneDrive/R/Master project/seqtab_metazoa_silva_processed.tsv")

#Original seqtab all

seqtab_18S_all=read.table("C:/Users/johan/OneDrive/R/Master project/seqtab_18S_silva_processed.tsv")

#Load libraries 
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)

#################################################################################

#Fix sample names 

colnames(norm_seqtab_18S) = gsub(colnames(norm_seqtab_18S), pattern = '^X', replacement = '')
colnames(norm_seqtab_18S_all) = gsub(colnames(seqtab_18S_all), pattern = '^X', replacement = '')

colnames(seqtab_18S) = gsub(colnames(seqtab_18S), pattern = '^X', replacement = '')
colnames(seqtab_18S_all) = gsub(colnames(seqtab_18S_all), pattern = '^X', replacement = '')

replace_xx = function(vector) {
  gsub("(_X|_XX)", "", vector)
}

taxa_18S[] = lapply(taxa_18S, replace_xx)

replace = function(vector) {
  gsub("_", " ", vector)
}


taxa_18S[] = lapply(taxa_18S, replace)

########################## #18S ALL ############################################

#Only keep ASVs found in both taxatable and seqtab

ix=rownames(seqtab_18S_all)
taxa_18S = taxa_18S[ix %in% rownames(taxa_18S), ]


#Number of reads for 18S (annotated)

all_reads=sum(seqtab_18S_all)
print(all_reads)

#Summarising counts at different taxonomic levels
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
  clade_counts[[i]]=matr
}

matr[is.na(matr)]=0
matr=matr[rowSums(matr !=0, na.rm=TRUE)>0 , , drop=FALSE]


# Total reads for each taxa (18S)


reads = lapply(1:7, function(i) sum(na.omit(clade_counts[[i]])))

reads_df = data.frame(
  clade = c("Domain", "Kingdom", "Class", "Order", 
             "Family", "Genus", "Species"),
  Reads = unlist(reads)
)

print(reads_df)

#Loop through clade counts

taxonomic_levels = list()
taxonomic_sums = list()

for (i in 1:7) {
  # Extract clade counts
  clade = clade_counts[[i]]
  
  # Convert to data frame and remove NAs
  clade_df = na.omit(data.frame(clade))
  
  # Calculate sum for each taxonomic level
  taxonomic_sums[[i]] = sum(unlist(clade_df), na.rm = TRUE)
  
  # Store taxonomic levels
  taxonomic_levels[[i]] = clade_df
}

# Calculate fractions

fractions = sapply(taxonomic_sums, function(category_sums) category_sums / taxonomic_sums[[1]] * 100)

#Make dataframe with fractions

fractions_df = data.frame(fractions)
rownames(fractions_df) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

print(fractions_df)


#################### 18S (METAZOA ONLY) #######################################

#18S metazoa only

taxa_metazoa = taxa_18S[taxa_18S[, 4] == "Metazoa (Animalia)", ]

taxa_metazoa = taxa_metazoa[!is.na(taxa_metazoa[, 1]), ]
colnames(taxa_metazoa) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

write.table(taxa_metazoa, "C:/Users/johan/OneDrive/R/Master project/taxa_metazoa_SILVA.tsv", sep="\t")

#Count unique ASVs and taxlevels

ASV_taxa = matrix(data = NA, 
                   nrow = 2,  
                   ncol = length(colnames(taxa_metazoa)),  
                   dimnames = list(c("#Taxa", "#ASV"), colnames(taxa_metazoa))) 

for (i in seq_along(colnames(taxa_metazoa))) {
  # Count unique taxa for the current column
  taxa = length(unique(taxa_metazoa[!is.na(taxa_metazoa[, i]), i]))
  ASV_taxa[1, i] = taxa
  
  # Count unique ASVs for the current column
  asvs = length(unique(rownames(taxa_metazoa[!is.na(taxa_metazoa[, i]), , drop = FALSE])))
  ASV_taxa[2, i] = asvs
}

print(ASV_taxa)

#Number of reads for metazoa

metazoa_reads=sum(seqtab_18S)

print(metazoa_reads)

#Non-normalized clade count for metazoa

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

#Reads at each taxonomic level

reads = lapply(1:7, function(i) sum(na.omit(meta_clade_counts[[i]])))

reads_df = data.frame(
  clade = c("Kingdom", "Phylum", "Class", "Order", 
            "Family", "Genus", "Species"),
  Reads = unlist(reads)
)


print(reads_df)


#Loop through meta_clade_counts and store the sums

taxonomic_levels_meta = list()
taxonomic_sums_meta = list()

for (i in 1:7) {
  # Extract clade counts
  clade_meta = meta_clade_counts[[i]]
  
  # Convert to data frame and remove NAs
  clade_df_meta = na.omit(data.frame(clade_meta))
  
  # Calculate sum for each taxonomic level
  taxonomic_sums_meta[[i]] = sum(unlist(clade_df_meta), na.rm = TRUE)
  
  # Store taxonomic levels
  taxonomic_levels_meta[[i]] = clade_df_meta
}

# Calculate sum of all reads
sum_domain_reads_meta = sum(unlist(taxonomic_sums_meta), na.rm = TRUE)

# Calculate fractions
fractions_meta = sapply(taxonomic_sums_meta, function(category_sums_meta) category_sums_meta / taxonomic_sums_meta[[1]] * 100)

#Create dataframe with fractions
fractions_df_meta = data.frame(fractions_meta)
rownames(fractions_df_meta) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

print(fractions_df_meta)


################## PLOT ########################################################

#Merge 18S fractions and metazoa only fractions

merged_fractions = merge(fractions_df, fractions_df_meta, by = "row.names", all = TRUE)
colnames(merged_fractions)[colnames(merged_fractions)=="Row.names"]="Taxa"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions"]="All 18S"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions_meta"]="Metazoa"


merged_df = gather(merged_fractions, key = "Category", value = "Value", -Taxa)
taxa_order = c("Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species")
merged_df$Taxa = factor(as.character(merged_df$Taxa), levels = taxa_order)

windows()
p=ggplot(merged_df, aes(x = Taxa, y = Value, fill = Category)) +
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
  scale_fill_manual(values = c("Metazoa" = '#636363', "All 18S" = '#bdbdbd'))
print(p)
ggsave = function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/reads_distribution_SILVA.png",p, width = 10, height = 6)


############## NMDS & Shannon for metazoa ASVs ###########################################################

matrix = t(norm_seqtab_18S)

#Calculate Bray-Curtis dissimilarity
bray_matrix = vegdist(matrix, method = "bray")

#Perform NMDS
nmds_result = metaMDS(bray_matrix)

# Extract coordinates and Shannon diversity from the NMDS result
nmds_coordinates = data.frame(nmds_result$points)
shannon_values = diversity(matrix, index = "shannon")

# Combine NMDS coordinates and Shannon diversity values
nmds_data = cbind(nmds_coordinates, Shannon_Diversity = shannon_values)

windows()
p = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Shannon_Diversity)) +
  geom_point() +
  scale_color_gradient(name = "Shannon Diversity", low = "lightblue" , high = "darkblue") +  # Grayscale color scale
  ggtitle("NMDS plot showing Shannon Diversity for ASVs in each sample (SILVA)") +
  theme_bw()  # Set theme to black and white (optional)

print(p)

ggsave = function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_shannon_all_samples_SILVA.png",p, width = 10, height = 6 )
