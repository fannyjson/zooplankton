# Comparison of reads from all 18S and only metazoa 18S metabarcoding data
# The files loaded are output files from read_merged_18S.R 

setwd("C:/Users/johan/OneDrive/R/Master project") 

#Taxa table

taxa_18S=read.table("C:/Users/johan/OneDrive/R/Master project/taxa_18S.tsv") 

#Normalized seqtab table metazoa only

norm_seqtab_18S=read.table("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_metazoa_240410.tsv")

#Normalized seqtab all 18S

norm_seqtab_18S_all=read.table("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_18S_240410.tsv")

#Original seqtab metazoa only

seqtab_18S=read.table("C:/Users/johan/OneDrive/R/Master project/seqtab_metazoa_240410.tsv")

#Original seqtab all 18S

seqtab_18S_all=read.table("C:/Users/johan/OneDrive/R/Master project/seqtab_18S_240410.tsv")

#Load libraries 

library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)

#Remove rows with NA on domain level

taxa_18S = taxa_18S[!is.na(taxa_18S[, 1]), ]

#Only keep ASVs found in both taxa_18S and seqtab

ix=rownames(seqtab_18S_all)

taxa_18S = taxa_18S[ix %in% rownames(taxa_18S), ]

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



########################## INVESTIGATE 18S DATA ################################

#Keep non-metazoa
taxa_18S = taxa_18S[!taxa_18S[, 4] == "Metazoa", ]

ix=rownames(seqtab_18S_all)
taxa_18S = taxa_18S[ix %in% rownames(taxa_18S), ]

#Number of total reads annotated for 18S

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


#Total reads for each taxa level

reads = lapply(1:9, function(i) sum(na.omit(clade_counts[[i]])))

reads_df = data.frame(
  clade = c("Domain", "Supergroup", "Division", "Subdivision", 
                  "Class", "Order", "Family", "Genus", "Species"),
  
  Reads = unlist(reads)
)

print(reads_df)


#Loop through clade counts to get the sum of reads for each taxonomic level 

taxonomic_levels = list()
taxonomic_sums =list()


for (i in 1:9) {
  
  clade = clade_counts[[i]]
  clade_df = na.omit(data.frame(clade))
  taxonomic_sums[[i]] = sum(unlist(clade_df), na.rm = TRUE)
  taxonomic_levels[[i]] = clade_df
}

#Calculate sum of all reads

sum_domain_reads=sum(unlist(taxonomic_sums), na.rm = TRUE)

#Calculate fractions

fractions = sapply(taxonomic_sums, function(category_sums) category_sums / taxonomic_sums[[1]] * 100)

#Create a data frame fractions for each taxonomic level 

fractions_df = data.frame(fractions)
rownames(fractions_df) = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
print(fractions_df)


#################### 18S (METAZOA ONLY) #######################################


#18S metazoa only

taxa_metazoa = taxa_18S[taxa_18S[, 4] == "Metazoa", ]
taxa_metazoa = taxa_metazoa[!is.na(taxa_metazoa[, 1]), ]

#Only keep taxa presemt om both seqtab and taxa table

ix=rownames(seqtab_18S)
taxa_metazoa = taxa_metazoa[ix %in% rownames(taxa_metazoa), ]

#Save metazoa taxa table 

write.table(taxa_metazoa, "C:/Users/johan/OneDrive/R/Master project/taxa_metazoa_processed.tsv", sep="\t")

#Unique taxa and ASVs found per taxlevel 

ASV_taxa = matrix(data = NA, 
                       nrow = 2,  
                       ncol = length(colnames(taxa_metazoa)),  
                       dimnames = list(c("#Taxa", "#ASV"), colnames(taxa_metazoa))) 

for (i in seq_along(colnames(taxa_metazoa))) {
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

#Normalized clade counts
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
meta_matr=meta_matr[rowSums(meta_matr !=0, na.rm=TRUE)>0 , , drop=FALSE]

#Total reads for each taxa level

reads = lapply(1:9, function(i) sum(na.omit(meta_clade_counts[[i]])))

reads_df = data.frame(
  clade = c("Domain", "Supergroup", "Division", "Subdivision", 
            "Class", "Order", "Family", "Genus", "Species"),
  
  Reads = unlist(reads)
)


print(reads_df)



#Loop through clade counts to get the sum of reads for each taxonomic level 

taxonomic_levels_meta = list()
taxonomic_sums_meta = list()

for (i in 1:9) {
  
  meta_clade = meta_clade_counts[[i]]
  
  meta_clade_df = na.omit(data.frame(meta_clade))
  
  taxonomic_sums_meta[[i]] = sum(unlist(meta_clade_df), na.rm = TRUE)

  taxonomic_levels_meta[[i]] = meta_clade_df
}


#Calculate fractions

fractions_meta = sapply(taxonomic_sums_meta, function(category_sums_meta) category_sums_meta/ taxonomic_sums_meta[[1]] * 100)

#Create a data frame fractions

fractions_df_meta = data.frame(fractions_meta)
rownames(fractions_df_meta) = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

print(fractions_df_meta)

################## Plot fractions  ########################################################

#Merge the fractopms for all 18S and metazoa only 

merged_fractions = merge(fractions_df, fractions_df_meta, by = "row.names", all = TRUE)
colnames(merged_fractions)[colnames(merged_fractions)=="Row.names"]="Taxa"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions"]="Non-metazoa"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions_meta"]="Metazoa"

merged_df = gather(merged_fractions, key = "Category", value = "Value", -Taxa)
taxa_order = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
merged_df$Taxa = factor(as.character(merged_df$Taxa), levels = taxa_order)


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
  scale_fill_manual(values = c("Metazoa" = 'lightblue', "Non-metazoa" = 'orange'))  

print(p)
ggsave = function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/reads_distribution_240516.png",p, width = 10, height = 6 )


############## NMDS & Shannon for ASVs (metazoa) ##############################

#Calculate similarity in ASV composition for samples for all metabarcoding samples

matrix = t(norm_seqtab_18S)

#Calculate Bray-Curtis dissimilarity

bray_matrix = vegdist(matrix, method = "bray")

#Perform NMDS

nmds_result = metaMDS(bray_matrix)

# Extract coordinates and Shannon diversity from the NMDS result

nmds_coordinates = data.frame(nmds_result$points)
shannon_values = diversity(matrix, index = "shannon")

#Combine NMDS coordinates and Shannon diversity values

nmds_data = cbind(nmds_coordinates, Shannon_Diversity = shannon_values)

windows()
p = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Shannon_Diversity)) +
  geom_point() +
  scale_color_gradient(name = "Shannon Diversity", low ="lightblue" , high = "darkblue") +  # Grayscale color scale
  ggtitle("NMDS plot showing Shannon Diversity for genera in each sample (PR2)") +
  theme_bw()  # Set theme to black and white (optional)

#Print the plot
print(p)

ggsave = function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)


#Save NMDS for ASVs

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_shannon_all_samples_asvs.png",final_nmds, width=10, height=6)

