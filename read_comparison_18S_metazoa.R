# Comparison of reads from 18S and only metazoan 18S 
# The files loaded are output files from read_merged_18S.R 

setwd("C:/Users/johan/OneDrive/R/Master project") 

#Taxa table
taxa_18S=read.table("C:/Users/johan/OneDrive/R/Master project/taxa_18S.tsv") 

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
clade_counts=list()

for (i in 1:ncol(taxa_18S)){
  matr=NULL
  clade=unique(taxa_18S[,i])
  clade=clade[!is.na(clade)]
  for (j in 1:length(clade)){
    ix=which(clade[j]==taxa_18S[,i])
    if (length(ix)>1){
      matr=rbind(matr, apply(norm_seqtab_18S_all[ix,],2,sum,na.rm=TRUE))
    } else {
      matr=rbind(matr, norm_seqtab_18S_all[ix,])
    }
  }
  rownames(matr)=clade
  #colnames(matr)=colnames(matr)=sample_id
  
  clade_counts[[i]]=matr
}


# Total reads for each taxa (18S)

domain_reads = clade_counts[[1]]
domain_reads = data.frame(domain_reads)
domain_reads = na.omit(domain_reads)
sum_domain_reads=sum(domain_reads)

supergroup_reads = clade_counts[[2]]
supergroup_reads = data.frame(supergroup_reads)
supergroup_reads = na.omit(supergroup_reads)
sum_supergroup_reads = sum(supergroup_reads)

division_reads = clade_counts[[3]]
division_reads = data.frame(division_reads)
division_reads = na.omit(division_reads)
sum_division_reads = sum(division_reads)

subdivision_reads = clade_counts[[4]]
subdivision_reads = data.frame(subdivision_reads)
subdivision_reads = na.omit(subdivision_reads)
sum_subdivision_reads = sum(subdivision_reads)

class_reads = clade_counts[[5]]
class_reads = data.frame(class_reads)
class_reads = na.omit(class_reads)
sum_class_reads = sum(class_reads)

order_reads = clade_counts[[6]]
order_reads = data.frame(order_reads)
order_reads = na.omit(order_reads)
sum_order_reads = sum(order_reads)

family_reads = clade_counts[[7]]
family_reads = data.frame(family_reads)
family_reads = na.omit(family_reads)
sum_family_reads = sum(family_reads)

genus_reads = clade_counts[[8]]
genus_reads = data.frame(genus_reads)
genus_reads = na.omit(genus_reads)
sum_genus_reads = sum(genus_reads)

species_reads = clade_counts[[9]]
species_reads = data.frame(species_reads)
species_reads = na.omit(species_reads)
species_reads_sum = sum(species_reads)

all_reads=sum_domain_reads+sum_supergroup_reads+sum_division_reads+sum_subdivision_reads+sum_class_reads+sum_order_reads+sum_family_reads+sum_genus_reads+sum_family_reads

taxonomic_levels = list(domain_reads, supergroup_reads, division_reads, subdivision_reads, class_reads, order_reads, family_reads, genus_reads, species_reads)

sums <- lapply(taxonomic_levels, function(df) sum(unlist(df), na.rm = TRUE))


all_reads_sum = colSums(do.call(cbind, sums), na.rm = TRUE)

# Calculate fractions
fractions = sapply(sums, function(category_sums) category_sums / sum_domain_reads * 100)
fractions_df=data.frame(fractions)
rownames(fractions_df)=c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
print(fractions)


#################### 18S (METAZOA ONLY) #######################################

#18S metazoa only
taxa_18S_metazoa = taxa_18S[!is.na(taxa_18S[, 4]), ]
taxa_metazoa = taxa_18S[taxa_18S[, 4] == "Metazoa", ]
#write.table(taxa_metazoa, "C:/Users/johan/OneDrive/R/Master project/taxa_metazoa")

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
meta_clade_counts=list()

for (i in 1:ncol(taxa_metazoa)){
  meta_matr=NULL
  clade=unique(taxa_metazoa[,i])
  clade=clade[!is.na(clade)]
  for (j in 1:length(clade)){
    ix=which(clade[j]==taxa_metazoa[,i])
    if (length(ix)>1){
      meta_matr=rbind(meta_matr, apply(norm_seqtab_18S[ix,],2,sum,na.rm=TRUE))
    } else {
      meta_matr=rbind(meta_matr, norm_seqtab_18S[ix,])
    }
  }
  
  rownames(meta_matr)=clade
  #colnames(matr)=colnames(norm_matr)=sample_id
  
  meta_clade_counts[[i]]=meta_matr
}


# Total reads for each taxa (18S)

domain_reads = meta_clade_counts[[1]]
domain_reads = data.frame(domain_reads)
domain_reads = na.omit(domain_reads)
sum_domain_reads=sum(domain_reads)

supergroup_reads = meta_clade_counts[[2]]
supergroup_reads = data.frame(supergroup_reads)
supergroup_reads = na.omit(supergroup_reads)
sum_supergroup_reads = sum(supergroup_reads)

division_reads = meta_clade_counts[[3]]
division_reads = data.frame(division_reads)
division_reads = na.omit(division_reads)
sum_division_reads = sum(division_reads)

subdivision_reads = meta_clade_counts[[4]]
subdivision_reads = data.frame(subdivision_reads)
subdivision_reads = na.omit(subdivision_reads)
sum_subdivision_reads = sum(subdivision_reads)

class_reads = meta_clade_counts[[5]]
class_reads = data.frame(class_reads)
class_reads = na.omit(class_reads)
sum_class_reads = sum(class_reads)

order_reads = meta_clade_counts[[6]]
order_reads = data.frame(order_reads)
order_reads = na.omit(order_reads)
sum_order_reads = sum(order_reads)

family_reads = meta_clade_counts[[7]]
family_reads = data.frame(family_reads)
family_reads = na.omit(family_reads)
sum_family_reads = sum(family_reads)

genus_reads = meta_clade_counts[[8]]
genus_reads = data.frame(genus_reads)
genus_reads = na.omit(genus_reads)
sum_genus_reads = sum(genus_reads)

species_reads = meta_clade_counts[[9]]
species_reads = data.frame(species_reads)
species_reads = na.omit(species_reads)
species_reads_sum = sum(species_reads)

all_reads=sum_domain_reads+sum_supergroup_reads+sum_division_reads+sum_subdivision_reads+sum_class_reads+sum_order_reads+sum_family_reads+sum_genus_reads+sum_family_reads


taxonomic_levels = list(domain_reads, supergroup_reads, division_reads, subdivision_reads, class_reads, order_reads, family_reads, genus_reads, species_reads)

sums <- lapply(taxonomic_levels, function(df) sum(unlist(df), na.rm = TRUE))


all_reads_sum = colSums(do.call(cbind, sums), na.rm = TRUE)

# Calculate fractions
fractions_met= sapply(sums, function(category_sums) category_sums / sum_domain_reads * 100)
fractions_df_met=data.frame(fractions_met)
rownames(fractions_df_met)=c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
print(fractions_met)


################## PLOT ########################################################

# Merge
merged_fractions = merge(fractions_df, fractions_df_met, by = "row.names", all = TRUE)
colnames(merged_fractions)[colnames(merged_fractions)=="Row.names"]="Taxa"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions"]="All 18S"
colnames(merged_fractions)[colnames(merged_fractions)=="fractions_met"]="Metazoa"

library(ggplot2)
library(tidyr)

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
