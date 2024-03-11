
#SET WD, LOAD FILES AND LIBRARIES

setwd("C:/Users/johan/OneDrive/R/Master project")

#Normalized seqtab table only common genera
new_seqtab_genus=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_common_genus.tsv", sep="\t")

#Merged df genus (all)
genus_meta_all=read.table("C:/Users/johan/OneDrive/R/Master project/merged_df_all.tsv", sep="\t")

#Merged df genus (only common with microscopy)
genus_meta_combined=read.table("C:/Users/johan/OneDrive/R/Master project/merged_df_common_genus.tsv", sep="\t")

#Count df microscopy
count_microscopy=read.table("C:/Users/johan/OneDrive/R/Master project/count_microscopy_commonsamples.tsv", sep="\t")

library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(ggplot2)


colnames(new_seqtab_genus) = gsub(colnames(new_seqtab_genus), pattern = '^X', replacement = '')
colnames(count_microscopy) = gsub(colnames(count_microscopy), pattern = '^X', replacement = '')
colnames(new_seqtab_genus) = gsub(colnames(new_seqtab_genus), pattern = '\\.', replacement = '-')
colnames(count_microscopy) = gsub(colnames(count_microscopy), pattern = '\\.', replacement = '-')

## Create dataframe with columns 'Sample', 'Genus', 'Microscopy counts', 'Reads', 'Rel abundance microscopy', 'Rel abundance reads' 

new_seqtab_genus = new_seqtab_genus[apply(new_seqtab_genus != 0, 1, any), , drop = FALSE]

count_microscopy$Genus = rownames(count_microscopy)
new_seqtab_genus$Genus = rownames(new_seqtab_genus)

df_microscopy_melted = count_microscopy %>%
  gather(key = "Sample", value = "Microscopy_count", -Genus)

df_reads_melted <- new_seqtab_genus %>%
  gather(key = "Sample", value = "Reads", -Genus)


merged_df <- full_join(df_microscopy_melted, df_reads_melted, by = c("Sample", "Genus"))
merged_df[is.na(merged_df)] <- 0


# Calculate relative abundance columns
merged_df <- merged_df %>%
  group_by(Sample) %>%
  mutate(Relative_abundance_microscopy = Microscopy_count / sum(Microscopy_count),
         Relative_abundance_reads = Reads / sum(Reads))


####################### BAR PLOTS (STATION) ##################################################


# Grouped stacked bar plot
my_df=merged_df

names(my_df)[names(my_df) == "Microscopy_count"] <- "Microscopy"
names(my_df)[names(my_df) == "Reads"] <- "Metabarcoding"

the_station='263856'
grep=grep(the_station, my_df$Sample)
my_df=my_df[grep,, drop=FALSE]
my_df=my_df %>%
  filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))

df_long <- pivot_longer(my_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy), 
                        names_to = "Method", values_to = "Relative_Abundance")

df_long <- df_long %>%
  mutate(counts = ifelse(Method == "Relative_Abundance", Metabarcoding, Microscopy))


df_long$Method <- ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                         ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding", df_long$Method))


df_long$Label <- paste(df_long$Sample, df_long$Method, sep = " - ")


df_long = subset(df_long, select = -c(Microscopy, Metabarcoding) )
df_long$Sample <- gsub(paste0(the_station, "_"), "", df_long$Sample)


myCol = c("#FFD700", "#FFA07A", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#FF6347", "#4CAF50","#FF9800","#673AB7", "#FF4081", "#9C27B0", "#03A9F4",'#02818a', "#3F51B5", "#CDDC39", "#3F51B5", '#e7298a', '#232356', 'darkred', '#fb8072', '#b2df8a', '#a63603', 'darkgreen', "#FFB6C1", '#fdbe6f', "#00FA9A", "#FFDAB9", "#D8BFD8", "#B0E0E6", "#20B2AA", "#9370DB", "#F0E68C", "#FFFACD", "#40E0D0", "#DA70D6","#B0C4DE", "#66CDAA",'#232356', 'blue','#f16913',"#D2B48C", "#CD853F","#BC8F8F",'#1f78b4', "#E6E6FA", "#DB7093", "#FF82AB", "#C0C0C0", "#C0D9D9", "#FDEE73", "#4682B4", "#B0E0E6", 'black', 'brown')

library(viridis)



# Plotting
p <- ggplot(df_long, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "RA2 (263856)",
       x = "", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 8)) +
  #scale_fill_viridis()
  scale_fill_manual(values = myCol)  # Use scale_fill_manual for fill colors

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_263856.png", p, width = 10, height = 5, bg = 'white')


################# SCATTER PLOTS ################################################

library(dplyr)
genus_counts <- merged_df %>%
  group_by(Genus) %>%
  summarize(counts = sum(Microscopy_count), reads = sum(Reads))

genus_counts= merged_df %>%
  

# Identify genera that appear more than 10 times
genera_to_keep = genus_counts %>%
  filter(counts > 20) %>%
  pull(Genus)


final_df <- merged_df %>%
  filter(Genus %in% genera_to_keep)

cor_values <- final_df %>%
  group_by(Genus) %>%
  summarize(correlation = cor(Relative_abundance_microscopy, Relative_abundance_reads))

cor_values=na.omit(cor_values)
final_df <- merge(final_df, cor_values, by = "Genus")

unique_correlations <- unique(final_df[, c("Genus", "correlation")])


# Plot
q <- ggplot(final_df, aes(x = Relative_abundance_microscopy, y = Relative_abundance_reads)) +
  geom_point(color = "black", size = 2, shape = 16) +
  geom_text(data = unique_correlations, 
            aes(label = sprintf("R = %.2f", correlation), x = 0.7, y = Inf, color = "blue"), 
            hjust = 1, vjust = 1, size = 3) +
  facet_wrap(~Genus, scales = "free") +
  ggtitle("Genera") +
  xlab("Relative abundance microscopy") +
  ylab("Relative abundance reads") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 14, face = "bold"))

print(q)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/scatter_genera.png", q, width = 30, height = 20, bg = 'white')