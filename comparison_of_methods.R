
#SET WD, LOAD FILES AND LIBRARIES

setwd("C:/Users/johan/OneDrive/R/Master project")

#Normalized seqtab table only common genera
new_seqtab_genus=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_common_genus.tsv", sep="\t")

#Normalized seqtab with all genera
new_seqtab=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_genus.tsv", sep="\t")
  
#Merged df
merged_metabar_df=read.table("C:/Users/johan/OneDrive/R/Master project/merged_metabar_df.tsv", sep="\t")

#Merged df genus (all)
all_merged_metabar_df=read.table("C:/Users/johan/OneDrive/R/Master project/merged_metabar_all_genera.tsv", sep="\t")

#Merged df genus (only common with microscopy)
genus_meta_combined=read.table("C:/Users/johan/OneDrive/R/Master project/merged_df_common_genus.tsv", sep="\t")

#Count df microscopy
count_microscopy=read.table("C:/Users/johan/OneDrive/R/Master project/count_microscopy_commonsamples.tsv", sep="\t")

library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(ggplot2)

myCol = c("#FFD700", "#FFA07A", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#FF6347", "#4CAF50","#FF9800","#673AB7", "#FF4081", "#9C27B0", "#03A9F4",'#02818a', "#3F51B5", "#CDDC39", "#3F51B5", '#e7298a', '#232356', 'darkred', '#fb8072', '#b2df8a', '#a63603', 'darkgreen', "#FFB6C1", '#fdbe6f', "#00FA9A", "#FFDAB9", "#D8BFD8", "#B0E0E6", "#20B2AA", "#9370DB", "#F0E68C", "#FFFACD", "#40E0D0", "#DA70D6","#B0C4DE", "#66CDAA",'#232356', 'blue','#f16913',"#D2B48C", "#CD853F","#BC8F8F",'#1f78b4', "#E6E6FA", "#DB7093", "#FF82AB", "#C0C0C0", "#C0D9D9", "#FDEE73", "#4682B4", "#B0E0E6", 'black', 'brown')

colnames(new_seqtab_genus) = gsub(colnames(new_seqtab_genus), pattern = '^X', replacement = '')
colnames(count_microscopy) = gsub(colnames(count_microscopy), pattern = '^X', replacement = '')
colnames(all_merged_metabar_df) = gsub(colnames(all_merged_metabar_df), pattern = '^X', replacement = '')
colnames(new_seqtab) = gsub(colnames(new_seqtab), pattern = '^X', replacement = '')

colnames(new_seqtab_genus) = gsub(colnames(new_seqtab_genus), pattern = '\\.', replacement = '-')
colnames(count_microscopy) = gsub(colnames(count_microscopy), pattern = '\\.', replacement = '-')
colnames(all_merged_metabar_df) = gsub(colnames(all_merged_metabar_df), pattern = '\\.', replacement = '-')
colnames(new_seqtab) = gsub(colnames(new_seqtab), pattern = '\\.', replacement = '-')

colnames(genus_meta_combined)[which(colnames(genus_meta_combined)=="sample")]="Sample"

#Rename months
rename_months <- data.frame(Original = c("januari", "februari", "mars", "april", "maj", "juni", "juli", "augusti", "september", "oktober", "november", "december"),
                                   Standardized = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"),
                                   stringsAsFactors = FALSE)

replace_months <- function(months) {
  mapping_row <- rename_months[rename_months$Original == months, ]
  if (nrow(mapping_row) > 0) {
    return(mapping_row$Standardized)
  } else {
    return(months)
  }
}

all_merged_metabar_df$month <- sapply(all_merged_metabar_df$month, replace_months)


## Create dataframe with columns 'Sample', 'Genus', 'Microscopy counts', 'Reads', 'Rel abundance microscopy', 'Rel abundance reads' 

new_seqtab = new_seqtab[apply(new_seqtab_genus != 0, 1, any), , drop = FALSE]

count_microscopy$Genus = rownames(count_microscopy)
new_seqtab$Genus = rownames(new_seqtab)

df_microscopy_melted = count_microscopy %>%
  gather(key = "Sample", value = "Microscopy_count", -Genus)

df_reads_melted <- new_seqtab %>%
  gather(key = "Sample", value = "Reads", -Genus)


merged_df=full_join(df_microscopy_melted, df_reads_melted, by = c("Sample", "Genus"))
merged_df[is.na(merged_df)] <- 0
merged_df$Station=substr(merged_df$Sample,1,6)
merged_df$Salinity=all_merged_metabar_df$Salinity[match(merged_df$Sample, all_merged_metabar_df$sample)]
merged_df$Month=all_merged_metabar_df$month[match(merged_df$Sample, all_merged_metabar_df$sample)]



####################### BAR PLOTS (STATION) ###################################################################

## stations=c("BY2 (263725)", "BY15 (263730)", "B1 (263738)", "BY31 (263732)", "REF M1V1 (263729)", "BY5 (263728)", "Å17 (264194)", "ANHOLT E (264876)", "N14 (263903)", "SLÄGGÖ (263628)", "GA1 (263757)", "C3 (190745)", "B7 (263629)", "RA2 (263856)")

# Grouped stacked bar plot
my_df=merged_df

# Calculate relative abundance columns
my_df <- my_df %>%
  group_by(Sample) %>%
  mutate(Relative_abundance_microscopy = Microscopy_count / sum(Microscopy_count),
         Relative_abundance_reads = Reads / sum(Reads))

names(my_df)[names(my_df) == "Microscopy_count"] <- "Microscopy"
names(my_df)[names(my_df) == "Reads"] <- "Metabarcoding"

the_station='264194'
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
df_long <- df_long %>%
  group_by(Month, Method, Sample, Genus, Station, Salinity) %>%
  summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
  filter(Relative_Abundance > 0.005)

library(viridis)

# Plotting
p <- ggplot(df_long, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Å17 (264194)",
       x = "", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 8)) +
  scale_fill_manual(values = myCol)

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_264194_0.005.png", p, width = 10, height = 5, bg = 'white')



################# SCATTER PLOT GENERA - NOT AVERAGED ################################################

merged_df <- merged_df %>%
  group_by(Sample) %>%
  mutate(Relative_abundance_microscopy = Microscopy_count / sum(Microscopy_count),
         Relative_abundance_reads = Reads / sum(Reads))

genus_counts <- merged_df %>%
  group_by(Genus) %>%
  summarize(counts = sum(Microscopy_count), reads = sum(Reads))

# Identify genera that appear more than 20 times
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

###################### SCATTER PLOT - MEAN RELATIVE ABUNANCE (ALL STATIONS) #######

merged_df=merged_df %>%
  filter(Genus %in% rownames(count_microscopy))

total_abundances=merged_df %>%
  group_by(Genus, Station) %>%
  summarise(Total_microscopy = sum(Relative_abundance_microscopy)/14,
            Total_metabarcoding = sum(Relative_abundance_reads)/14)


overall_correlation=total_abundances %>%
  group_by(Genus) %>%
  summarise(Correlation = cor(Total_microscopy, Total_metabarcoding))

overall_correlation=na.omit(overall_correlation)

q <- ggplot(total_abundances, aes(x =Total_microscopy , y = Total_metabarcoding)) +
  geom_point(color = "black", size = 2, shape = 16) +
  geom_text(data = overall_correlation, 
            aes(label = sprintf("R = %.2f", Correlation), x = Inf, y =Inf , color = "blue"), 
            hjust = 1, vjust = 1, size = 3) +
  facet_wrap(~Genus, scales = "free") +
  ggtitle("Genera") +
  xlab("Mean Relative abundance microscopy") +
  ylab("Mean Relative abundance reads") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 14, face = "bold"))

print(q)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/scatter_genera_avg123.png", q, width = 35, height = 25, bg = 'white')

################################## SCATTER PLOTS - ONE PLOT FOR EACH GENUS ###############
genera_list = unique(scatter_df$Genus)

# Loop through every genus, create a plot, and save it
for (genus in genera_list) {
  # Subset data for the current genus
  genus_data = subset(scatter_df, Genus == genus)
  
  # Plot
  q <- ggplot(genus_data, aes(x = Relative_abundance_microscopy, y = Relative_abundance_reads)) +
    geom_point(color = "black", size = 2, shape = 16) +
    geom_text(data = subset(unique_correlations, Genus == genus), 
              aes(label = sprintf("R = %.2f", Correlation), x = Inf, y =Inf , color = "blue"), 
              hjust = 1, vjust = 1, size = 3) +
    ggtitle(paste("Genus:", genus)) +
    xlab("Relative abundance microscopy") +
    ylab("Relative abundance reads") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 6),
          plot.title = element_text(size = 14, face = "bold")) +
    scale_x_continuous(limits = c(0, 1))
  
  # Save the plot with a unique name based on genus
  ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/scatter_", gsub(" ", "_", genus), "_avg.png"), q, width = 10, height = 10, bg = 'white')
}
###################### BAR PLOT - ONE BAR PER STATION (ALL YEARS/ONE YEAR) ##########################

#Only keep one year
#merged_df <- merged_df %>%
  #filter(str_detect(Sample, "2020"))

merged_df <- merged_df %>%
  group_by(Station) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy_count / sum(Microscopy_count),
    Relative_abundance_reads = Reads / sum(Reads)
  ) %>%
  filter(
    Relative_abundance_microscopy >= 0.01,
    Relative_abundance_reads >=0.01
  ) 

my_df=merged_df

names(my_df)[names(my_df) == "Microscopy_count"] <- "Microscopy"
names(my_df)[names(my_df) == "Reads"] <- "Metabarcoding"

#the_station='263856'
#grep=grep(the_station, my_df$Sample)
#my_df=my_df[grep,, drop=FALSE]
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


  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Community composition at common stations",
       x = "Station", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 8)) +
  scale_fill_manual(values = myCol)  # Use scale_fill_manual for fill colors

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_stations_0.01_both.png", p, width = 15, height = 7, bg = 'white')


################### VENN DIAGRAM MICROSCOPY AND METABARCODING ALL SAMPLES ######

genus_list_meta=unique(all_merged_metabar_df$Genus)
genus_list_micro=rownames(count_microscopy)
genus_list_common=intersect(genus_list_meta, genus_list_micro)

library(VennDiagram)

windows()
venn.plot <- venn.diagram(
  x = list(Metabarcoding = genus_list_meta, Microscopy = genus_list_micro), 
  category.names = c("Metabarcoding", "Microscopy"),
  filename = NULL,
  main = "Number of genera identified",
  col = "transparent", # Set circle colors to transparent to change them individually
  fill = c("#FF9800", "#87CEEB"), # Set fill colors for circles
  cat.pos = c(0, 0), # Set category positions to place names outside the circles
  cat.dist = c(0.1, 0.1), # Set distance of category names from the circles
  cat.fontfamily = "sans", # Set font family for category names
  cat.fontface = "bold", # Set font face for category names
  cat.col = c("#FF9800", "#87CEEB"), # Set colors for category names
  cat.cex = 1, # Set size of category names
  title.fontfamily = "sans", # Set font family for title
  title.fontface = "bold", # Set font face for title
  title.col = "black", # Set color for title
  title.cex = 1.5 # Set size of title
)

grid.draw(venn.plot)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/venn_diagram.png", venn.plot, width = 5, height = 5, bg = 'white')

############## TOTAL MICROSCOPY COUNTS IN COMMON SAMPLES #######################

count_microscopy[] <- lapply(count_microscopy, as.numeric)
sums_count = rowSums(count_microscopy, na.rm=TRUE)
total_count = data.frame(row.names = rownames(count_microscopy), Sum=sums_count )
total_count = total_count[rowSums(total_count) != 0, , drop=FALSE]
total_count = total_count %>%
  filter(Sum>100)


p<-ggplot(total_count, aes(x = reorder(row.names(total_count), -Sum), y = Sum)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Most common genera across all samples (count>100)",
       x = "Genus",
       y = "Counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10))

print(p)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/bar_count.png", p, width = 10, height = 5, bg = 'white')


###################### BAR PLOT - MONTHS ########################################


# Calculate relative abundance columns
months_df = merged_df %>%
  group_by(Month) %>%
  mutate(Relative_abundance_microscopy = Microscopy_count / sum(Microscopy_count),
         Relative_abundance_reads = Reads / sum(Reads))

names(months_df)[names(months_df) == "Microscopy_count"] <- "Microscopy"
names(months_df)[names(months_df) == "Reads"] <- "Metabarcoding"

months_df=months_df %>%
  filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))

df_long <- pivot_longer(months_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy), 
                        names_to = "Method", values_to = "Relative_Abundance")

df_long <- df_long %>%
  mutate(counts = ifelse(Method == "Relative_Abundance", Metabarcoding, Microscopy))


df_long$Method <- ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                         ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding", df_long$Method))


df_long$Label <- paste(df_long$Sample, df_long$Method, sep = " - ")


df_long = subset(df_long, select = -c(Microscopy, Metabarcoding) )

df_long$Month <- factor(df_long$Month, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))

df_long <- df_long %>%
  group_by(Month, Method, Genus, Station, Salinity) %>%
  summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
  filter(Relative_Abundance > 0.005)

p <- ggplot(df_long, aes(x = Month, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Temporal variation of microbial composition across samples",
       x = "Month", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 8)) +
  scale_fill_manual(values = myCol)  # Use scale_fill_manual for fill colors

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_months_0.005.png", p, width = 15, height = 7, bg = 'white')


############## NMDS FOR STATIONS (NOT DONE)####################################################################

merged_df$Station <- as.numeric(merged_df$Station)

# Create separate dataframes for metabarcoding and microscopy data
metabarcoding_data <- merged_df[, c("Sample", "Genus", "Station", "Reads", "Month")]
microscopy_data <- merged_df[, c("Sample", "Genus", "Station", "Microscopy_count", "Month")]

month_names <- c("January", "February", "March", "April", "May", "June", 
                 "July", "August", "September", "October", "November", "December")

# Create a mapping between month names and numeric values
month_mapping <- setNames(1:12, month_names)

metabarcoding_data <- metabarcoding_data[metabarcoding_data$Reads > 0, ]

bray_curtis_metabarcoding <- vegdist(metabarcoding_data[, c("Reads", "Station")], method = "bray")

# Run NMDS for metabarcoding data
metabarcoding_nmds <- metaMDS(bray_curtis_metabarcoding)

# Convert NMDS results to data frame
metabarcoding_plot_data <- data.frame(metabarcoding_nmds$points, Station = metabarcoding_data$Station)

# Plot NMDS ordination for metabarcoding data using ggplot2
ggplot(metabarcoding_plot_data, aes(x = NMDS1, y = NMDS2, color = Genus, shape = as.factor(Station))) +
  geom_point() +
  labs(color = "Genus", shape = "Station") +
  theme_minimal() +
  ggtitle("NMDS Plot for Metabarcoding Data (Bray-Curtis)")

# Assessing the NMDS results
stress_value <- metabarcoding_nmds$stress
print(paste("Stress value:", stress_value))





