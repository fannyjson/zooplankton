
#SET WD, LOAD FILES AND LIBRARIES

setwd("C:/Users/johan/OneDrive/R/Master project")

#Normalized seqtab table only common genera
new_seqtab_genus=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_common_genera_processed.tsv", sep="\t")

#Normalized seqtab with all genera
new_seqtab=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_genus_processed.tsv", sep="\t")

#Merged df genus (all)
all_merged_metabar_df=read.table("C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_PR2_processed.tsv", sep="\t")

#Merged df genus (only common with microscopy)
genus_meta_combined=read.table("C:/Users/johan/OneDrive/R/Master project/merged_metabar_common_genera_processed.tsv", sep="\t")

#Count df microscopy
count_microscopy=read.table("C:/Users/johan/OneDrive/R/Master project/count_microscopy_commonsamples.tsv", sep="\t")

#Taxa metazoa 
taxa=read.delim("C:/Users/johan/Downloads/microscopy_taxa.tsv")


library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(ggplot2)

myCol = c("#FFD700", "#FFA07A", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#FF6347", "#4CAF50","#FF9800","#673AB7", "#FF4081", "#9C27B0", "#03A9F4",'#02818a', "#CDDC39", "#3F51B5", '#e7298a', '#232356', 'darkred', '#fb8072', '#b2df8a', '#a63603', 'darkgreen', "#FFB6C1", '#fdbe6f', "#00FA9A", "#FFDAB9", "#D8BFD8", "#B0E0E6", "#20B2AA", "#9370DB", "#F0E68C", "#FFFACD", "#40E0D0", "#DA70D6","#B0C4DE", "#66CDAA",'#232356', 'blue','#f16913',"#D2B48C", "#CD853F","#BC8F8F",'#1f78b4', "#E6E6FA", "#DB7093", "#FF82AB", "#C0C0C0", "#C0D9D9", "#FDEE73", "#4682B4", "#B0E0E6", 'black'
          ,"limegreen","hotpink","brown","turquoise", "#012345")

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
rename_months = data.frame(Original = c("januari", "februari", "mars", "april", "maj", "juni", "juli", "augusti", "september", "oktober", "november", "december"),
                                   Standardized = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"),
                                   stringsAsFactors = FALSE)

replace_months = function(months) {
  mapping_row = rename_months[rename_months$Original == months, ]
  if (nrow(mapping_row) > 0) {
    return(mapping_row$Standardized)
  } else {
    return(months)
  }
}

all_merged_metabar_df$month = sapply(all_merged_metabar_df$month, replace_months)

#Create dataframe with columns 'Sample', 'Genus', 'Microscopy counts', 'Reads',
#'Relative abundance microscopy', 'Relative abundance reads' 

count_microscopy$Genus = rownames(count_microscopy)
new_seqtab$Genus = rownames(new_seqtab)

df_microscopy_melted = count_microscopy %>%
  gather(key = "Sample", value = "Microscopy_count", -Genus)

df_reads_melted = new_seqtab %>%
  gather(key = "Sample", value = "Reads", -Genus)


merged_df=full_join(df_microscopy_melted, df_reads_melted, by = c("Sample", "Genus"))
merged_df[is.na(merged_df)] = 0
merged_df$Station=substr(merged_df$Sample,1,6)
merged_df$Salinity=all_merged_metabar_df$Salinity[match(merged_df$Sample, all_merged_metabar_df$sample)]
merged_df$Month=all_merged_metabar_df$month[match(merged_df$Sample, all_merged_metabar_df$sample)]

names(merged_df)[names(merged_df) == "Microscopy_count"] = "Microscopy"
names(merged_df)[names(merged_df) == "Reads"] = "Metabarcoding"
merged_df=na.omit(merged_df)

#Use for comparison between stations/months
stations_df=merged_df
months_df=merged_df

#Use for other plots
merged_df = merged_df %>%
  group_by(Sample) %>%
  mutate(Relative_abundance_microscopy = Microscopy / sum(Microscopy),
         Relative_abundance_reads = Metabarcoding / sum(Metabarcoding))

merged_df$Relative_abundance_reads[is.nan(merged_df$Relative_abundance_reads)]=0


#Check what genera is unique for metabarcoding
metabarcoding_genera=merged_df %>%
  filter(Genus %in% Genus[Relative_abundance_reads > 0] & 
           !(Genus %in% Genus[Relative_abundance_microscopy > 0]))

unique_genera=unique(metabarcoding_genera$Genus)
print(unique_genera)


####################### BAR PLOT - ONE PER STATION ##################################################

#Station IDs and names

unique_stations=unique(merged_df$Station)
print(unique_stations)

merged_df = merged_df %>%
  group_by(Sample) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding)
  ) %>% #Use for filtering out relative abundances >=0.5% in each sample for both methods
  filter(
    Relative_abundance_microscopy >= 0.005,
    Relative_abundance_reads >=0.005
  ) 


station_ids = c("264194","263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")

station_names = c("Å17 (264194)", "BY31 (263732)", "BY2 (263725)", "B1 (263738)", "REF M1V1 (263729)", "BY5 (263728)", "ANHOLT E (264876)", "N14 (263903)", "SLÄGGÖ (263628)", "GA1 (263757)","C3 (190745)", "B7 (263629)", "RA2 (263856)", "BY15 (263730)")

#Loop over each station
for (i in seq_along(station_ids)) {
  the_station = station_ids[i]
  station_name = station_names[i]
  
  #Filter data for the station
  grep = grep(the_station, merged_df$Sample)
  filtered_df = merged_df[grep,, drop = FALSE]
  
  #Remove rows with 0 abundance for both methods
  filtered_df = filtered_df %>%
    filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))
  
  #Pivot data
  df_long = pivot_longer(filtered_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy), 
                          names_to = "Method", values_to = "Relative_Abundance")
  
  #Add counts column
  df_long = df_long %>%
    mutate(counts = ifelse(Method == "Relative_Abundance", Metabarcoding, Microscopy))
  
  # Rename Method values
  df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                           ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding", df_long$Method))
  
  #Label
  df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")
  
  #Calculate new relative abundances for each sample
  df_long = df_long %>%
    group_by(Month, Method, Sample, Genus, Station, Salinity) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance))
  
  #Run when including everything
  # Label genera with low abundance as "Other genera"
  #df_long = df_long %>%
    #mutate(Genus = ifelse(Relative_Abundance > 0.05, as.character(Genus), "Other genera"))

  
  # Plotting
  p = ggplot(df_long, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Method, scales = "free_y") +
    labs(title = paste("Comparison", station_name, sep = " - "),
         x = "", y = "Relative abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 8)) +
    scale_fill_manual(values = myCol)  # Use scale_fill_manual for fill colors
  
  # Print the plot
  print(p)
  
  # Save the plot
  filename = paste("C:/Users/johan/OneDrive/R/Master project/plots/comparison_", the_station, "_5.png", sep = "")
  ggsave(filename, p, width = 15, height = 6, bg = 'white')


}

################# SEPARATE SCATTER PLOTS RELATIVE ABUNDACE SAMPLES ################################################

genera_to_keep=c("Oithona", "Calanus", "Pseudocalanus", "Acartia", "Eurytemora", "Paracalanus", "Temora", "Centropages")

final_df = merged_df %>%
  filter(Genus %in% genera_to_keep)

cor_values = final_df %>%
  group_by(Genus) %>%
  summarize(Correlation = cor(Relative_abundance_microscopy, Relative_abundance_reads))

cor_values = na.omit(cor_values)
final_df = merge(final_df, cor_values, by = "Genus")

unique_correlations = unique(final_df[, c("Genus", "Correlation")])

for (genus in genera_to_keep) {
  # Subset data for the current genus
  genus_data = subset(final_df, Genus == genus)
  
  # Plot
  q = ggplot(genus_data, aes(x = Relative_abundance_microscopy, y = Relative_abundance_reads)) +
    geom_point(color = "black", size = 2, shape = 16) +
    ggtitle(paste(genus)) +
    xlab("Relative abundance microscopy") +
    ylab("Relative abundance reads") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 20, face = "bold")) +
    scale_x_continuous(limits = c(0, max(genus_data$Relative_abundance_microscopy))) + 
    scale_y_continuous(limits = c(0, max(genus_data$Relative_abundance_reads)))+ 
    geom_smooth(method = "lm", se = TRUE) + 
    geom_ribbon(aes(ymin = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,2],
                    ymax = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,3]),
                fill = "gray", alpha = 0.2)   # Add shaded area around the trend line
    
  # Add Pearson value to plot
  corr_value = unique_correlations$Correlation[which(unique_correlations$Genus == genus)]
  q = q +
    geom_text(aes(label = paste0("R = ", round(corr_value, 2))), 
              x = Inf, y = Inf, color = "blue", hjust = 1, vjust = 1, size = 3)
  
  # Save the plot with a unique name based on genus
    ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/scatter_", gsub(" ", "_", genus), "_avg_samples.png"), q, width = 10, height = 8, bg = 'white')
}


###################### SCATTER PLOT - MEAN RELATIVE ABUNDANCE (ALL STATIONS) #######

# Filter the merged dataframe
scatter_df = merged_df %>%
  filter(Genus %in% rownames(count_microscopy))

# Calculate total abundances
total_abundances = scatter_df %>%
  group_by(Genus, Station) %>%
  summarise(Total_microscopy = sum(Relative_abundance_microscopy) / 14,
            Total_metabarcoding = sum(Relative_abundance_reads) / 14)

# Calculate correlations
overall_correlation = total_abundances %>%
  group_by(Genus) %>%
  summarise(Correlation = cor(Total_microscopy, Total_metabarcoding, use = "pairwise.complete.obs"))

# Remove rows with NA values
overall_correlation = na.omit(overall_correlation)

# Create a plot
q = ggplot(total_abundances, aes(x = Total_microscopy, y = Total_metabarcoding)) +
  geom_point(color = "black", size = 2, shape = 16) +
  facet_wrap(~ Genus, scales = "free") +
  ggtitle("Genera") +
  xlab("Mean Relative abundance microscopy") +
  ylab("Mean Relative abundance reads") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold")) +
  scale_x_continuous(limits = c(0, max(total_abundances$Total_microscopy))) + 
  scale_y_continuous(limits = c(0, max(total_abundances$Total_metabarcoding)))+ 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_ribbon(aes(ymin = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,2],
                  ymax = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,3]),
              fill = "gray", alpha = 0.2)   # Add shaded area around the trend line

#Add Pearson values to plots
q = q + 
  geom_text(data = overall_correlation, 
            aes(label = paste0("R = ", round(Correlation, 2)), 
                x = Inf, y = Inf, color = "blue"), 
            hjust = 1, vjust = 1, size = 3)

# Print the plot
print(q)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/scatter_genera_240415.png", q, width = 35, height = 25, bg = 'white')


################################## FAMILY - GENUS - PEARSON table to report ###############################

# Trim whitespace from genus column in family_genus_df
family_genus_df$Genus = trimws(family_genus_df$Genus)

#Trim whitespace from genus column in taxa DataFrame
taxa$Genus = trimws(taxa$Genus)

#Perform left join to retain all genus names from family_genus_df
merge_family_genus = left_join(family_genus_df, taxa[, c("Genus", "Family")], by = "Genus")

# Keep only unique genus in merge_family_genus
merge_family_genus = merge_family_genus[!duplicated(merge_family_genus$Genus), ]
merge_family_genus = merge_family_genus[, c("Family", setdiff(names(merge_family_genus), "Family"))]

library(xtable)

italic_function = function(x) {
  paste("\\textit{", x, "}", sep = "")
}


latex_table = xtable(merge_family_genus)


print(latex_table, sanitize.text.function = italic_function, include.rownames = FALSE)


################################## SCATTER PLOT - SAVE SEPARATELY - MEAN RELATIVE ABUNDANCE STATION  ###############

genus_spec=c("Oithona", "Calanus", "Pseudocalanus", "Acartia", "Eurytemora", "Paracalanus", "Temora", "Centropages")

# Loop through each genus
for (genus in genus_spec) {
  # Filter data for the current genus
  genus_df = total_abundances %>%
    filter(Genus == genus)

  q = ggplot(genus_df, aes(x = Total_microscopy, y = Total_metabarcoding)) +
    geom_point(color = "black", size = 5, shape = 16) +
    ggtitle(paste(genus)) +
    xlab("Mean relative abundance microscopy") +
    ylab("Mean relative abundance reads") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 20, face = "bold")) +
    scale_x_continuous(limits = c(0, max(genus_df$Total_microscopy))) + 
    scale_y_continuous(limits = c(0, max(genus_df$Total_metabarcoding)))+ 
    geom_smooth(method = "lm", se = TRUE) + 
    geom_ribbon(aes(ymin = predict(lm(Total_metabarcoding ~ Total_microscopy, data = genus_df), interval = "confidence")[,2],
                    ymax = predict(lm(Total_metabarcoding ~ Total_microscopy, data = genus_df), interval = "confidence")[,3]),
                fill = "gray", alpha = 0.2)   # Add shaded area around the trend line
  
  # Add Pearson value to plot
  corr_value = overall_correlation$Correlation[which(overall_correlation$Genus == genus)]
  q = q +
    geom_text(aes(label = paste0("R = ", round(corr_value, 2))), 
              x = Inf, y = Inf, color = "blue", hjust = 1, vjust = 1, size = 3)

  
  ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/", "scatter_", genus, "_avgstation.png"), q, width = 10, height = 8, bg = 'white')
}



genera_list = unique(scatter_df$Genus)

# Loop through each genus, create a plot, and save it
for (genus in genera_list) {
  # Subset data for the current genus
  genus_data = subset(scatter_df, Genus == genus)
  
  # Plot
  q = ggplot(genus_data, aes(x = Relative_abundance_microscopy, y = Relative_abundance_reads)) +
    geom_point(color = "black", size = 2, shape = 16) +
    geom_text(data = subset(unique_correlations, Genus == genus), 
              aes(label = sprintf("R = %.2f", Correlation), x = Inf, y =Inf , color = "blue"), 
              hjust = 1, vjust = 1, size = 3) +
    ggtitle(paste("Genus:", genus)) +
    xlab("Relative abundance microscopy") +
    ylab("Relative abundance reads") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 20, face = "bold")) +
    scale_x_continuous(limits = c(0, max(genus_data$Relative_abundance_microscopy))) + 
    scale_y_continuous(limits = c(0, max(genus_data$Relative_abundance_reads)))+ 
    geom_smooth(method = "lm", se = FALSE) + 
    geom_ribbon(aes(ymin = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,2],
                    ymax = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,3]),
                fill = "gray", alpha = 0.2)   # Add shaded area around the trend line
  
  # Save the plot with a unique name based on genus
  ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/scatter_", gsub(" ", "_", genus), "_avg.png"), q, width = 10, height = 10, bg = 'white')
}
###################### BAR PLOT (ONE BAR PER STATION, ALL YEARS) ##########################

#Only keep one year
#merged_df = merged_df %>%
  #filter(str_detect(Sample, "2020"))


stations_df = stations_df %>%
  group_by(Station) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding)
  ) %>%
  filter(
    Relative_abundance_microscopy >= 0.005,
    Relative_abundance_reads >=0.005
  ) 



#the_station='263856'
#grep=grep(the_station, my_df$Sample)
#my_df=my_df[grep,, drop=FALSE]
#stations_df=stations_df %>%
  #filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))

df_long = pivot_longer(stations_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy), 
                        names_to = "Method", values_to = "Relative_Abundance")

df_long = df_long %>%
  mutate(counts = ifelse(Method == "Relative_Abundance", Metabarcoding, Microscopy))


df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                         ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding", df_long$Method))


df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")


df_long = subset(df_long, select = -c(Microscopy, Metabarcoding) )


df_long = df_long %>%
  group_by(Month, Method, Genus, Station, Salinity) %>%
  summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
  filter(Relative_Abundance > 0.005)

#df_long = df_long %>%
  #mutate(Genus = ifelse(Relative_Abundance > 0.005, as.character(Genus), "Other genera"))

p = ggplot(df_long, aes(x = reorder(Station, Salinity), y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Spatial variation of microbial composition across stations",
       x = "Station", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 8)) +
  scale_fill_manual(values = myCol)  # Use scale_fill_manual for fill colors

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_stations_240415_0.5.png", p, width = 15, height = 7, bg = 'white')


################### VENN DIAGRAM MICROSCOPY AND METABARCODING ALL SAMPLES ######

genus_list_meta=unique(all_merged_metabar_df$Genus)
genus_list_micro=rownames(count_microscopy)
genus_list_common=intersect(genus_list_meta, genus_list_micro)

library(VennDiagram)

windows()
venn.plot = venn.diagram(
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
count_microscopy[] = lapply(count_microscopy, as.numeric)
sums_count = rowSums(count_microscopy, na.rm=TRUE)
total_count = data.frame(row.names = rownames(count_microscopy), Sum=sums_count )
total_count = total_count[rowSums(total_count) != 0, , drop=FALSE]
total_count = total_count %>%
  filter(Sum>100)


p=ggplot(total_count, aes(x = reorder(row.names(total_count), -Sum), y = Sum)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Most common genera across all samples (count>100)",
       x = "Genus",
       y = "Counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10))+


print(p)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/bar_count.png", p, width = 10, height = 5, bg = 'white')

###################### BAR PLOT MONTH ########################################

# Calculate relative abundance columns
months_df = months_df %>%
  group_by(Month) %>%
  mutate(Relative_abundance_microscopy = Microscopy/sum(Microscopy),
         Relative_abundance_reads = Metabarcoding/sum(Metabarcoding))

months_df=months_df %>%
  filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))

months_df=na.omit(months_df)

df_long = pivot_longer(months_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy), 
                        names_to = "Method", values_to = "Relative_Abundance")

df_long = df_long %>%
  mutate(counts = ifelse(Method == "Relative_Abundance", Metabarcoding, Microscopy))


df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                         ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding", df_long$Method))


df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")


df_long = subset(df_long, select = -c(Microscopy, Metabarcoding) )

df_long$Month = factor(df_long$Month, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))

#df_long = df_long %>%
  #group_by(Month, Method, Genus, Station, Salinity) %>%
  #summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
  #filter(Relative_Abundance > 0.005)

df_long = df_long %>%
mutate(Genus = ifelse(Relative_Abundance > 0.005, as.character(Genus), "Other genera"))

df_long
p = ggplot(df_long, aes(x = Month, y = Relative_Abundance, fill = Genus)) +
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
  scale_fill_manual(values = myCol) 

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_months_0.5_all.png", p, width = 15, height = 7, bg = 'white')



