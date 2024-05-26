
#SET WD, LOAD FILES AND LIBRARIES

setwd("C:/Users/johan/OneDrive/R/Master project")


#Normalized seqtab table only common genera
#PR2
new_seqtab_genus=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_common_genera.tsv", sep="\t")


#Normalized seqtab with all genera
#PR2
new_seqtab=read.table("C:/Users/johan/OneDrive/R/Master project/new_seqtab_genus.tsv", sep="\t")

#Merged df genus (all)
#PR2
all_merged_metabar_df=read.table("C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_PR2.tsv", sep="\t")


#Merged df genus (only common with microscopy)
#PR2
genus_meta_combined=read.table("C:/Users/johan/OneDrive/R/Master project/merged_metabar_common_genera.tsv", sep="\t")

#Count df microscopy
#PR2
count_microscopy=read.table("C:/Users/johan/OneDrive/R/Master project/count_microscopy_commonsamples.tsv", sep="\t")

#Actual and predicted genera from random forest using 5-fold cross-validation 
rf_data=read.delim("C:/Users/johan/OneDrive/R/ML/RF/sample_predicted.tsv")

#Not cross-validated

#rf_data=read.delim("C:/Users/johan/OneDrive/R/ML/RF/sample_predicted_240504.tsv")
#rf_data$Predicted[rf_data$Predicted<1] = 0

#Taxa metazoa 
taxa=read.delim("C:/Users/johan/Downloads/microscopy_taxa.tsv")

#Intersect
match_genus=intersect(rf_data$Genus, all_merged_metabar_df$Genus)

library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(ggplot2)
library(VennDiagram)

myCol = c("#FFD700", "#FFA07A", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#FF6347", "#4CAF50","#FF9800","#673AB7", "#FF4081", "#9C27B0", "#03A9F4",'#02818a', "#CDDC39", "#3F51B5", '#e7298a', '#232356', 'darkred', '#fb8072', '#b2df8a', '#a63603', 'darkgreen', "#FFB6C1", '#fdbe6f', "#00FA9A", "#FFDAB9", "#D8BFD8", "#B0E0E6", "#20B2AA", "#9370DB", "#F0E68C", "#FFFACD", "#40E0D0", "#DA70D6","#B0C4DE", "#66CDAA",'#232356', 'blue','#f16913',"#D2B48C", "#CD853F","#BC8F8F",'#1f78b4', "#E6E6FA", "#DB7093", "#FF82AB", "#C0C0C0", "#C0D9D9", "#FDEE73", "#4682B4", "#B0E0E6", 'black'
          ,"limegreen","hotpink","brown","turquoise", "#012345", "#555555", "#432321",  "#333333", "#B22222", "#556B2F", "#8A2BE2", "#20B2AA", "#800080", "#4682B4", "#A52A2A","#FFA500")

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
common_samples=intersect(merged_df$Sample, rf_data$Sample)
rf_data=rf_data[rf_data$Sample %in% common_samples,]
merged_df=merged_df[merged_df$Sample %in% common_samples,]
rf_data=rf_data[,c("Genus", "Sample", "Predicted"), drop=FALSE]
merged_df=merge(rf_data, merged_df, by=c("Sample", "Genus"), all.x=TRUE)
merged_df[is.na(merged_df)] = 0
merged_df$Station=substr(merged_df$Sample,1,6)
merged_df$Salinity=all_merged_metabar_df$Salinity[match(merged_df$Sample, all_merged_metabar_df$sample)]
merged_df$Month=all_merged_metabar_df$month[match(merged_df$Sample, all_merged_metabar_df$sample)]
merged_df$Depth=station_mean_depth$Mean_depth[match(merged_df$Station, station_mean_depth$station_id)]
merged_df$Basin=genus_meta_combined$sea_basin[match(merged_df$Sample, genus_meta_combined$Sample)]

names(merged_df)[names(merged_df) == "Microscopy_count"] = "Microscopy"
names(merged_df)[names(merged_df) == "Reads"] = "Metabarcoding"
#cols=c("Microscopy", "Predicted")
#merged_df=merged_df[-5]
#merged_df=merged_df[rowSums(merged_df[cols] == 0) == 0 , ]
merged_df=na.omit(merged_df)

#Use for comparison between stations/months
stations_df=merged_df
months_df=merged_df

#Use for other plots
merged_df = merged_df %>%
  group_by(Sample) %>%
  mutate(Relative_abundance_microscopy = Microscopy / sum(Microscopy),
         Relative_abundance_reads = Metabarcoding / sum(Metabarcoding), 
         Relative_abundance_rf = Predicted / sum(Predicted))


#######################BAR PLOT - ONE PER STATION ##################################################

#Station IDs and names

unique_stations=unique(merged_df$Station)
print(unique_stations)

merged_df = merged_df %>%
  group_by(Sample) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding), 
    Relative_abundance_rf = Predicted / sum(Predicted)
  ) %>%
  filter(
    Relative_abundance_microscopy>0,
    Relative_abundance_reads >0, 
    Relative_abundance_rf >0
  ) 


merged_df = merged_df %>%
  group_by(Sample) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding), 
    Relative_abundance_rf = Predicted / sum(Predicted)
  )



station_ids = c("264194","263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")

station_names = c("Å17", "BY31", "BY2", "B1", "REF M1V1", "BY5", "ANHOLT E", "N14", "SLÄGGÖ", "GA1","C3", "B7", "RA2", "BY15")

#Loop over each station
for (i in seq_along(station_ids)) {
  the_station = station_ids[i]
  station_name = station_names[i]
  
  #Filter data for the station
  grep = grep(the_station, merged_df$Sample)
  filtered_df = merged_df[grep,, drop = FALSE]
  
  #Remove rows with 0 abundance for both methods
  filtered_df = filtered_df %>%
    filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0 & Relative_abundance_rf ==0))
  
  #Pivot data
  df_long = pivot_longer(filtered_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy, Relative_abundance_rf), 
                          names_to = "Method", values_to = "Relative_Abundance")
  
  #Add counts column
  df_long = df_long %>%
    mutate(counts = case_when(
      Method == "Relative_abundance_reads" ~ Metabarcoding,
      Method == "Relative_abundance_microscopy" ~ Microscopy,
      Method == "Relative_abundance_rf" ~ Predicted
    ))
  
  #Rename Method
  df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                           ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding",
                                  ifelse(df_long$Method == "Relative_abundance_rf", "Predicted", df_long$Method)))
  
  #Label
  df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")
  
  #Calculate new relative abundances for each sample
  df_long = df_long %>%
    group_by(Month, Method, Sample, Genus, Station, Salinity) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance)) #%>%
    #filter(Relative_Abundance > 0.05)

  
  #Run when including everything
  #Label genera with low abundance as "Other genera"
  #df_long = df_long %>%
    #mutate(Genus = ifelse(Relative_Abundance > 0.05, as.character(Genus), "Other genera"))

  
  #Plotting
  p = ggplot(df_long, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Method, scales = "free_y") +
    labs(title = paste("Comparison", station_name, sep = " - "),
         x = "", y = "Relative abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 12)) +
    scale_fill_manual(values = myCol)  #Use scale_fill_manual for fill colors
  
  #Print the plot
  print(p)
  
  #Save the plot
  filename = paste("C:/Users/johan/OneDrive/R/Master project/plots/comparison_", the_station, "_unfiltered_finale.png", sep = "")
  ggsave(filename, p, width = 15, height = 6, bg = 'white')


}


#######LOOP OVER EACH BASIN ##############

basins = unique(merged_df$Basin) 

#Loop over each station
for (i in seq_along(basins)) {
  the_basin = basins[i]

  
  #Filter data for the basin
  grep = grep(the_basin, merged_df$Sample)
  filtered_df = merged_df[grep,, drop = FALSE]
  
  #Remove rows with 0 abundance for both methods
  #filtered_df = filtered_df %>%
    #filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0 & Relative_abundance_rf ==0))
  
  #Pivot data
  df_long = pivot_longer(filtered_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy, Relative_abundance_rf), 
                         names_to = "Method", values_to = "Relative_Abundance")
  
  #Add counts column
  df_long = df_long %>%
    mutate(counts = case_when(
      Method == "Relative_abundance_reads" ~ Metabarcoding,
      Method == "Relative_abundance_microscopy" ~ Microscopy,
      Method == "Relative_abundance_rf" ~ Predicted
    ))
  
  #Rename Method
  df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                          ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding",
                                 ifelse(df_long$Method == "Relative_abundance_rf", "Predicted", df_long$Method)))
  
  #Label
  df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")
  
  #Calculate new relative abundances for each sample
  df_long = df_long %>%
    group_by(Method, Sample, Genus, Basin) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance)) #%>%
  #filter(Relative_Abundance > 0.05)
  
  
  #Run when including everything
  #Label genera with low abundance as "Other genera"
  #df_long = df_long %>%
  #mutate(Genus = ifelse(Relative_Abundance > 0.05, as.character(Genus), "Other genera"))
  
  
  #Plotting
  p = ggplot(df_long, aes(x = Sample, y = Relative_Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Method, scales = "free_y") +
    labs(title = paste("Comparison", the_basin, sep = " - "),
         x = "", y = "Relative abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color = "black"),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 8)) +
    scale_fill_manual(values = myCol)  #Use scale_fill_manual for fill colors
  
  #Print the plot
  print(p)
  
  #Save the plot
  filename = paste("C:/Users/johan/OneDrive/R/Master project/plots/comparison_", the_basin, "_unfiltered.png", sep = "")
  ggsave(filename, p, width = 15, height = 6, bg = 'white')
  
  
}

####ONLY  PREDICTED AND MICROSCOPY ##########
unique_stations = unique(merged_df$Station)
print(unique_stations)

merged_df = merged_df %>%
  group_by(Sample) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_rf = Predicted / sum(Predicted)
  )

station_ids = c("264194", "263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")
station_names = c("Å17 (264194)", "BY31 (263732)", "BY2 (263725)", "B1 (263738)", "REF M1V1 (263729)", "BY5 (263728)", "ANHOLT E (264876)", "N14 (263903)", "SLÄGGÖ (263628)", "GA1 (263757)", "C3 (190745)", "B7 (263629)", "RA2 (263856)", "BY15 (263730)")

#Loop over each station
for (i in seq_along(station_ids)) {
  the_station = station_ids[i]
  station_name = station_names[i]
  
  #Filter data for the station
  grep = grep(the_station, merged_df$Sample)
  filtered_df = merged_df[grep,, drop = FALSE]
  
  #Remove rows with 0 abundance for microscopy method
  filtered_df = filtered_df %>%
    filter(!(Relative_abundance_microscopy == 0))
  
  #Pivot data
  df_long = pivot_longer(filtered_df, cols = c(Relative_abundance_microscopy, Relative_abundance_rf), 
                          names_to = "Method", values_to = "Relative_Abundance")
  
  #Add counts column
  df_long = df_long %>%
    mutate(counts = case_when(
      Method == "Relative_abundance_microscopy" ~ Microscopy,
      Method == "Relative_abundance_rf" ~ Predicted
    ))
  
  #Rename Method
  df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                           ifelse(df_long$Method == "Relative_abundance_rf", "Predicted", df_long$Method))
  
  #Label
  df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")
  
  #Calculate new relative abundances for each sample
  df_long = df_long %>%
    group_by(Month, Method, Sample, Genus, Station, Salinity) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance))
  
 
  tryCatch({
    #Plotting
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
      scale_fill_manual(values = myCol)  #Use scale_fill_manual for fill colors
    
    #Print the plot
    print(p)
    
    #Save the plot
    filename = paste("C:/Users/johan/OneDrive/R/Master project/plots/comparison_", the_station, "_rf.png", sep = "")
    ggsave(filename, p, width = 15, height = 6, bg = 'white')
  }, error = function(e) {
    cat("Error occurred for station:", the_station, "\n", conditionMessage(e), "\n")
  })
}
#################SEPARATE SCATTER PLOTS RELATIVE ABUNDACE SAMPLES ################################################


genera_to_keep = c("Oithona", "Calanus", "Pseudocalanus", "Acartia", "Eurytemora", "Paracalanus", "Temora", "Centropages")

#Filter the dataframe to keep only specified genera
final_df = merged_df %>%
  filter(Genus %in% genera_to_keep)

#Calculate correlation and p-values for each genus
cor_values = final_df %>%
  group_by(Genus) %>%
  summarize(
    Correlation = cor(Relative_abundance_microscopy, Relative_abundance_reads, use = "complete.obs"),
    P_Value = cor.test(Relative_abundance_microscopy, Relative_abundance_reads)$p.value
  )

#Remove rows with NA values (if any)
cor_values = na.omit(cor_values)

#Merge the correlation and p-values with the final dataframe
final_df = merge(final_df, cor_values, by = "Genus")

#Extract unique correlations for adding to the plot
unique_correlations = unique(final_df[, c("Genus", "Correlation", "P_Value")])

#Plotting for each genus
for (genus in genera_to_keep) {
  #Subset data for the current genus
  genus_data = subset(final_df, Genus == genus)
  
  #Create the plot
  q = ggplot(genus_data, aes(x = Relative_abundance_microscopy, y = Relative_abundance_reads)) +
    geom_point(color = "black", size = 2, shape = 16) +
    ggtitle(paste(genus)) +
    xlab("Relative abundance microscopy") +
    ylab("Relative abundance metabarcoding") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 50, face = "bold")) +
    scale_x_continuous(limits = c(0, max(genus_data$Relative_abundance_microscopy, na.rm = TRUE))) + 
    scale_y_continuous(limits = c(0, max(genus_data$Relative_abundance_reads, na.rm = TRUE))) + 
    geom_smooth(method = "lm", se = TRUE) + 
    geom_ribbon(aes(ymin = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,2],
                    ymax = predict(lm(Relative_abundance_reads ~ Relative_abundance_microscopy, data = genus_data), interval = "confidence")[,3]),
                fill = "gray", alpha = 0.2)   #Add shaded area around the trend line
  
  #Add Pearson correlation value and p-value to the plot
  corr_value = unique_correlations$Correlation[unique_correlations$Genus == genus]
  p_value = unique_correlations$P_Value[unique_correlations$Genus == genus]
  
  q = q +
    geom_text(aes(label = paste0("R = ", round(corr_value, 2), ", p = ", format.pval(p_value, digits = 2))), 
              x = Inf, y = Inf, color = "blue", hjust = 1, vjust = 1, size = 10)
  
  #Save the plot with a unique name based on genus
  ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/scatter_", gsub(" ", "_", genus), "_avg_samples_finale.png"), q, width = 10, height = 10, bg = 'white')
}


######################SCATTER PLOT - MEAN RELATIVE ABUNDANCE (ALL STATIONS) #######

#Filter the merged dataframe
scatter_df = merged_df %>%
  filter(Genus %in% rownames(count_microscopy))

#Calculate total abundances
total_abundances = scatter_df %>%
  group_by(Genus, Station) %>%
  summarise(Total_microscopy = sum(Relative_abundance_microscopy) / 14,
            Total_metabarcoding = sum(Relative_abundance_reads) / 14)

#Calculate correlations
overall_correlation = total_abundances %>%
  group_by(Genus) %>%
  summarise(Correlation = cor(Total_microscopy, Total_metabarcoding, use = "pairwise.complete.obs"))

#Remove rows with NA values
overall_correlation = na.omit(overall_correlation)

#Create a plot
q = ggplot(total_abundances, aes(x = Total_microscopy, y = Total_metabarcoding)) +
  geom_point(color = "black", size = 2, shape = 16) +
  facet_wrap(~ Genus, scales = "free") +
  ggtitle("Genera") +
  xlab("Mean Relative abundance microscopy") +
  ylab("Mean Relative abundance metabarcoding") +
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
              fill = "gray", alpha = 0.2)   #Add shaded area around the trend line

#Add Pearson values to plots
q = q + 
  geom_text(data = overall_correlation, 
            aes(label = paste0("R = ", round(Correlation, 2)), 
                x = Inf, y = Inf, color = "blue"), 
            hjust = 1, vjust = 1, size = 15)

#Print the plot
print(q)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/scatter_genera_finale.png", q, width = 35, height = 25, bg = 'white')


##################################FAMILY - GENUS - PEARSON table to report ###############################

#Trim whitespace from genus column in family_genus_df
family_genus_df$Genus = trimws(family_genus_df$Genus)

#Trim whitespace from genus column in taxa DataFrame
taxa$Genus = trimws(taxa$Genus)

#Perform left join to retain all genus names from family_genus_df
merge_family_genus = left_join(family_genus_df, taxa[, c("Genus", "Family")], by = "Genus")

#Keep only unique genus in merge_family_genus
merge_family_genus = merge_family_genus[!duplicated(merge_family_genus$Genus), ]
merge_family_genus = merge_family_genus[, c("Family", setdiff(names(merge_family_genus), "Family"))]

library(xtable)

italic_function = function(x) {
  paste("\\textit{", x, "}", sep = "")
}


latex_table = xtable(merge_family_genus)


print(latex_table, sanitize.text.function = italic_function, include.rownames = FALSE)


##################################SCATTER PLOT - SAVE SEPARATELY - MEAN RELATIVE ABUNDANCE STATION  ###############

genera_to_keep = c("Oithona", "Calanus", "Pseudocalanus", "Acartia", "Eurytemora", "Paracalanus", "Temora", "Centropages", "Penilia", "Oikopleura")

#Filter the dataframe to keep only specified genera
scatter_df = merged_df %>%
  filter(Genus %in% rownames(count_microscopy))

#Calculate total abundances
total_abundances = scatter_df %>%
  group_by(Genus, Station) %>%
  summarise(Total_microscopy = sum(Relative_abundance_microscopy) / 14,
            Total_metabarcoding = sum(Relative_abundance_reads) / 14)

#Calculate correlations and p-values
overall_correlation = total_abundances %>%
  group_by(Genus) %>%
  summarize(
    Correlation = cor(Total_microscopy, Total_metabarcoding, use = "pairwise.complete.obs"),
    P_Value = cor.test(Total_microscopy, Total_metabarcoding)$p.value
  )

#Remove rows with NA values
overall_correlation = na.omit(overall_correlation)

#Loop through each genus
for (genus in genera_to_keep) {
  #Filter data for the current genus
  genus_df = total_abundances %>%
    filter(Genus == genus)
  
  #Create the plot
  q = ggplot(genus_df, aes(x = Total_microscopy, y = Total_metabarcoding)) +
    geom_point(color = "black", size = 5, shape = 16) +
    ggtitle(paste(genus)) +
    xlab("Mean relative abundance microscopy") +
    ylab("Mean relative abundance metabarcoding") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 50, face = "bold")) +
    scale_x_continuous(limits = c(0, max(genus_df$Total_microscopy))) + 
    scale_y_continuous(limits = c(0, max(genus_df$Total_metabarcoding)))+ 
    geom_smooth(method = "lm", se = TRUE) + 
    geom_ribbon(aes(ymin = predict(lm(Total_metabarcoding ~ Total_microscopy, data = genus_df), interval = "confidence")[,2],
                    ymax = predict(lm(Total_metabarcoding ~ Total_microscopy, data = genus_df), interval = "confidence")[,3]),
                fill = "gray", alpha = 0.2)   #Add shaded area around the trend line
  
  #Add Pearson correlation value and p-value to plot
  corr_value = overall_correlation$Correlation[which(overall_correlation$Genus == genus)]
  p_value = overall_correlation$P_Value[which(overall_correlation$Genus == genus)]
  
  q = q +
    geom_text(aes(label = paste0("R = ", round(corr_value, 2), ", p = ", format.pval(p_value, digits = 2))), 
              x = Inf, y = Inf, color = "blue", hjust = 1, vjust = 1, size = 10)
  
  #Save the plot with a unique name based on genus
  ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/scatter_", gsub(" ", "_", genus), "_avgstation_finale.png"), q, width = 10, height = 10, bg = 'white')
}

genera_list = unique(scatter_df$Genus)

#Loop through each genus, create a plot, and save it
for (genus in genera_list) {
  #Subset data for the current genus
  genus_data = subset(scatter_df, Genus == genus)
  
  #Plot
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
                fill = "gray", alpha = 0.2)   #Add shaded area around the trend line
  
  #Save the plot with a unique name based on genus
  ggsave(paste0("C:/Users/johan/OneDrive/R/Master project/plots/scatter_", gsub(" ", "_", genus), "_avg.png"), q, width = 10, height = 10, bg = 'white')
}
######################BAR PLOT (ONE BAR PER STATION, ALL YEARS) ##########################

#Only keep one year
#merged_df = merged_df %>%
  #filter(str_detect(Sample, "2020"))


station_ids = c("264194", "263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")
station_names = c("Å17", "BY31", "BY2", "B1", "REF M1V1", "BY5", "ANHOLT E", "N14", "SLÄGGÖ ", "GA1", "C3", "B7", "RA2", "BY15")

#Replace station IDs with station names in stations_df
stations_df = stations_df %>%
  mutate(Station = ifelse(Station %in% station_ids, station_names[match(Station, station_ids)], Station))

stations_df = stations_df %>%
  group_by(Station) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding),
    Relative_abundance_rf = Predicted / sum(Predicted)
  ) %>%
  filter(
    Relative_abundance_microscopy>0,
    Relative_abundance_reads >0, 
    Relative_abundance_rf >0
  ) 


stations_df = stations_df %>%
  group_by(Station) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding),
    Relative_abundance_rf = Predicted / sum(Predicted))


#the_station='263856'
#grep=grep(the_station, my_df$Sample)
#my_df=my_df[grep,, drop=FALSE]
#stations_df=stations_df %>%
  #filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))

df_long = pivot_longer(stations_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy, Relative_abundance_rf), 
                        names_to = "Method", values_to = "Relative_Abundance")



#Rename Method
df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                        ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding",
                               ifelse(df_long$Method == "Relative_abundance_rf", "Predicted", df_long$Method)))
                        
df_long  = df_long %>%
  filter(!is.nan(Relative_Abundance))


df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")


df_long = subset(df_long, select = -c(Microscopy, Metabarcoding, Predicted) )

df_long = df_long %>%
  group_by(Station)%>%
  mutate(Relative_Abundance=Relative_Abundance/sum(Relative_Abundance))

#df_long = df_long %>%
  #mutate(Genus = ifelse(Relative_Abundance > 0.01, as.character(Genus), "Other genera"))

p = ggplot(df_long, aes(x = reorder(Station, Salinity), y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Spatial variation of microbial composition across stations",
       x = "Station", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = myCol)  #Use scale_fill_manual for fill colors

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_stations_filtered.png", p, width = 15, height = 7, bg = 'white')




###################VENN DIAGRAM MICROSCOPY AND METABARCODING ALL SAMPLES ######

genus_list_meta=unique(all_merged_metabar_df$Genus)
genus_list_micro=rownames(count_microscopy)
genus_list_meta=na.omit(genus_list_meta)
genus_list_micro=na.omit(genus_list_micro)
genus_list_common=intersect(genus_list_meta, genus_list_micro)

library(VennDiagram)

windows()
venn.plot = venn.diagram(
  x = list(Metabarcoding = genus_list_meta, Microscopy = genus_list_micro), 
  category.names = c("Metabarcoding", "Microscopy"),
  filename = NULL,
  main = "Number of genera identified",
  col = "transparent", #Set circle colors to transparent to change them individually
  fill = c("#FF9800", "#87CEEB"), #Set fill colors for circles
  cat.pos = c(0, 0), #Set category positions to place names outside the circles
  cat.dist = c(0.1, 0.1), #Set distance of category names from the circles
  cat.fontfamily = "sans", #Set font family for category names
  cat.fontface = "bold", #Set font face for category names
  cat.col = c("#FF9800", "#87CEEB"), #Set colors for category names
  cat.cex = 1, #Set size of category names
  title.fontfamily = "sans", #Set font family for title
  title.fontface = "bold", #Set font face for title
  title.col = "black", #Set color for title
  title.cex = 1.5 #Set size of title
)

grid.draw(venn.plot)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/venn_diagram_240517.png", venn.plot, width = 5, height = 5, bg = 'white')

##############TOTAL MICROSCOPY COUNTS IN COMMON SAMPLES #######################

count_microscopy[] = lapply(count_microscopy, as.numeric)
sums_count = rowSums(count_microscopy, na.rm = TRUE)
total_count = data.frame(row.names = rownames(count_microscopy), Sum = sums_count)
total_count = total_count[rowSums(total_count) != 0, , drop = FALSE]
total_count = total_count %>%
  mutate(Relative_abundance = Sum / sum(Sum))

total_count = total_count %>%
  filter(Relative_abundance >= 0.01)
windows()

p = ggplot(total_count, aes(x = reorder(row.names(total_count), -Relative_abundance), y = Relative_abundance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Most observed genera in microscopy",
       x = "Genera",
       y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 30),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 30))

#Print the plot
print(p)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/bar_counts_rel_ab.png", p, width = 20, height = 6, bg = 'white')


#############TOTAL METABARCODING IN COMMON SAMPLES #########################

sum_reads = rowSums(new_seqtab_genus, na.rm=TRUE)
total_reads = data.frame(row.names = rownames(new_seqtab_genus), Sum=sum_reads )
total_reads = total_reads %>%
  mutate(Relative_abundance = Sum / sum(Sum))

total_reads = total_reads %>%
  filter(Relative_abundance >= 0.05)


windows()
p=ggplot(total_reads, aes(x = reorder(row.names(total_reads), -Relative_abundance), y = Relative_abundance)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  labs(title = "Most observed genera in metabarcoding",
       x = "Genera",
       y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, color = "black"),
        legend.position = "right",
        plot.title = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 30))

print(p)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/bar_meta_rel_ab.png", p, width = 10, height = 6, bg = 'white')



###TEST

library(ggplot2)

generate_plot = function(data, title, filename) {
  p = ggplot(data, aes(x = reorder(row.names(data), -Relative_abundance), y = Relative_abundance)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    ggtitle(title) +  #Set the plot title directly
    xlab("Genera") +  #Set the x-axis label
    ylab("Relative Abundance") +  #Set the y-axis label
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15, color = "black"),
      legend.position = "right",
      plot.title = element_text(size = 25),  #Adjust title size here
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 30)
    )
  
  print(p)
  ggsave(filename, p, width = 12, height = 6, bg = 'white')
}

#Call the function to generate the first plot (microscopy)
generate_plot(total_count, "Genera in microscopy", "C:/Users/johan/OneDrive/R/Master project/plots/bar_microscopy.png")

#Call the function to generate the second plot (metabarcoding)
generate_plot(total_reads, "Genera in metabarcoding", "C:/Users/johan/OneDrive/R/Master project/plots/bar_metabarcoding.png")

######################BAR PLOT MONTH ########################################

#Calculate relative abundance columns
months_df = months_df %>%
  group_by(Month) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding),
    Relative_abundance_rf = Predicted / sum(Predicted)
  ) %>%
  filter(
    Relative_abundance_microscopy > 0,
    Relative_abundance_reads > 0,
    Relative_abundance_rf > 0
  ) 

months_df = months_df %>%
  group_by(Month) %>%
  mutate(
    Relative_abundance_microscopy = Microscopy / sum(Microscopy),
    Relative_abundance_reads = Metabarcoding / sum(Metabarcoding),
    Relative_abundance_rf = Predicted / sum(Predicted)
  )

#Filter out rows with zero abundance for both Microscopy and Metabarcoding
#months_df = months_df %>%
  #filter(!(Relative_abundance_reads == 0 & Relative_abundance_microscopy == 0))

#Omit rows with NA values
months_df = na.omit(months_df)

#Reshape data to long format
df_long = pivot_longer(months_df, cols = c(Relative_abundance_reads, Relative_abundance_microscopy, Relative_abundance_rf), 
                        names_to = "Method", values_to = "Relative_Abundance")

#Create a new column for counts based on the selected Method
df_long = df_long %>%
  mutate(
    counts = case_when(
      Method == "Relative_abundance_reads" ~ Metabarcoding,
      Method == "Relative_abundance_microscopy" ~ Microscopy,
      Method == "Relative_abundance_rf" ~ Predicted
    )
  )

#Rename Method column values
df_long$Method = ifelse(df_long$Method == "Relative_abundance_microscopy", "Microscopy",
                         ifelse(df_long$Method == "Relative_abundance_reads", "Metabarcoding", "Predicted"))

#Create Label column
df_long$Label = paste(df_long$Sample, df_long$Method, sep = " - ")

#Subset unnecessary columns
df_long = subset(df_long, select = -c(Microscopy, Metabarcoding, Predicted))

#Factorize Month column
df_long$Month = factor(df_long$Month, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))

#Modify Genus column based on abundance threshold
#df_long = df_long %>%
  #mutate(Genus = ifelse(Relative_Abundance > 0.005, as.character(Genus), "Other genera"))

#Plot
windows()
p = ggplot(df_long, aes(x = Month, y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Temporal variation of microbial composition across samples",
       x = "Month", y = "Relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = myCol)

#Print plot
print(p)

#Save plot
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/comparison_months_filtered.png", p, width = 15, height = 7, bg = 'white')


####CALCULATE HOW MANY SAMPLES EACH GENUS IS FOUND IN ##########

long_df=merged_df %>%
  pivot_longer(cols = c(Microscopy, Metabarcoding, Predicted),
               names_to = "Method",
               values_to = "Presence") %>%
  filter(Presence > 0)

#Count the number of samples for each genus and method
sample_counts=long_df %>%
  group_by(Genus, Method) %>%
  summarise(Count = n_distinct(Sample))


