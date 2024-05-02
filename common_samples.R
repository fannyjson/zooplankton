
#Set working directory

setwd("C:/Users/johan/OneDrive/R/Master project")

#Load files and libraries

#Metazoa only
taxa_metazoa=read.delim("C:/Users/johan/OneDrive/R/Master project/taxa_metazoa_processed.tsv")

#Normalized seqtab metazoa
seqtab_processed=read.delim("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_18S_processed.tsv")

#Metadata
metadata_processed=read.delim("C:/Users/johan/OneDrive/R/Master project/metadata_processed.tsv")

#SMHI data
smhi_data=read.delim("C:/Users/johan/OneDrive/R/Master project/zooplankton_2015_2020_2024-01-25_utf8.txt")

#Count table 
count_microscopy=read.delim("C:/Users/johan/OneDrive/R/Master project/count_table_zooplankton_genus.tsv")

colnames(count_microscopy) = gsub(colnames(count_microscopy), pattern = '^X', replacement = '')
colnames(count_microscopy) = gsub(colnames(count_microscopy), pattern = '\\.', replacement = '-')

colnames(seqtab_processed) = gsub(colnames(seqtab_processed), pattern = '^X', replacement = '')
colnames(seqtab_processed) = gsub(colnames(seqtab_processed), pattern = '\\.', replacement = '-')

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(vegan)

#Make dataframes with unique taxa for every taxlevel

taxa_class=taxa_metazoa[,5, drop=FALSE]
taxa_class=na.omit(taxa_class)
taxa_class_unique=unique(taxa_class)

taxa_order=taxa_metazoa[,6, drop=FALSE]
taxa_order=na.omit(taxa_order)
taxa_order_unique=unique(taxa_order)

taxa_family=taxa_metazoa[,7, drop=FALSE]
taxa_family=na.omit(taxa_family)
taxa_family_unique=unique(taxa_family)
                         
taxa_genus= taxa_metazoa[,8, drop=FALSE]
taxa_genus= na.omit(taxa_genus)
taxa_genus_unique=unique(taxa_genus)

taxa_species=taxa_metazoa[,9, drop=FALSE]
taxa_species=na.omit(taxa_species)
taxa_species_unique=unique(taxa_species)

#Loop through each column of seqtab_processed and store cells>0. 
#The corresponding ASVs and taxa for each sample are stored in taxa_metabar

extract_taxa_list = function(seqtab_processed, taxa_metazoa, taxa_level) {
  taxa_list = list()
  
  for (i in 1:ncol(seqtab_processed)) {
    present_indices = which(seqtab_processed[, i] > 0)
    
    if (length(present_indices) > 0) {
      present_asvs = rownames(seqtab_processed)[present_indices]
      
      #Check if present ASVs are at the specified taxonomic level
      
      present_asvs_taxa = intersect(present_asvs, rownames(taxa_metazoa))
      
      if (length(present_asvs_taxa) > 0) {
        present_taxa = taxa_metazoa[rownames(taxa_metazoa) %in% present_asvs_taxa, ]
        
        taxa_list[[i]] = data.frame(
          sample_id = colnames(seqtab_processed)[i],
          ASV = present_asvs_taxa,
          Taxa = present_taxa
        )
      }
    }
  }
  taxa_metabar=do.call(rbind, taxa_list)
}

taxa_metabar=extract_taxa_list(seqtab_processed, taxa_metazoa, 5)
rownames(taxa_metabar)=NULL
colnames(taxa_metabar)=gsub("^Taxa\\.", "", colnames(taxa_metabar))

#Combines taxa found in the metabar data with the metadata and only keep unique rows

combined_df = merge(metadata_processed, taxa_metabar, by = "sample_id")
names(combined_df)[names(combined_df)=="station_name"] = "Station"
names(combined_df)[names(combined_df)=="date"] = "Date"

combined_df$Date=ymd(combined_df$Date)
combined_df$Station = gsub("\\s{2,}", " ", combined_df$Station)

#write.table(combined_df, "C:/Users/johan/OneDrive/R/Master project/merged_df_all_genus.tsv", sep="\t")


################# STANDARDIZE STATION NAMES AND SAMPLES#########################


#Coherent column names for both datasets

year_smhi=2015:2020
smhi_data=smhi_data[year(combined_df$Date)==year_smhi,]
names(smhi_data)[names(smhi_data) == "sample_date"] = "Date"
names(smhi_data)[names(smhi_data) == "reported_station_name"] = "Station"


#Standardize stations in smhi data

smhi_data$Date=ymd(smhi_data$Date)
smhi_data$Station = gsub("\\s{2,}", " ", smhi_data$Station)

standardize_stations = data.frame(Original = c("BY5", "BY2", "Släggö", "BY15", "Ref M1V1"),
                                   Standardized = c("BY5 BORNHOLMSDJ", "BY2 ARKONA", "SLÄGGÖ", "BY15 GOTLANDSDJ", "REF M1V1"),
                                   stringsAsFactors = FALSE)

replace_station_names = function(station_name) {
  mapping_row = standardize_stations[standardize_stations$Original == station_name, ]
  if (nrow(mapping_row) > 0) {
    return(mapping_row$Standardized)
  } else {
    return(station_name)
  }
}

smhi_data$Station = sapply(smhi_data$Station, replace_station_names)


#Match stations in smhi_data with stations in metabar data

#Standardize stations in combined_df

standardize_stations = data.frame(Original = c("GAVIK-1", "RÅNEÅ-2"),
                                   Standardized = c("GA1", "RA2"),
                                   stringsAsFactors = FALSE)

combined_df$Station = sapply(combined_df$Station, replace_station_names)

# Create common sample names (same as microscopy count table sample names)
station_ids=smhi_data[,c("Station", "station_id")]
station_ids=unique(station_ids)
metabar_df=merge(combined_df, station_ids, by = "Station") #all.x=TRUE will give all stations
metabar_df$sample=paste0(metabar_df$station_id, "_", metabar_df$Date)

#Translation dataframe

translate_sample_name=metabar_df[,c("sample_id", "sample")]

#Rename colnames for seqtab to new sample names using translation table

for (old_name in colnames(seqtab_processed)) {
  index = match(old_name, translate_sample_name$sample_id)
  if (!is.na(index)) {
    new_name = translate_sample_name$sample[index]
    colnames(seqtab_processed)[colnames(seqtab_processed) == old_name] = new_name
  }
}

#Extract samples present in both microscopy and metabarcoding

smhi_data$sample=paste0(smhi_data$station_id, "_", smhi_data$Date)
smhi_data$sample = as.character(smhi_data$sample)
metabar_df$sample = as.character(metabar_df$sample)
common_samples=intersect(smhi_data$sample, metabar_df$sample)
metabar_df = metabar_df[metabar_df$sample %in% common_samples, ]
smhi_data = smhi_data[smhi_data$sample %in% common_samples, ]

#Check what stations have common station id:s

common_stations = intersect(smhi_data$station_id, metabar_df$station_id)
print(unique(common_stations))

#Save smhi_data and metabar_df

write.table(smhi_data, "C:/Users/johan/OneDrive/R/Master project/smhi_data_common_samples_240411.tsv", sep="\t")
write.table(metabar_df, "C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_PR2_240411.tsv", sep="\t")

#Only keep rows that have annotations on genus level

genus_col=taxa_metazoa$Genus
genus_df=data.frame(Genus = genus_col, row.names = rownames(taxa_metazoa))
matching_asvs=intersect(rownames(seqtab_processed), rownames(genus_df))
genus_df=genus_df[matching_asvs, , drop = FALSE]
genus_df=na.omit(genus_df)
matching_asvs=intersect(rownames(seqtab_processed), rownames(genus_df))

#Subset seqtab_processed based on matching ASVs

seqtab_processed_subset=seqtab_processed[matching_asvs, ]

# Combine data frames

new_seqtab=cbind(genus_df, seqtab_processed_subset)

# Make seqtab with genus instead of asvs and only for common samples
#common_columns=intersect(colnames(new_seqtab), colnames(count_microscopy))

keep_col=c("Genus",common_samples)
new_seqtab=new_seqtab[,keep_col, drop=FALSE]
new_seqtab_genus = aggregate(. ~ Genus, data = new_seqtab, FUN = sum)
rownames(new_seqtab_genus) = new_seqtab_genus$Genus
new_seqtab_genus=new_seqtab_genus[,-1]
count_microscopy=count_microscopy[,common_samples,drop=FALSE]
count_microscopy = count_microscopy %>%
  filter(rowSums(. != 0, na.rm = TRUE) > 0)


#write.table(metabar_df, "C:/Users/johan/OneDrive/R/Master project/merged_metabar_all_genera_240411.tsv", sep="\t")
write.table(count_microscopy, "C:/Users/johan/OneDrive/R/Master project/count_microscopy_commonsamples_240411.tsv", sep="\t")
write.table(new_seqtab_genus, "C:/Users/johan/OneDrive/R/Master project/new_seqtab_genus_240411.tsv", sep="\t")

#Run if only common genera should be included

metabar_df=metabar_df[metabar_df$Genus %in% rownames(count_microscopy), , drop=FALSE] 
write.table(metabar_df, "C:/Users/johan/OneDrive/R/Master project/merged_metabar_common_genera_240411.tsv", sep="\t")
new_seqtab_common_genus=new_seqtab_genus[rownames(new_seqtab_genus) %in% rownames(count_microscopy), , drop=FALSE]
write.table(new_seqtab_common_genus, "C:/Users/johan/OneDrive/R/Master project/new_seqtab_common_genera_240411.tsv", sep="\t")


############## NMDS for genera in common samples ##############################

######## MICROSCOPY NMDS ##########################

matrix = t(count_microscopy)

#Extract station IDs from sample names

sample=colnames(new_seqtab_genus)
season = metabar_df$season[match(sample, metabar_df$sample)]
salinity = metabar_df$Salinity[match(sample,metabar_df$sample)]
station_ids = substr(rownames(matrix), 1, 6)
stations = c("264194","263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")
station_names = c("Å17 (264194)", "BY31 (263732)", "BY2 (263725)", "B1 (263738)", "REF M1V1 (263729)", "BY5 (263728)", "ANHOLT E (264876)", "N14 (263903)", "SLÄGGÖ (263628)", "GA1 (263757)","C3 (190745)", "B7 (263629)", "RA2 (263856)", "BY15 (263730)")


translate_table = setNames(station_names, stations)

#Replace station IDs with names
station_ids = translate_table[station_ids]

#Calculate Bray-Curtis dissimilarity
bray_matrix = vegdist(matrix, method = "bray")

#Perform NMDS
nmds_result = metaMDS(bray_matrix)

#Set colors
colors = c("#232356", "#FFD700", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#4CAF50", "#FF82AB","#673AB7", "#FF4081",'#1f78b4','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')

#Extract coordinates from NMDS result
nmds_coordinates = data.frame(nmds_result$points)

#shannon_values = diversity(matrix, index = "shannon")

# Combine NMDS coordinates with metadata
nmds_data = cbind(nmds_coordinates, Salinity=salinity, Station_ID = station_ids, Season=season)
#nmds_data = cbind(nmds_coordinates, Shannon_Diversity = shannon_values, Station_ID = station_ids, Season=season)

#Define the range of point sizes based on salinity
max_size=10  
min_size=1  
range_salinity=range(nmds_data$Salinity)
#range_shannon=range(nmds_data$Shannon_Diversity)

#Scale salinity to point sizes
scaled_sizes = rescale(nmds_data$Salinity, to = c(min_size, max_size))

windows()

#Plot NMDS with different colors for different stations and different sizes for salinity
final_nmds = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Station_ID, size = scaled_sizes)) +
  geom_point() +
  ggtitle("NMDS microscopy") +
  theme_bw() +  
  scale_color_manual(name = "Station", values = colors) +
  scale_size_continuous(name = "Salinity", range = c(min_size, max_size), guide=FALSE)


print(final_nmds)

#Save NMDS plot
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_microscopy_salinity.png",final_nmds, width=10, height=6)


#Plot seasonal variations

windows()

final_nmds = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Station_ID, shape = Season)) +
  geom_point(size=4) +
  ggtitle("NMDS microscopy") +
  theme_bw() + 
  scale_color_manual(name = "Station", values = colors) +
  scale_shape_manual(name = "Season", 
                     values = c("Spring" = 0, "Summer" = 1, "Fall" = 2, "Winter" = 3))
print(final_nmds)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_microscopy_season.png", final_nmds, width = 10, height = 6)


################ METABARCODING NMDS GENERA ############################################

matrix = t(new_seqtab_genus)

#Extract station IDs from sample names
sample=colnames(new_seqtab_genus)
season = metabar_df$season[match(sample, metabar_df$sample)]
salinity = metabar_df$Salinity[match(sample,metabar_df$sample)]
station_ids = substr(rownames(matrix), 1, 6)
stations = c("264194","263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")
station_names = c("Å17 (264194)", "BY31 (263732)", "BY2 (263725)", "B1 (263738)", "REF M1V1 (263729)", "BY5 (263728)", "ANHOLT E (264876)", "N14 (263903)", "SLÄGGÖ (263628)", "GA1 (263757)","C3 (190745)", "B7 (263629)", "RA2 (263856)", "BY15 (263730)")

#Translate table for stations
translate_table = setNames(station_names, stations)

#Replace station IDs with names

station_ids = translate_table[station_ids]

#Calculate Bray-Curtis dissimilarity

bray_matrix = vegdist(matrix, method = "bray")

#Perform NMDS

nmds_result = metaMDS(bray_matrix)

colors = c("#232356", "#FFD700", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#4CAF50", "#FF82AB","#673AB7", "#FF4081",'#1f78b4','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')

# Extract coordinates and Shannon diversity from the NMDS result

nmds_coordinates = data.frame(nmds_result$points)

#shannon_values = diversity(matrix, index = "shannon")

#Combine NMDS coordinates and metadata

nmds_data = cbind(nmds_coordinates, Salinity=salinity, Station_ID = station_ids, Season=season)
#nmds_data = cbind(nmds_coordinates, Shannon_Diversity = shannon_values, Station_ID = station_ids, Season=season)

#Define the range of point sizes based on Shannon diversity

max_size=10 
min_size=1   
range_salinity=range(nmds_data$Salinity)

#Scale salinity to point sizes

scaled_sizes=rescale(nmds_data$Salinity, to = c(min_size, max_size))

windows()

#Plot NMDS with different colors for different stations and different sizes for salinity

final_nmds = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Station_ID, size = scaled_sizes)) +
  geom_point() +
  ggtitle("NMDS metabarcoding") +
  theme_bw() +  
  scale_color_manual(name = "Station", values = colors) +
  scale_size_continuous(name = "Salinity", range = c(min_size, max_size), guide=FALSE)


print(final_nmds)
ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_metabarcoding_salinity.png",final_nmds, width=10, height=6)

#Seasonal NMDS metabarcoding

windows()

final_nmds = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Station_ID, shape = Season)) +
  geom_point(size=5) +
  ggtitle("NMDS metabarcoding") +
  theme_bw() + 
  scale_color_manual(name = "Station", values = colors) +
  scale_shape_manual(name = "Season", 
                     values = c("Spring" = 0, "Summer" = 1, "Fall" = 2, "Winter" = 3)) 
print(final_nmds)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_metabarcoding_seasonal.png", final_nmds, width = 10, height = 6)



################ METABARCODING NMDS ASVs ############################################

seqtab_asvs=seqtab_processed[, colnames(seqtab_processed) %in% common_samples , drop=FALSE]

matrix = t(seqtab_asvs)

sample=colnames(new_seqtab_genus)
season = metabar_df$season[match(sample, metabar_df$sample)]
station_ids = substr(rownames(matrix), 1, 6)
salinity = metabar_df$Salinity[match(sample,metabar_df$sample)]
stations = c("264194","263732", "263725", "263738", "263729", "263728", "264876", "263903", "263628", "263757", "190745", "263629", "263856", "263730")
station_names = c("Å17 (264194)", "BY31 (263732)", "BY2 (263725)", "B1 (263738)", "REF M1V1 (263729)", "BY5 (263728)", "ANHOLT E (264876)", "N14 (263903)", "SLÄGGÖ (263628)", "GA1 (263757)","C3 (190745)", "B7 (263629)", "RA2 (263856)", "BY15 (263730)")

#Translate stations
translate_table = setNames(station_names, stations)

#Replace station IDs with names
station_ids = translate_table[station_ids]

#Calculate Bray-Curtis dissimilarity
bray_matrix = vegdist(matrix, method = "bray")

#Perform NMDS
nmds_result = metaMDS(bray_matrix)

colors = c("#232356", "#FFD700", "#87CEEB", "#98FB98", "#F08080", "#DDA0DD", "#4CAF50", "#FF82AB","#673AB7", "#FF4081",'#1f78b4','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')

#Extract coordinates from nmds

nmds_coordinates = data.frame(nmds_result$points)

#shannon_values = diversity(matrix, index = "shannon")

# Combine NMDS coordinates and Shannon diversity values

nmds_data = cbind(nmds_coordinates, Salinity = salinity, Station_ID = station_ids, Season=season)

#Range of point sizes based on salinity

max_size = 10
min_size = 1   
range_salinity = range(nmds_data$Salinity)

#Scale salinity to point sizes

scaled_sizes = rescale(nmds_data$Salinity, to = c(min_size, max_size))

#Plot NMDS with different colors for different stations and salinity

windows()

final_nmds = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Station_ID, size = scaled_sizes)) +
  geom_point() +
  ggtitle("NMDS ASVs PR2") +
  theme_bw() +  # Set theme to black and white (optional)
  scale_color_manual(name = "Station", values = colors) +
  scale_size_continuous(name = "Salinity", range = c(min_size, max_size), guide=FALSE) # Add more shapes as needed

print(final_nmds)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_PR2_salinity_asvs.png",final_nmds, width=10, height=6)

#Seasonal NMDS metabarcoding ASVs

windows()

final_nmds = ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = Station_ID, shape = Season)) +
  geom_point(size=5) +
  ggtitle("NMDS ASVs PR2") +
  theme_bw() + 
  scale_color_manual(name = "Station", values = colors) +
  scale_shape_manual(name = "Season", 
                     values = c("Spring" = 0, "Summer" = 1, "Fall" = 2, "Winter" = 3)) 
print(final_nmds)

ggsave("C:/Users/johan/OneDrive/R/Master project/plots/NMDS_PR2_seasonal_asvs.png", final_nmds, width = 10, height = 6)
