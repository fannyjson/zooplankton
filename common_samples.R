
#Set working directory
setwd("C:/Users/johan/OneDrive/R/Master project")

#Load files and libraries

#Metazoa only
taxa_metazoa=read.delim("C:/Users/johan/OneDrive/R/Master project/taxa_metazoa.tsv")

#Normalized seqtab metazoa
seqtab_processed=read.delim("C:/Users/johan/OneDrive/R/Master project/norm_seqtab_18S_240308.tsv")

#Metadata
metadata_processed=read.delim("C:/Users/johan/OneDrive/R/Master project/metadata_240308.tsv")

#SMHI data
smhi_data=read.table("C:/Users/johan/OneDrive/R/Master project/zooplankton_2015_2020_2024-01-25_utf8.txt", header = TRUE, sep = "\t")

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

############### MERGE METADATA AND GENUS #######################################
#Only keep one taxa level 
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

#Loop through each column of seqtab_processed and store indices>0. 
#The corresponding ASVs for each sample and genera are stored in taxa_metabar. 
#If the ASV is present but not on genus level, it is not added. 

taxa_list_genus=list()

for (i in 1:ncol(seqtab_processed)) {
  present_indices=which(seqtab_processed[,i] > 0)
  

  if (length(present_indices)>0) {
    present_asvs=rownames(seqtab_processed)[present_indices]
    
    if(any(present_asvs %in% rownames(taxa_genus))){
      present_asvs=present_asvs[present_asvs %in% rownames(taxa_genus)]
    
      present_genus=taxa_genus[rownames(taxa_genus) %in% present_asvs,]
      taxa_list_genus[[i]]=data.frame(sample_id=colnames(seqtab_processed)[i],ASV=present_asvs,
      Genus=present_genus)
   }
  }
}

taxa_metabar = do.call(rbind, taxa_list_genus)

#Combines the taxa_metabar with the metadata and only keep unique rows. 
taxa_metabar$sample_id = sub("^X", "", taxa_metabar$sample_id)
combined_df = merge(metadata_processed, taxa_metabar, by = "sample_id")
combined_df = combined_df[, !(names(combined_df) %in% c("ASV"))]
combined_df = unique(combined_df)
names(combined_df)[names(combined_df)=="station_name"] <- "Station"
names(combined_df)[names(combined_df)=="date"] <- "Date"

combined_df$Date=ymd(combined_df$Date)
combined_df$Station <- gsub("\\s{2,}", " ", combined_df$Station)

write.table(combined_df, "C:/Users/johan/OneDrive/R/Master project/merged_df_genus.tsv", sep="\t")

################# STANDARDIZE STATION NAMES AND SAMPLES#########################

genus_meta_combined=combined_df

#Coherent column names for both datasets
year_smhi=2015:2020
smhi_data=smhi_data[year(genus_meta_combined$Date)==year_smhi,]
names(smhi_data)[names(smhi_data) == "sample_date"] <- "Date"
names(smhi_data)[names(smhi_data) == "reported_station_name"] <- "Station"
names(smhi_data)[names(smhi_data) == "scientific_name"] <- "Genus"

#Standardize stations
smhi_data$Date=ymd(smhi_data$Date)
smhi_data$Station <- gsub("\\s{2,}", " ", smhi_data$Station)

standardize_stations <- data.frame(Original = c("BY5", "BY2", "Släggö", "BY15", "Ref M1V1"),
                                   Standardized = c("BY5 BORNHOLMSDJ", "BY2 ARKONA", "SLÄGGÖ", "BY15 GOTLANDSDJ", "REF M1V1"),
                                   stringsAsFactors = FALSE)

replace_station_names <- function(station_name) {
  mapping_row <- standardize_stations[standardize_stations$Original == station_name, ]
  if (nrow(mapping_row) > 0) {
    return(mapping_row$Standardized)
  } else {
    return(station_name)
  }
}

smhi_data$Station <- sapply(smhi_data$Station, replace_station_names)

#Match stations in smhi_data with stations in metabar data

standardize_stations <- data.frame(Original = c("GAVIK-1", "RÅNEÅ-2"),
                                   Standardized = c("GA1", "RA2"),
                                   stringsAsFactors = FALSE)

genus_meta_combined$Station = sapply(genus_meta_combined$Station, replace_station_names)

# Create common sample names (same as microscopy count table sample names)
station_ids=smhi_data[,c("Station", "station_id")]
station_ids=unique(station_ids)
metabar_df=merge(genus_meta_combined, station_ids, by = "Station") #all.x=TRUE will give all stations
metabar_df$sample=paste0(metabar_df$station_id, "_", metabar_df$Date)

#Translation df
translate_sample_name=metabar_df[,c("sample_id", "sample")]

#Rename colnames for seqtab 
for (old_name in colnames(seqtab_processed)) {
  index = match(old_name, translate_sample_name$sample_id)
  if (!is.na(index)) {
    new_name = translate_sample_name$sample[index]
    colnames(seqtab_processed)[colnames(seqtab_processed) == old_name] = new_name
  }
}


#Only keep rows with stations that are common for both datasets
stations_smhi = unique(smhi_data$Station)
stations_combined = unique(genus_meta_combined$Station)
common_stations = intersect(stations_smhi, stations_combined)
smhi_data = smhi_data[smhi_data$Station %in% common_stations, ]
genus_meta_combined = genus_meta_combined[genus_meta_combined$Station %in% common_stations, ]

print(common_stations)

#Extract samples present in both df:s
smhi_data$sample=paste0(smhi_data$station_id, "_", smhi_data$Date)
common_samples=intersect(smhi_data$sample, metabar_df$sample)
metabar_df = metabar_df[metabar_df$sample %in% common_samples, ]
smhi_data = smhi_data[smhi_data$sample %in% common_samples, ]


write.table(smhi_data, "C:/Users/johan/OneDrive/R/Master project/smhi_data_commonsamples", sep="\t")


# genus_df with genus and asv column
genus_col=taxa_metazoa$Genus
genus_df=data.frame(Genus = genus_col, row.names = rownames(taxa_metazoa))
matching_asvs=intersect(rownames(seqtab_processed), rownames(genus_df))
genus_df=genus_df[matching_asvs, , drop = FALSE]
genus_df=na.omit(genus_df)
matching_asvs=intersect(rownames(seqtab_processed), rownames(genus_df))

# Subset seqtab_processed based on matching ASVs
seqtab_processed_subset=seqtab_processed[matching_asvs, ]

# Combine data frames without row names causing issues
new_seqtab=cbind(genus_df, seqtab_processed_subset)

# Make seqtab with genus instead of asvs and only for common samples
common_columns=intersect(colnames(new_seqtab), colnames(count_microscopy))
keep_col=c("Genus",common_columns)
new_seqtab=new_seqtab[,keep_col, drop=FALSE]
new_seqtab_genus = aggregate(. ~ Genus, data = new_seqtab, FUN = sum)
rownames(new_seqtab_genus) = new_seqtab_genus$Genus
new_seqtab_genus=new_seqtab_genus[,-1]
count_microscopy=count_microscopy[,common_columns,drop=FALSE]
count_microscopy = count_microscopy %>%
  filter(rowSums(. != 0, na.rm = TRUE) > 0)

#Save combined metadata/genus df and relative abundance df for genera 
write.table(metabar_df, "C:/Users/johan/OneDrive/R/Master project/merged_metabar_all_genera.tsv", sep="\t")
write.table(count_microscopy, "C:/Users/johan/OneDrive/R/Master project/count_microscopy_commonsamples.tsv", sep="\t")
write.table(new_seqtab_genus, "C:/Users/johan/OneDrive/R/Master project/new_seqtab_genus.tsv", sep="\t")

#Only common genera included

metabar_df=metabar_df[metabar_df$Genus %in% rownames(count_microscopy), , drop=FALSE] 
write.table(metabar_df, "C:/Users/johan/OneDrive/R/Master project/merged_metabar_common_genera.tsv", sep="\t")
new_seqtab_common_genus=new_seqtab_genus[rownames(new_seqtab_genus) %in% rownames(count_microscopy), , drop=FALSE]
write.table(new_seqtab_common_genus, "C:/Users/johan/OneDrive/R/Master project/new_seqtab_common_genus.tsv", sep="\t")
