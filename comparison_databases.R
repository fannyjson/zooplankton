
#Set wd, load libraries and files
setwd="C:/Users/johan/OneDrive/R/Master project"

library(readxl)
library(tidyr)
library(tidyverse)

smhi_data=read.delim("C:/Users/johan/OneDrive/R/Master project/zooplankton_2015_2020_2024-01-25_utf8.txt")
pr2_taxa=read_excel("C:/Users/johan/Downloads/pr2_version_5.0.0_taxonomy.xlsx")
smhi_data_common=read.delim("C:/Users/johan/OneDrive/R/Master project/smhi_data_common_samples_240516.tsv")
pr2_taxa_common=read.delim("C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_PR2_240516x.tsv")
silva_taxa_common=read.delim("C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_silva_240520.tsv")

silva_taxa=read.delim("C:/Users/johan/OneDrive/silva_headers_132_metazoa.csv", sep=";")


#################### ALL TAXA #######################################

#################### CLEAN PR2 FILE AND GENERATE DF WITH UNIQUE TAXLEVELS#######

pr2_taxa=pr2_taxa[pr2_taxa$subdivision=="Metazoa", , drop=FALSE]
pr2_taxa=pr2_taxa[, 1:9, drop=FALSE]
capitalize=function(string) {
  paste(toupper(substring(string, 1, 1)), substring(string, 2), sep = "")
}
colnames(pr2_taxa) = sapply(colnames(pr2_taxa), capitalize)


replace_patterns = function(vector) {
  # Replace _X and _XX with nothing
  vector = gsub("_X|_XX", "", vector)
  # Remove underscores
  vector = gsub("_", " ", vector)
  
  return(vector)
  
}

pr2_taxa = lapply(pr2_taxa, replace_patterns)
pr2_taxa = as.data.frame(pr2_taxa)

#Unique taxlevel PR2


pr2_taxa$Subdivision[pr2_taxa$Subdivision == ""] = NA
subdivision_pr2 = na.omit(unique(pr2_taxa$Subdivision))

pr2_taxa$Class[pr2_taxa$Class == ""] = NA
class_pr2 = na.omit(unique(pr2_taxa$Class))

pr2_taxa$Order[pr2_taxa$Order == ""] = NA
order_pr2 = na.omit(unique(pr2_taxa$Order))

pr2_taxa$Family[pr2_taxa$Family == ""] = NA
family_pr2 = na.omit(unique(pr2_taxa$Family))

pr2_taxa$Genus[pr2_taxa$Genus == ""] = NA
genus_pr2 = na.omit(unique(pr2_taxa$Genus))

pr2_taxa$Species[pr2_taxa$Species == ""] = NA
species_pr2 = na.omit(unique(pr2_taxa$Species))

# Find the maximum length
max_length = max(length(subdivision_pr2), 
                  length(class_pr2), length(order_pr2), length(family_pr2), length(genus_pr2), length(species_pr2))

# Make equal length by padding with NAs

subdivision_pr2 = c(subdivision_pr2, rep(NA, max_length - length(subdivision_pr2)))
class_pr2 = c(class_pr2, rep(NA, max_length - length(class_pr2)))
order_pr2 = c(order_pr2, rep(NA, max_length - length(order_pr2)))
family_pr2 = c(family_pr2, rep(NA, max_length - length(family_pr2)))
genus_pr2 = c(genus_pr2, rep(NA, max_length - length(genus_pr2)))
species_pr2 = c(species_pr2, rep(NA, max_length - length(species_pr2)))

# Combine into a single data frame

unique_pr2 = data.frame(
  Subdivision = subdivision_pr2,
  Class = class_pr2,
  Order = order_pr2,
  Family = family_pr2,
  Genus = genus_pr2,
  Species = species_pr2
)


################### CLEAN SILVA FILE AND GENERATE DF WITH UNIQUE TAXLEVELS ####

#Clean up SILVA reference file

silva_taxa = separate(silva_taxa, 
                      col = colnames(silva_taxa), 
                      into = paste0("col", 1:15), 
                      sep = ";")



silva_taxa=silva_taxa[,5:13, drop=FALSE]

colnames(silva_taxa)= c("Clade1", "Clade2", "Phylum", "Superclass", "Class", "Subclass", "Infraclass","Order", "Species")
#col=c("Class", "Order", "Species")
#silva_taxa=silva_taxa[,col, drop=FALSE]

#Move species in other column to last column
for (i in 1:(ncol(silva_taxa) - 1)) {
  two_words = grepl("^\\S+\\s+\\S+$", silva_taxa[[i]])
  silva_taxa[two_words, ncol(silva_taxa)] = silva_taxa[two_words, i]
}


silva_taxa=silva_taxa %>%
  mutate(Genus = str_extract(Species, "^[^ ]+")) %>%
  select(Clade1, Clade2, Phylum, Superclass, Class, Subclass, Infraclass, Order, Genus, Species)



replace_patterns = function(vector) {
  # Replace _X and _XX with nothing
  vector = gsub("_X|_XX", "", vector)

  # Remove underscores
  vector =gsub("_", " ", vector)
 
  
  vector =gsub(" \\([^()]+\\)", "", vector)
  
  vector =gsub(">", "", vector)
  
  
  return(vector)
}

silva_taxa = lapply(silva_taxa, replace_patterns)

###NEW

tax_levels = c("Clade1", "Clade2", "Phylum", "Superclass", "Class", "Subclass", "Infraclass", "Order", "Genus", "Species")

# Initialize an empty list to store unique values for each taxonomic level
unique_taxa = list()

# Loop through each taxonomic level, replace empty strings with NA and get unique non-NA values
for (level in tax_levels) {
  silva_taxa[[level]][silva_taxa[[level]] == ""] = NA
  unique_taxa[[level]] = na.omit(unique(silva_taxa[[level]]))
}

# Find the maximum length among all unique taxonomic levels
max_length = max(sapply(unique_taxa, length))

# Add NAs to make vectors of equal lengths
for (level in tax_levels) {
  unique_taxa[[level]] = c(unique_taxa[[level]], rep(NA, max_length - length(unique_taxa[[level]])))
}

# Create the final dataframe with unique taxonomic levels
unique_silva = data.frame(unique_taxa)

#Rename columns for clarity
colnames(unique_silva) = tax_levels


############# GENERATE DF WITH UNIQUE TAXLEVELS FOR MICROSCOPY (COMMON SAMPLES) ###########################################

#Only keep columns with taxonomic levels
cols_to_keep=c("taxon_class", "taxon_order", "taxon_family", "taxon_genus","taxon_species")
microscopy_taxa=smhi_data_common[ ,cols_to_keep, drop=FALSE] 

colnames(microscopy_taxa)[c(1:5)] = c("Class","Order", "Family", "Genus", "Species")

#Unique taxlevels

microscopy_taxa$Class[microscopy_taxa$Class == ""] = NA
class = na.omit(trimws(unique(microscopy_taxa$Class)))

microscopy_taxa$Order[microscopy_taxa$Order == ""] = NA
order = na.omit(trimws(unique(microscopy_taxa$Order)))

microscopy_taxa$Family[microscopy_taxa$Family == ""] = NA
family = na.omit(trimws(unique(microscopy_taxa$Family)))

microscopy_taxa$Genus[microscopy_taxa$Genus == ""] = NA
genus = na.omit(trimws(unique(microscopy_taxa$Genus)))

microscopy_taxa$Species[microscopy_taxa$Species == ""] = NA
species = na.omit(trimws(unique(microscopy_taxa$Species)))


max_length =  max(length(class), length(order), length(family), length(genus), length(species))


# Pad the shorter vectors with NA values to make them equal in length

class = c(class, rep(NA, max_length - length(class)))
order = c(order, rep(NA, max_length - length(order)))
family = c(family, rep(NA, max_length - length(family)))
genus = c(genus, rep(NA, max_length - length(genus)))
species = c(species, rep(NA, max_length - length(species)))

#Final df with unique taxlevels for microscopy 

unique_microscopy = data.frame(
                         Class=class,
                         Order = order,
                         Family = family,
                         Genus = genus,
                         Species = species)



unique_microscopy=mutate_all(unique_microscopy, ~na_if(., " "))

####### MATCH MICROSCOPY TAXA WITH PR2 TAXA ###################################


#For every column in the unique_microscopy df loop through the entire unique pr2 and find matches/
#mismatches


match_counts = matrix(0, nrow = ncol(unique_pr2), ncol = ncol(unique_microscopy), 
                       dimnames = list(colnames(unique_pr2), colnames(unique_microscopy)))

# Loop through each column in unique_silva, extract values from pr2 for comparison,
#count the number of matching cells between columns in microscopy and pr2

for (i in colnames(unique_pr2)) {
  for (j in colnames(unique_microscopy)) {
    pr2_values = as.character(unique_pr2[[i]])
    microscopy_values = as.character(unique_microscopy[[j]])
    matches = sum(pr2_values %in% microscopy_values & !is.na(pr2_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_microscopy=sapply(unique_microscopy, function(x) length(unique(na.omit(x))))
total_unique_pr2=sapply(unique_pr2, function(x) length(unique(na.omit(x))))

#Calculate unmatched counts
total_matches_microscopy=colSums(match_counts)
total_matches_pr2=rowSums(match_counts)
total_unmatched_microscopy=total_unique_microscopy - total_matches_microscopy
total_unmatched_pr2=total_unique_pr2 - total_matches_pr2

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched_microscopy)
match_counts=rbind(match_counts, total_unmatched_pr2)

rownames(match_counts)[7] = "Unmatched microscopy"
rownames(match_counts)[8] = "Unmatched PR2"

#Final match table

print(match_counts)



############### MATCH MICROSCOPY TAXA WITH SILVA TAXA #########################


# Create an empty match_counts matrix
match_counts = matrix(0, nrow = ncol(unique_silva), ncol = ncol(unique_microscopy), 
                      dimnames = list(colnames(unique_silva), colnames(unique_microscopy)))

# Loop through the columns of unique_silva and unique_microscopy
for (i in colnames(unique_silva)) {
  for (j in colnames(unique_microscopy)) {
    # Convert column values to character vectors
    silva_values = as.character(unique_silva[[i]])
    microscopy_values = as.character(unique_microscopy[[j]])
    # Count matches
    matches = sum(silva_values %in% microscopy_values & !is.na(silva_values) & !is.na(microscopy_values))
    # Store match count in the match_counts matrix
    match_counts[i, j] = matches
  }
}

# Calculate total unique values in each column of unique_microscopy and unique_silva
total_unique_microscopy = sapply(unique_microscopy, function(x) length(unique(na.omit(x))))
total_unique_silva = sapply(unique_silva, function(x) length(unique(na.omit(x))))

# Calculate total matches for each column in unique_microscopy
total_matches_silva = rowSums(match_counts)
total_matches_microscopy = colSums(match_counts)

# Calculate unmatched counts
total_unmatched_microscopy= total_unique_microscopy - total_matches_microscopy
total_unmatched_SILVA = total_unique_silva - total_matches_silva



# Add unmatched counts to the match_counts matrix
match_counts = rbind(match_counts, total_unmatched_microscopy)
match_counts = rbind(match_counts, total_unmatched_SILVA)

# Rename row names for clarity
rownames(match_counts)[nrow(match_counts)-1] = "Unmatched Microscopy"
rownames(match_counts)[nrow(match_counts)] = "Unmatched SILVA"

# Print the final match_counts table
print(match_counts)


################ MATCH SILVA WITH PR2 TAXA #####################################

match_counts = matrix(0, nrow = ncol(unique_silva), ncol = ncol(unique_pr2), 
                       dimnames = list(colnames(unique_silva), colnames(unique_pr2)))

# Loop through the columns of unique_silva and unique_pr2
for (i in colnames(unique_silva)) {
  for (j in colnames(unique_pr2)) {
    # Convert column values to character vectors
    silva_values = as.character(unique_silva[[i]])
    pr2_values = as.character(unique_pr2[[j]])
    # Count matches
    matches = sum(silva_values %in% pr2_values & !is.na(silva_values) & !is.na(pr2_values))
    # Store match count in the match_counts matrix
    match_counts[i, j] = matches
  }
}

# Calculate total unique values in each column of unique_pr2 and unique_silva
total_unique_pr2 = sapply(unique_pr2, function(x) length(unique(na.omit(x))))
total_unique_silva = sapply(unique_silva, function(x) length(unique(na.omit(x))))

# Calculate total matches for each column in unique_pr2
total_matches_silva = colSums(match_counts)
total_matches_pr2 = rowSums(match_counts)

# Calculate unmatched counts
total_unmatched_pr2 = total_unique_pr2 - total_matches_pr2
total_unmatched_SILVA = total_unique_silva - total_matches_silva

# Set negative unmatched counts to zero (optional, if needed)
total_unmatched_pr2[total_unmatched_pr2 < 0] = 0
total_unmatched_SILVA[total_unmatched_SILVA < 0] = 0

# Add unmatched counts to the match_counts matrix
match_counts = rbind(match_counts, total_unmatched_pr2)
match_counts = rbind(match_counts, total_unmatched_SILVA)

# Rename row names for clarity
rownames(match_counts)[nrow(match_counts)-1] = "Unmatched PR2"
rownames(match_counts)[nrow(match_counts)] = "Unmatched SILVA"

# Print the final match_counts table
print(match_counts)
###OLD 
match_counts = matrix(0, nrow = ncol(unique_pr2), ncol = ncol(unique_silva), 
                       dimnames = list(colnames(unique_pr2), colnames(unique_silva)))


for (i in colnames(unique_pr2)) {
  for (j in colnames(unique_silva)) {
    pr2_values = as.character(unique_pr2[[i]])
    silva_values = as.character(unique_silva[[j]])  # Corrected to silva_values
    matches = sum(pr2_values %in% silva_values & !is.na(pr2_values) & !is.na(silva_values))  # Corrected to silva_values
    match_counts[i, j] = matches
  }
}

total_unique_silva = sapply(unique_silva, function(x) length(unique(na.omit(x))))  # Corrected to sapply

# Calculate unmatched counts
total_matches = colSums(match_counts)
total_unmatched_sil = total_unique_silva - total_matches  # Corrected to total_unique_silva

total_unique_pr2 = sapply(unique_pr2, function(x) length(unique(na.omit(x))))

total_unmatched_pr = total_unique_pr2 - total_matches

# Add unmatched counts to table
total_unmatched_pr[total_unmatched_pr < 0] = 0
total_unmatched_sil[total_unmatched_sil < 0] = 0

match_counts = rbind(match_counts, total_unmatched_pr)
match_counts = rbind(match_counts, total_unmatched_sil)


rownames(match_counts)[6] = "Unmatched PR2"
rownames(match_counts)[7] = "Unmatched SILVA"

print(match_counts)

#################### ONLY TAXA IDENTIFIED IN COMMON SAMPLES ############################

########### PR2 TAXA FOR COMMON SAMPLES #################################################

pr2_taxa_common=pr2_taxa_common[,29:37, drop=FALSE]

pr2_taxa_common$Class[pr2_taxa_common$Class == ""] = NA
class_pr2_common = na.omit(unique(pr2_taxa_common$Class))

pr2_taxa_common$Order[pr2_taxa_common$Order == ""] = NA
order_pr2_common = na.omit(unique(pr2_taxa_common$Order))

pr2_taxa_common$Family[pr2_taxa_common$Family == ""] = NA
family_pr2_common = na.omit(unique(pr2_taxa_common$Family))

pr2_taxa_common$Genus[pr2_taxa_common$Genus == ""] = NA
genus_pr2_common = na.omit(unique(pr2_taxa_common$Genus))

pr2_taxa_common$Species[pr2_taxa_common$Species == ""] = NA
species_pr2_common = na.omit(unique(pr2_taxa_common$Species))


max_length = max(length(class_pr2_common), length(order_pr2_common), length(family_pr2_common), length(genus_pr2_common), length(species_pr2_common))

# Pad the shorter vectors with NA values to make them equal in length
class_pr2_common = c(class_pr2_common, rep(NA, max_length - length(class_pr2_common)))
order_pr2_common = c(order_pr2_common, rep(NA, max_length - length(order_pr2_common)))
family_pr2_common = c(family_pr2_common, rep(NA, max_length - length(family_pr2_common)))
genus_pr2_common = c(genus_pr2_common, rep(NA, max_length - length(genus_pr2_common)))
species_pr2_common = c(species_pr2_common, rep(NA, max_length - length(species_pr2_common)))

#Df with all the unique taxa (columns are not related)

unique_pr2_common = data.frame(Class = class_pr2_common,
                        Order = order_pr2_common,
                        Family = family_pr2_common,
                        Genus = genus_pr2_common,
                        Species = species_pr2_common)


unique_pr2_common=mutate_all(unique_pr2_common, ~na_if(., " "))

##PR2 MICROSCOPY FOR COMMON SAMPLES

match_counts = matrix(0, nrow = ncol(unique_pr2_common), ncol = ncol(unique_microscopy), 
                      dimnames = list(colnames(unique_pr2_common), colnames(unique_microscopy)))

for (i in colnames(unique_pr2_common)) {
  for (j in colnames(unique_microscopy)) {
    pr2_values = as.character(unique_pr2_common[[i]])
    microscopy_values = as.character(unique_microscopy[[j]])
    matches = sum(pr2_values %in% microscopy_values & !is.na(pr2_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_microscopy=sapply(unique_microscopy, function(x) length(unique(na.omit(x))))
total_unique_pr2=sapply(unique_pr2_common, function(x) length(unique(na.omit(x))))

#Calculate unmatched counts
total_matches=colSums(match_counts)
total_unmatched=total_unique_microscopy - total_matches
total_unmatched_pr2=total_unique_pr2 - total_matches

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched)
match_counts=rbind(match_counts, total_unmatched_pr2)

rownames(match_counts)[6] = "Unmatched microscopy"
rownames(match_counts)[7] = "Unmatched PR2"

print(match_counts)

###GENERATE LATEX TABLES 
t(match_counts)
table = kable(match_counts, format = "html", align = "c", escape = FALSE) %>%
  kable_styling(full_width = FALSE, font_size = 10, latex_options = "scale_down") %>%
  column_spec(1:4, width = "auto")

print(table)
