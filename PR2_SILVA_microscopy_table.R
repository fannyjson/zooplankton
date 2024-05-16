
#Set wd, load libraries and files
setwd="C:/Users/johan/OneDrive/R/Master project"

library(readxl)
library(tidyr)
library(tidyverse)

smhi_data=read.delim("C:/Users/johan/OneDrive/R/Master project/zooplankton_2015_2020_2024-01-25_utf8.txt")
pr2_taxa=read_excel("C:/Users/johan/Downloads/pr2_version_5.0.0_taxonomy.xlsx")
silva_taxa=read.delim("C:/Users/johan/Downloads/silva_132_headers.txt")
smhi_data_common=read.delim("C:/Users/johan/OneDrive/R/Master project/smhi_data_common_samples_240411.tsv")
pr2_taxa_common=read.delim("C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_PR2_240411.tsv")
silva_taxa_common=read.delim("C:/Users/johan/OneDrive/R/Master project/metabar_df_common_samples_silva_240413.tsv")

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

class_pr2=unique(pr2_taxa$Class)
class_pr2=na.omit(class_pr2)

order_pr2=unique(pr2_taxa$Order)
order_pr2=na.omit(order_pr2)

family_pr2=unique(pr2_taxa$Family)
family_pr2=na.omit(family_pr2)

genus_pr2=unique(pr2_taxa$Genus)
genus_pr2=na.omit(genus_pr2)

species_pr2=unique(pr2_taxa$Species)
species_pr2=na.omit(species_pr2)

max_length = max(length(class_pr2), length(order_pr2), length(family_pr2), length(genus_pr2), length(species_pr2))

# Pad the shorter vectors with NA values to make them equal in length
class_pr2 = c(class_pr2, rep(NA, max_length - length(class_pr2)))
order_pr2 = c(order_pr2, rep(NA, max_length - length(order_pr2)))
family_pr2 = c(family_pr2, rep(NA, max_length - length(family_pr2)))
genus_pr2 = c(genus_pr2, rep(NA, max_length - length(genus_pr2)))
species_pr2 = c(species_pr2, rep(NA, max_length - length(species_pr2)))

#df with all the unique taxa (columns are NOT related)

unique_pr2 = data.frame(Class = class_pr2,
                        Order = order_pr2,
                        Family = family_pr2,
                        Genus = genus_pr2,
                        Species = species_pr2)

################### CLEAN SILVA FILE AND GENERATE DF WITH UNIQUE TAXLEVELS ####

silva_taxa = separate(silva_taxa, 
                      col = colnames(silva_taxa), 
                      into = paste0("col", 1:7), 
                      sep = ";")

colnames(silva_taxa)[c(1, 2, 3, 4, 5, 6, 7)] = c("Kingdom" , "Phylum", "Class", "Order", "Family", "Genus", "Species")

silva_taxa=silva_taxa[,3:7, drop=FALSE]

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
silva_taxa=as.data.frame(silva_taxa)
silva_taxa=silva_taxa[silva_taxa$Order=="Metazoa", , drop=FALSE]

#Unique taxlevels 

class_silva=unique(silva_taxa$Class)
class_silva=na.omit(class_silva)

order_silva=unique(silva_taxa$Order)
order_silva=na.omit(order_silva)

family_silva=unique(silva_taxa$Family)
family_silva=na.omit(family_silva)

genus_silva=unique(silva_taxa$Genus)
genus_silva=na.omit(genus_silva)

species_silva=unique(silva_taxa$Species)
species_silva=na.omit(species_silva)

max_length = max(length(class_silva), length(order_silva), length(family_silva), length(genus_silva), length(species_silva))


#Insert NAs to make vectors of equal lengths

class_silva = c(class_silva, rep(NA, max_length - length(class_silva)))
order_silva = c(order_silva, rep(NA, max_length - length(order_silva)))
family_silva = c(family_silva, rep(NA, max_length - length(family_silva)))
genus_silva = c(genus_silva, rep(NA, max_length - length(genus_silva)))
species_silva = c(species_silva, rep(NA, max_length - length(species_silva)))

#Final dataframe with unique taxlevels

unique_silva = data.frame(Class = class_silva,
                          Order = order_silva,
                          Family = family_silva,
                          Genus = genus_silva,
                          Species = species_silva)



############# GENERATE DF WITH UNIQUE TAXLEVELS FOR MICROSCOPY ###########################################


#Only keep columns with taxonomic levels
cols_to_keep=c("taxon_class", "taxon_order", "taxon_family", "taxon_genus","taxon_species")
microscopy_taxa=smhi_data[ ,cols_to_keep, drop=FALSE]

colnames(microscopy_taxa)[c(1:5)] = c("Class","Order", "Family", "Genus", "Species")

#Unique taxlevels

class=trimws(unique(microscopy_taxa$Class))
class=na.omit(class)

order=trimws(unique(microscopy_taxa$Order))
order=na.omit(order)

family=trimws(unique(microscopy_taxa$Family))
family=na.omit(family)

genus=trimws(unique(microscopy_taxa$Genus))
genus=na.omit(genus)

species=trimws(unique(microscopy_taxa$Species))
species=na.omit(species)

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

total_unique_microscopy=sapply(unique_microscopy, function(x) length(unique(x)))
total_unique_pr2=sapply(unique_pr2, function(x) length(unique(x)))

#Calculate unmatched counts
total_matches=colSums(match_counts)
total_unmatched=total_unique_microscopy - total_matches
total_unmatched_pr2=total_unique_pr2 - total_matches

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched)
match_counts=rbind(match_counts, total_unmatched_pr2)

rownames(match_counts)[6] = "Unmatched microscopy"
rownames(match_counts)[7] = "Unmatched PR2"

#Final match table

print(match_counts)


############### MATCH MICROSCOPY TAXA WITH SILVA TAXA #########################


match_counts = matrix(0, nrow = ncol(unique_silva), ncol = ncol(unique_microscopy), 
                       dimnames = list(colnames(unique_silva), colnames(unique_microscopy)))

for (i in colnames(unique_silva)) {
  for (j in colnames(unique_microscopy)) {
    silva_values = as.character(unique_silva[[i]])
    microscopy_values = as.character(unique_microscopy[[j]])
    matches = sum(silva_values %in% microscopy_values & !is.na(silva_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_microscopy=sapply(unique_microscopy, function(x) length(unique(x)))
total_unique_silva=sapply(unique_silva, function(x) length(unique(x)))

#Calculate unmatched counts
total_matches=colSums(match_counts)
total_unmatched=total_unique_microscopy - total_matches
unmatched_SILVA=total_unique_silva - total_matches

unmatched_SILVA <- total_unique_silva - total_matches
unmatched_SILVA[unmatched_SILVA < 0] = 0

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched)
match_counts=rbind(match_counts, unmatched_SILVA)

rownames(match_counts)[6] = "Unmatched Microscopy"
rownames(match_counts)[7] = "Unmatched SILVA"

#Final match table

print(match_counts)

################ MATCH SILVA WITH PR2 TAXA #####################################

match_counts = matrix(0, nrow = ncol(unique_pr2), ncol = ncol(unique_silva), 
                       dimnames = list(colnames(unique_pr2), colnames(unique_silva)))

for (i in colnames(unique_pr2)) {
  for (j in colnames(unique_silva)) {
    silva_values = as.character(unique_pr2[[i]])
    microscopy_values = as.character(unique_silva[[j]])
    matches = sum(silva_values %in% microscopy_values & !is.na(silva_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_silva=sapply(unique_silva, function(x) length(unique(x)))

#Calculate unmatched counts
total_matches=colSums(match_counts)
total_unmatched_sil=total_unique_silva - total_matches
total_unmatched_pr=total_unique_pr2 - total_matches

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched_sil)
match_counts=rbind(match_counts, total_unmatched_pr)

rownames(match_counts)[6] = "Unmatched PR2"
rownames(match_counts)[7] = "Unmatched SILVA"

#Final match table 
print(match_counts)


#################### ONLY TAXA IDENTIFIED IN COMMON SAMPLES ############################

########### PR2 TAXA FOR COMMON SAMPLES #################################################

pr2_taxa_common=pr2_taxa_common[,29:37, drop=FALSE]

class_pr2_common=unique(pr2_taxa_common$Class)
class_pr2_common=na.omit(class_pr2_common)

order_pr2_common=unique(pr2_taxa_common$Order)
order_pr2_common=na.omit(order_pr2_common)

family_pr2_common=unique(pr2_taxa_common$Family)
family_pr2_common=na.omit(family_pr2_common)

genus_pr2_common=unique(pr2_taxa_common$Genus)
genus_pr2_common=na.omit(genus_pr2_common)

species_pr2_common=unique(pr2_taxa_common$Species)
species_pr2_common=na.omit(species_pr2_common)

max_length = max(length(class_pr2_common), length(order_pr2_common), length(family_pr2_common), length(genus_pr2_common), length(species_pr2_common))

# Pad the shorter vectors with NA values to make them equal in length
class_pr2_common = c(class_pr2_common, rep(NA, max_length - length(class_pr2_common)))
order_pr2_common = c(order_pr2_common, rep(NA, max_length - length(order_pr2_common)))
family_pr2_common = c(family_pr2_common, rep(NA, max_length - length(family_pr2_common)))
genus_pr2_common = c(genus_pr2_common, rep(NA, max_length - length(genus_pr2_common)))
species_pr2_common = c(species_pr2_common, rep(NA, max_length - length(species_pr2_common)))

#df with all the unique taxa (columns are NOT related)

unique_pr2_common = data.frame(Class = class_pr2_common,
                        Order = order_pr2_common,
                        Family = family_pr2_common,
                        Genus = genus_pr2_common,
                        Species = species_pr2_common)


########### MICROSCOPY TAXA FOR COMMON SAMPLES #################################################

cols_to_keep=c("taxon_class", "taxon_order", "taxon_family", "taxon_genus","taxon_species")
microscopy_taxa_common=smhi_data_common[ ,cols_to_keep, drop=FALSE]

colnames(microscopy_taxa_common)[c(1:5)] = c("Class","Order", "Family", "Genus", "Species")

class_common=trimws(unique(microscopy_taxa_common$Class))
class_common=na.omit(class_common)

order_common=trimws(unique(microscopy_taxa_common$Order))
order_common=na.omit(order_common)

family_common=trimws(unique(microscopy_taxa_common$Family))
family_common=na.omit(family_common)

genus_common=trimws(unique(microscopy_taxa_common$Genus))
genus_common=na.omit(genus_common)

species_common=trimws(unique(microscopy_taxa_common$Species))
species_common=na.omit(species_common)

max_length =  max(length(class_common), length(order_common), length(family_common), length(genus_common), length(species_common))

# Pad the shorter vectors with NA values to make them equal in length

class_common = c(class_common, rep(NA, max_length - length(class_common)))
order_common = c(order_common, rep(NA, max_length - length(order_common)))
family_common = c(family_common, rep(NA, max_length - length(family_common)))
genus_common = c(genus_common, rep(NA, max_length - length(genus_common)))
species_common = c(species_common, rep(NA, max_length - length(species_common)))

#Final df with unique taxlevels for microscopy 

unique_microscopy_common= data.frame(
  Class=class_common,
  Order = order_common,
  Family = family_common,
  Genus = genus_common,
  Species = species_common)

match_counts = matrix(0, nrow = ncol(unique_pr2_common), ncol = ncol(unique_microscopy_common), 
                      dimnames = list(colnames(unique_pr2_common), colnames(unique_microscopy_common)))

# Loop through each column in unique_silva, extract values from pr2 for comparison,
#count the number of matching cells between columns in microscopy and pr2

for (i in colnames(unique_pr2_common)) {
  for (j in colnames(unique_microscopy_common)) {
    pr2_values = as.character(unique_pr2_common[[i]])
    microscopy_values = as.character(unique_microscopy_common[[j]])
    matches = sum(pr2_values %in% microscopy_values & !is.na(pr2_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_microscopy=sapply(unique_microscopy_common, function(x) length(unique(x)))
total_unique_pr2=sapply(unique_pr2_common, function(x) length(unique(x)))

#Calculate unmatched counts
total_matches=colSums(match_counts)
total_unmatched=total_unique_microscopy - total_matches
total_unmatched_pr2=total_unique_pr2 - total_matches

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched)
match_counts=rbind(match_counts, total_unmatched_pr2)

rownames(match_counts)[6] = "Unmatched microscopy"
rownames(match_counts)[7] = "Unmatched PR2"

#Final match table

print(match_counts)

########### SILVA TAXA FOR COMMON SAMPLES #################################################

silva_taxa_common$Order= gsub("Metazoa (Animalia)", "Metazoa", silva_taxa_common$Order, fixed = TRUE)

#Unique taxlevels 

class_silva_common=unique(silva_taxa_common$Class)
class_silva_common=na.omit(class_silva_common)

order_silva_common=unique(silva_taxa_common$Order)
order_silva_common=na.omit(order_silva_common)

family_silva_common=unique(silva_taxa_common$Family)
family_silva_common=na.omit(family_silva_common)

genus_silva_common=unique(silva_taxa_common$Genus)
genus_silva_common=na.omit(genus_silva_common)

species_silva_common=unique(silva_taxa_common$Species)
species_silva_common=na.omit(species_silva_common)

max_length = max(length(class_silva_common), length(order_silva_common), length(family_silva_common), length(genus_silva_common), length(species_silva_common))

# Pad the shorter vectors with NA values to make them equal in length

class_silva_common = c(class_silva_common, rep(NA, max_length - length(class_silva_common)))
order_silva_common = c(order_silva_common, rep(NA, max_length - length(order_silva_common)))
family_silva_common = c(family_silva_common, rep(NA, max_length - length(family_silva_common)))
genus_silva_common = c(genus_silva_common, rep(NA, max_length - length(genus_silva_common)))
species_silva_common = c(species_silva_common, rep(NA, max_length - length(species_silva_common)))

#Final df with unique taxlevels

unique_silva_common = data.frame(Class = class_silva_common,
                          Order = order_silva_common,
                          Family = family_silva_common,
                          Genus = genus_silva_common,
                          Species = species_silva_common)

#Make match table

match_counts = matrix(0, nrow = ncol(unique_pr2_common), ncol = ncol(unique_silva_common), 
                      dimnames = list(colnames(unique_pr2_common), colnames(unique_silva_common)))

for (i in colnames(unique_pr2_common)) {
  for (j in colnames(unique_silva_common)) {
    silva_values = as.character(unique_silva_common[[i]])
    pr2_values = as.character(unique_pr2_common[[j]])
    matches = sum(silva_values %in% pr2_values & !is.na(pr2_values) & !is.na(silva_values))
    match_counts[i, j] = matches
  }
}

total_unique_silva_common=sapply(unique_silva_common, function(x) length(unique(x)))
total_unique_pr2_common=sapply(unique_pr2_common, function(y) length(unique(y)))

#Calculate unmatched counts
total_matches=colSums(match_counts)

unmatched_SILVA <- total_unique_silva_common - total_matches
unmatched_SILVA[unmatched_SILVA < 0] = 0

unmatched_pr <- total_unique_pr2_common - total_matches
unmatched_SILVA[unmatched_pr < 0] = 0




#Add unmatched counts to table
match_counts=rbind(match_counts, unmatched_pr)
match_counts=rbind(match_counts, unmatched_SILVA)

rownames(match_counts)[6] = "Unmatched PR2"
rownames(match_counts)[7] = "Unmatched SILVA"

#Final match table 
print(match_counts)


#PR2 SILVA FOR COMMON SAMPLES
match_counts = matrix(0, nrow = ncol(unique_pr2_common), ncol = ncol(unique_microscopy_common), 
                      dimnames = list(colnames(unique_pr2_common), colnames(unique_microscopy_common)))

# Loop through each column in unique_silva, extract values from pr2 for comparison,
#count the number of matching cells between columns in microscopy and pr2

for (i in colnames(unique_pr2_common)) {
  for (j in colnames(unique_microscopy_common)) {
    pr2_values = as.character(unique_pr2_common[[i]])
    microscopy_values = as.character(unique_microscopy_common[[j]])
    matches = sum(pr2_values %in% microscopy_values & !is.na(pr2_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_microscopy=sapply(unique_microscopy_common, function(x) length(unique(x)))
total_unique_pr2=sapply(unique_pr2_common, function(x) length(unique(x)))

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

#SILVA AND MICROSCOPY COMMON SAMPLES 
match_counts = matrix(0, nrow = ncol(unique_silva_common), ncol = ncol(unique_microscopy_common), 
                      dimnames = list(colnames(unique_silva_common), colnames(unique_microscopy_common)))

for (i in colnames(unique_silva_common)) {
  for (j in colnames(unique_microscopy_common)) {
    silva_values = as.character(unique_silva_common[[i]])
    microscopy_values = as.character(unique_microscopy_common[[j]])
    matches = sum(silva_values %in% microscopy_values & !is.na(silva_values) & !is.na(microscopy_values))
    match_counts[i, j] = matches
  }
}

total_unique_microscopy=sapply(unique_microscopy_common, function(x) length(unique(x)))
total_unique_silva=sapply(unique_silva_common, function(x) length(unique(x)))

#Calculate unmatched counts
total_matches=colSums(match_counts)
total_unmatched=total_unique_microscopy - total_matches
unmatched_SILVA=total_unique_silva - total_matches

unmatched_SILVA <- total_unique_silva_common - total_matches
unmatched_SILVA[unmatched_SILVA < 0] = 0

#Add unmatched counts to table
match_counts=rbind(match_counts, total_unmatched)
match_counts=rbind(match_counts, unmatched_SILVA)

rownames(match_counts)[6] = "Unmatched Microscopy"
rownames(match_counts)[7] = "Unmatched SILVA"

#Final match table

print(match_counts)