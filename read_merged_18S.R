# PRE-PROCESSING OF DATA (KRZYZTOFS MODIFIED CODE)
# SET WORKING DIRECTORY, LOAD LIBRARIES, SET IN FILES

setwd("C:/Users/johan/OneDrive/R")

library(tidyverse)
library(vegan)
library(geosphere)


in_folder = paste('C:/Users/johan/OneDrive/R/Master project', '/', sep = '')
asv_seq_file = paste(in_folder, 'asvs_18S', '.txt', sep = '') ## ASV sequences
seqtab_file =  paste(in_folder, 'OTUs_18S','.txt', sep = '') ## count matrix for each ASV & sample
taxa_file =  paste(in_folder, 'taxa_18S', '.tsv', sep = '') ## taxonomic annotation for each ASV
metadata_file = paste(in_folder,'metadata.tsv', sep ='')  ## metadata file INCLUDING physiochemical data
plot_folder <- paste0(in_folder, 'plots/') ## folder to store the plots
if(!dir.exists(plot_folder)){dir.create(plot_folder)}

# READ 18S DATA

asv = as.vector(read.csv(asv_seq_file, sep = '\t')[[2]])
names(asv) = as.vector(read.csv(asv_seq_file, sep = '\t')[[1]])
asv = toupper(asv)

seqtab = as.matrix(read.delim(seqtab_file, row.names = 'OTU_ID'))

taxa = as.matrix(read.delim(taxa_file, sep = '\t', row.names = 1)) # modified using correct_taxa.R (from taxa_original.txt)

metadata = read_tsv(metadata_file)

## Pick the taxonomic annotation and sequences only of the ASVs representative of dbOTUs

ix_taxa = match(rownames(seqtab), rownames(taxa))
taxa = taxa[ix_taxa,]
ix_taxa = match(rownames(seqtab), names(asv))
asv = asv[ix_taxa]

## Adjust seqtab colnames

colnames(seqtab) = gsub(colnames(seqtab), pattern = '_18S', replacement = '')
colnames(seqtab) = gsub(colnames(seqtab), pattern = '^X', replacement = '')

## Pick metadata only for the monitoring samples which were not technical tests
## But also volume test 500ml, since that corresponds to the sampling for monioting

ix_monitoring = which(metadata$group == 'Monitoring')
ix_volume_test_500ml = intersect(which(metadata$group == 'volume_test'),
                                 which(metadata$sampled_volume == 500))

ix = union(ix_monitoring, ix_volume_test_500ml)
metadata = metadata[ix,]

## Keep only high-quality samples
## Also: adjust seqtab and metadata (after picking montioring samples) to keep only the samples present in both
## Also; remove replicates

ix = which(metadata$sample_id %in% colnames(seqtab))
metadata = metadata[ix,]

## Pick only columns from seqtab which correspond to the monitoring samples

iy = which(colnames(seqtab) %in% metadata$sample_id)
seqtab = seqtab[,iy]

# And set the same order of samples in metadata and seqtab
iy = match(metadata$sample_id, colnames(seqtab))
seqtab = seqtab[,iy]

#check
all(colnames(seqtab) == metadata$sample_id)

## Exclude samples which are not high-quality

samples_to_exclude = c('20191015_263953_1') ## suspected contamination
samples_to_exclude = c(samples_to_exclude, '20190820_264450_1') ## undersequenced samples with replicates available

## Choose only one replicate from each date and station 
dates_stations = c()
meta_ix = setdiff(1:nrow(metadata), which(metadata$sample_id %in% samples_to_exclude))

for(i in meta_ix){
  ds = paste(metadata$date[i], metadata$station_name[i])
  if(! ds %in% dates_stations){
    dates_stations = c(dates_stations, ds)
  }else{
    meta_ix = meta_ix[! meta_ix %in% i]
  }
}


# Pick only selected stations for further diversity analyses
metadata = metadata[meta_ix,]
ix = metadata$sample_id
seqtab = seqtab[,ix]

## Change NB1 / 3 to B7
## Since they are practically the same sampling location

metadata$station_name[metadata$station_name == 'NB1 / B3'] = 'B7'

## Remove spike 
spike = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'

spike_ix = agrep(spike, asv)

ix_taxa =  setdiff(1:nrow(taxa), spike_ix)

## Excluding non-anotated ASVs
ix_taxa = intersect(ix_taxa, which(complete.cases(taxa[,2])))

## Include only metazoa
#ix_taxa = which(taxa[, 4] == 'Metazoa')
#taxa = taxa[ix_taxa,]
#seqtab = seqtab[ix_taxa,]
#asv = asv[ix_taxa]



### Analyze nr of reads per sample ###
### Exclude under/overssequenced samples ###
## Different thresholds across the years, since different sequencing technologies were used (wore reads in general in 2015-2017 samples)

samples_to_exclude = c()

ix_2019 = which(metadata$date > as.Date('2019-01-01'))

png(paste0(plot_folder, 'tot_reads_sorted_18S_2019.png'),
    width = 5000, height = 3000, res = 500, units = 'px')
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjusted margin for title
plot(sort(colSums(seqtab)[ix_2019]), xlab = 'Samples', ylab = 'Reads', main = 'Total reads for samples after 2019 (All 18S)')

dev.off()

sort(colSums(seqtab)[ix_2019])

ix = which(colSums(seqtab)[ix_2019] > 200000)

samples_to_exclude = c(samples_to_exclude, names(colSums(seqtab)[ix_2019][ix]))

ix_2015 = which(metadata$date < as.Date('2019-01-01'))

png(paste0(plot_folder, 'tot_reads_sorted_18S_2015.png'),
    width = 5000, height = 3000, res = 500, units = 'px')
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjusted margin for title
plot(sort(colSums(seqtab)[ix_2015]), xlab = 'Samples', ylab = 'Reads', main = 'Total reads for samples before 2019 (All 18S')
dev.off()

sort(colSums(seqtab)[ix_2015])

## Remove taxa with 0 reads across all finally chosen samples

which(rowSums(seqtab) == 0)

ix_taxa = which(rowSums(seqtab) != 0)

taxa = taxa[ix_taxa,]
seqtab = seqtab[ix_taxa,]
asv = asv[ix_taxa]

## Get normalised count table - i.e. relative abundances
## Used for relative abundances of taxa

norm_seqtab = seqtab
for (i in 1:ncol(seqtab)) {
  norm_seqtab[,i] = seqtab[,i]/sum(seqtab[,i])
}


### Rarefaction ###

## Specify parameters
## Sample size

m = min(colSums(seqtab))
print(m)

## NR of interations

n_it = 100

## Plot rarefaction curves
ix_2019 = which(metadata$date > as.Date('2019-01-01'))
  
  png(paste0(plot_folder, 'rarecurve_2019.png'),
      width = 5000, height = 3000, res = 500, units = 'px')
  # pdf(paste0(plot_folder, 'rarecurve_', data_type, '.pdf'),
  #     width = 10, height = 6)
  
  par(mar = c(4,4,2,2))
  rarecurve(t(seqtab[,ix_2019]), step = 1000, label = FALSE, xlab = 'Sequencing depth', ylab = 'Samples', main='Rarefaction for samples after 2019 (All 18S)')
  abline(v = m, col = "red", lty = 2)
  
  dev.off()
  
  print(sort(rareslope(t(seqtab[,ix_2019]), m)))
  
  ix_2015 = which(metadata$date > as.Date('2015-01-01'))
  
  png(paste0(plot_folder, 'rarecurve_2015.png'),
      width = 5000, height = 3000, res = 500, units = 'px')
  # pdf(paste0(plot_folder, 'rarecurve_', data_type, '.pdf'),
  #     width = 10, height = 6)
  
  par(mar = c(4,4,2,2))
  rarecurve(t(seqtab[,ix_2015]), step = 1000, label = FALSE, xlab = 'Sequencing depth', ylab = 'Samples' , main='Rarefaction for samples before 2019 (All 18S)')
  abline(v = m, col = "red", lty = 2)
  
  dev.off()
  
  print(sort(rareslope(t(seqtab[,ix_2015]), m)))
  print(sort(rareslope(t(seqtab), m)))


  ### Longitude, Latitude, and Salinity -- assume station's mean if not available ###
  
  df_sum_up = aggregate.data.frame(metadata, list(metadata$station_name), mean, na.rm = TRUE)
  
  ## Order stations by Salinity (like other groups)
  
  df_sum_up = df_sum_up[order(-df_sum_up$Salinity),]
  
  ## Get missing Latitude/Longitude data from other datapoints
  ## for the same station (assume stations do not move)
  
  for(ind in seq(1,nrow(metadata))){
    if(is.na(metadata$Longitude[ind])){
      ix = which(df_sum_up$Group.1 == metadata$station_name[ind])
      metadata$Longitude[ind] = df_sum_up$Longitude[ix]
    }
    if(is.na(metadata$Latitude[ind])){
      ix = which(df_sum_up$Group.1 == metadata$station_name[ind])
      metadata$Latitude[ind] = df_sum_up$Latitude[ix]
    }
  }
  
  
  ## Get average Salinity for each station, fill Salinity gaps for lacking datapoints
  
  for(ind in seq(1,nrow(metadata))){
    if(is.na(metadata$Salinity[ind])){
      ix = which(df_sum_up$Group.1 == metadata$station_name[ind])
      metadata$Salinity[ind] = df_sum_up$Salinity[ix]
    }
  }
  
  # Add Sea basin to metadata, code provided by Bengt Karlsson
  
  metadata <- metadata %>%
    mutate(sea_basin = station_name)
  
  metadata$sea_basin<-gsub("RÅNEÅ-1", "Bothnian Bay", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("RÅNEÅ-2", "Bothnian Bay", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("F9 / A13", "Bothnian Bay", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("B7", "Bothnian Sea", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("NB1 / B3", "Bothnian Sea", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("C3", "Bothnian Sea", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("GAVIK-1", "Bothnian Sea", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("SR3", "Bothnian Sea", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("B1", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("H4", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BY29 / LL19", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BY31 LANDSORTSDJ", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BY15 GOTLANDSDJ", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("REF M1V1", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BY2 ARKONA", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BY5 BORNHOLMSDJ", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BCS III-10", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("BY38 KARLSÖDJ", "Baltic Proper", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("ANHOLT E", "Kattegat", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("N14 FALKENBERG", "Kattegat", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("SLÄGGÖ", "Skagerrak", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  metadata$sea_basin<-gsub("Å17", "Skagerrak", metadata$sea_basin, ignore.case = FALSE, perl = FALSE,
                           fixed = TRUE, useBytes = FALSE)
  
  ### Get Sea Basin stats table
  
  sea_basins = c("Bothnian Bay", "Bothnian Sea", "Baltic Proper", "Kattegat", "Skagerrak")

    sea_basin_stats = data.frame(sea_basin = sea_basins, Salinity = NA, Latitude = NA, Temperature = NA)
    
    for(sb in sea_basins){
      ix = which(metadata$sea_basin == sb)
      mn_sal = mean(metadata$Salinity[ix])
      mn_lat = mean(metadata$Latitude[ix])
      mn_temp = mean(na.omit(metadata$Temperature[ix]))
      i = which(sea_basin_stats$sea_basin == sb)
      sea_basin_stats$Salinity[i] = mn_sal
      sea_basin_stats$Latitude[i] = mn_lat
      sea_basin_stats$Temperature[i] = mn_temp
    }
    
  
  
  ## Annotate the samples to respective seasons
  
  ## adopted from https://stackoverflow.com/questions/36502140/determine-season-from-date-using-lubridate-in-r
  
  seasons = c("Winter","Spring","Summer","Fall")
  
  getSeason <- function(input.date){
    numeric.date <- 100*month(input.date)+day(input.date)
    # numeric.date <- month(input.date)
    ## input Seasons upper limits in the form MMDD in the "break =" option:
    cuts <- base::cut(numeric.date, breaks = c(0,319,0620,0921,1220,1231))
    # cuts <- base::cut(numeric.date, breaks = c(0,2,5,8,11,12))
    # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
    levels(cuts) <- c(seasons, seasons[1])
    return(cuts)
  }
  
  metadata$season = getSeason(metadata$date)
  
  ## Annotate the samples to respective months
  
  metadata$month = months(metadata$date)
  
  ## Get daylength
  
  metadata$day_length = daylength(metadata$Latitude, yday(metadata$date))
  
  ## Set ploting theme and colors
  
  ## By default the plots will now have white background
  my_theme = theme_bw() +
    theme(
      text = element_text(size = 16),  # Adjust font size for all text
      plot.title = element_text(size = 18),  # Adjust font size for plot titles
      axis.title = element_text(size = 16)  # Adjust font size for axis labels
    )
  
 
    max_tax_l = ncol(taxa)
    bar_col = '#1c9099'
  
  basin_colors = c('#4d004b', '#3690c0', '#238443', '#fec44f', '#f768a1')
  sea_basin_shape_scale = c(25, 24, 23, 22, 21)
  names(sea_basin_shape_scale) = sea_basins
  names(basin_colors) = sea_basins
  
  ## Save alfa- and beta-diversity data
  
#write.table(metadata,
              #file = paste(in_folder, 'metadata_240308.tsv', sep = ''),
              #sep = '\t', row.names = FALSE, col.names = TRUE)
 
write.table(norm_seqtab,
                file = paste(in_folder, 'norm_seqtab_18S_all_240308.tsv', sep = ''),
                sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(seqtab, 
              file = paste(in_folder, 'seqtab_18S_all_240308.tsv', sep = ''),
              sep = '\t', row.names = TRUE, col.names = TRUE)

save.image(paste('read_merged_18S_240308.RData', sep = '')) 
  
