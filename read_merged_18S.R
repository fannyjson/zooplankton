# PRE-PROCESSING OF DATA (KRZYZTOFS MODIFIED CODE)
# SET WORKING DIRECTORY, LOAD LIBRARIES, SET IN FILES

setwd("C:/Users/johan/OneDrive/R")

library(tidyverse)
library(vegan)
library(geosphere)


in_folder = paste('C:/Users/johan/OneDrive/R/Test preprocessing', '/', sep = '')
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
ix_taxa = which(taxa[, 4] == 'Metazoa')
taxa = taxa[ix_taxa,]
seqtab = seqtab[ix_taxa,]
asv = asv[ix_taxa]



### Analyze nr of reads per sample ###
### Exclude under/overssequenced samples ###
## Different thresholds across the years, since different sequencing technologies were used (wore reads in general in 2015-2017 samples)

samples_to_exclude = c()

ix_2019 = which(metadata$date > as.Date('2019-01-01'))

png(paste0(plot_folder, 'tot_reads_sorted_18S', '_2019.png'),
    width = 5000, height = 3000, res = 500, units = 'px')
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjusted margin for title
plot(sort(colSums(seqtab)[ix_2019]), xlab = 'Samples', ylab = 'Reads', main = 'Total reads for samples after 2019')

dev.off()

sort(colSums(seqtab)[ix_2019])

ix = which(colSums(seqtab)[ix_2019] > 200000)

samples_to_exclude = c(samples_to_exclude, names(colSums(seqtab)[ix_2019][ix]))

ix_2015 = which(metadata$date < as.Date('2019-01-01'))

png(paste0(plot_folder, 'tot_reads_sorted_18S', '_2015.png'),
    width = 5000, height = 3000, res = 500, units = 'px')
par(mar = c(5, 4, 4, 2) + 0.1)  # Adjusted margin for title
plot(sort(colSums(seqtab)[ix_2015]), xlab = 'Samples', ylab = 'Reads', main = 'Total reads for samples before 2019')
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


## Sum up clade counts at each taxonomic level

clade_counts = list()
norm_clade_counts = list()

for (i in 1:ncol(taxa)) {
  matr = norm_matr = NULL
  clade = unique(taxa[,i])
  clade = clade[!is.na(clade)]
  for (j in 1:length(clade)) {
    ix = which(clade[j]==taxa[,i])
    if (length(ix) > 1) {
      matr = rbind(matr, apply(seqtab[ix,], 2, sum, na.rm=TRUE))
      norm_matr = rbind(norm_matr, apply(norm_seqtab[ix,], 2, sum, na.rm=TRUE))
    } else {
      matr = rbind(matr, seqtab[ix,])
      norm_matr = rbind(norm_matr, norm_seqtab[ix,])
    }
  }
  rownames(matr) = rownames(norm_matr) = clade
  colnames(matr) = colnames(norm_matr) = metadata$sample_id
  clade_counts[[i]] = matr
  norm_clade_counts[[i]] = norm_matr
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
  rarecurve(t(seqtab[,ix_2019]), step = 1000, label = FALSE, xlab = 'Sequencing depth', ylab = "Rarefied No. of dbOTUs")
  abline(v = m, col = "red", lty = 2)
  
  dev.off()
  
  print(sort(rareslope(t(seqtab[,ix_2019]), m)))
  
  ix_2015 = which(metadata$date > as.Date('2015-01-01'))
  
  png(paste0(plot_folder, 'rarecurve_2015.png'),
      width = 5000, height = 3000, res = 500, units = 'px')
  # pdf(paste0(plot_folder, 'rarecurve_', data_type, '.pdf'),
  #     width = 10, height = 6)
  
  par(mar = c(4,4,2,2))
  rarecurve(t(seqtab[,ix_2015]), step = 1000, label = FALSE, xlab = 'Sequencing depth', ylab = 'Samples')
  abline(v = m, col = "red", lty = 2)
  
  dev.off()
  
  print(sort(rareslope(t(seqtab[,ix_2015]), m)))
  print(sort(rareslope(t(seqtab), m)))
}




