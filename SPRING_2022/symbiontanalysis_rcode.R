BiocManager::install("Rcolor")


library("phyloseq")
library("ggplot2")
library("readxl")
library("dplyr")
library("tibble")
library("microbiome")
library("ggpubr")

setwd("C:/Users/alexv/My Drive/Alex_Veglia/RICE/FLsctld/its2/2022-05-23_09-52-30.853221/its2_type_profiles")

#read in excel sheets with otu counts, otu taxonomy and sample info

otu_mat<- read_excel("flsctldgrant_phyloseq_abscounts.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("flsctldgrant_phyloseq_abscounts.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("flsctldgrant_phyloseq_abscounts.xlsx", sheet = "Samples")

#define row names for phyloseq
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu")

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample")

#transform otu and tax table into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#tranform to phyloseq obejcts
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
flsctldphy <- phyloseq(OTU, TAX, samples)
flsctldphy
#saving instacne to call phyloseq object later
saveRDS(flsctldphy, file="flsctld_grantsampsITS2.rds")
#call phyloseq object that was saved previously
flsctldphy <- readRDS("flsctld_grantsampsITS2.rds")

#####ABOVE CREATES/SAVES/LOADS THE PHYLOSEQ OBJECT
### creating named-color palette for types
psmelted <- psmelt(flsctldphy) #melt the giant phyloseq object into a dataframe to use with ggplot
fsdf <- psmelted
fsdf$Type ##check the list bruv
unique(fsdf$Type) ## checking the creation of a unique list of types
typeU <- unique(fsdf$Type)
typeU <- sort(typeU) #order the type alphabetically
#c('C34A2C', 'C04000', 'D4AF37', '16F529', 'C36241', '8B8000', 'C47451', '7E3517', 'FF0000', 'C11B17', '990000', '7E191B', '7D0552', 'DC381F', '872657', '66FF00', '2B65EC', 'A70D2A') ### hard code the colors for each type
types_palette = setNames(object = c('blue', 'indianred', 'red', 'red1', 'red2', 'red3', 'red4', 'tomato1', 'tomato2', 'tomato3', 'tomato4', 'yellow', 'yellow2', 'springgreen2', 'springgreen3', 'purple1', 'purple2'), nm = typeU)
print(types_palette)
### creating named-color palette for Genus
fsdf$Genus
unique(fsdf$Genus) ## create a unique list of types
genusU <- unique(fsdf$Genus)
#c('red', 'yellow', 'green', 'purple', 'blue') ### hard code the colors for each type
genus_palette = setNames(object = c('red', 'yellow', 'green', 'purple', 'blue'), nm = genusU)
print(genus_palette)

####CNAT PLOTS, B
#subsampling cnat samples that were diseased
cnat <- subset_samples(flsctldphy, sample_data(flsctldphy)$host_genus=="Colpophyllia")
cnatd <- subset_samples(cnat, sample_data(flsctldphy)$Health!="healthy")
cnatd <- prune_taxa(taxa_sums(cnatd)>0, cnat) #24 ITS2 profiles
#cnat <- prune_samples(sample_sums(cnat) >= 100, cnat) ## select samples with > 100 reads
cnatd  = transform_sample_counts(cnatd, function(x) x / sum(x) ) # only use for relative abundance
psmelted <- psmelt(cnatd) #melt the giant phyloseq object into a dataframe to use with ggplot
cnatd <- psmelted

#cnat relative abundance plot genus level for disease samples
n <- cnatd %>%
  ggplot(aes(x=Sample, y=Abundance, fill = Type)) +
  geom_bar(stat="identity", width = .7) +
  facet_grid(~host_genus, scales="free", space="free") +
  ylab("Relative Abundance") +  xlab("") +
  theme(axis.text.x= element_text(angle = 90, size=9, vjust=0.3),
        axis.text.y= element_text(size=12),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size=15),
        strip.background = element_blank()) +
  theme(strip.placement = "inside",
        strip.text = element_text(size=13),
        legend.text=element_text(size=11),
        panel.border=element_rect(fill=NA, colour="black",size=.7),
        panel.spacing = unit(1.0, "lines"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = types_palette)
n




flsctldphy <-subset_samples(flsctldphy, sample_data(flsctldphy)$TEM=="yes")
flsctldphy <- prune_samples(sample_sums(flsctldphy) >= 100, flsctldphy)
flsctldphy = transform_sample_counts(flsctldphy, function(x) x / sum(x) )









##convert phylo object to deseq object
diagdds = phyloseq_to_deseq2(cnat, ~sample)
diagdds = estimateSizeFactors(diagdds, type=c("poscounts"))


#boxplot
plot_bar(flsctldphy, fill= "Genus")
plot_bar(flsctldphy, x="sample", fill="Genus")
plot_bar(flsctldphy, "coral", fill="Genus", facet_grid=~host_genus)



###laurens code below here

library(ggplot2); packageVersion("ggplot2") #version‘3.3.5’
library(phyloseq); packageVersion("phyloseq") #version1.36.0
library(dplyr); packageVersion("dplyr") #version1.0.7
library(ggpubr); packageVersion("ggpubr") #version 0.4.0
library(DESeq2); packageVersion('DESeq2') #version 1.32.0
library(RColorBrewer)
library(vegan); packageVersion('vegan') #version 2.5.7
library(rstatix); packageVersion('rstatix') #version 0.7.0
library(patchwork); packageVersion('patchwork') #version 1.1.1
library(tidyr)


################################################################################
##Make phyloseq object from Symportal Genus-level absolute abundance output
################################################################################

setwd("/Users/laurenhowe-kerr/Dropbox/LaurenMVP/ACR/MVP-ACR-ITS2-SE-THESIS")
taxa <- as.matrix(read.csv("taxa_premed.csv", row.names = 1))
samdf <- read.csv("sample_info_healthdata_updated_csv.csv", row.names = 1)
counts <- read.csv("counts_uid_abs.csv", row.names=1, header = TRUE, check.names=FALSE)

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(counts, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))

ps

saveRDS(ps, file="ACR_MVP_symportal-genus.rds") #save as rds so that you can read in phyloseq object at any time

################################################################################
###################RELATIVE ABUNDANCE BAR PLOTS- GENUS##########################
################################################################################

ps <- readRDS("ACR_MVP_symportal-genus.rds")
ps <- subset_samples(ps, sample_data(ps)$Coral=="ACR")
ps <- subset_samples(ps, sample_data(ps)$Coral!="MESS") #remove MESS samples
ps <- subset_samples(ps, sample_data(ps)$Remove_Duplicate !="x")
ps <- prune_taxa(taxa_sums(ps)>0, ps) #24 ITS2 profiles

#to look at only samples with Aug 19 health data
CTs <- c("26", "27", "30", "31", "43", "44", "50", "145", "150", "157", "158", "159", "180", "188", "201", "204", "208", "209", "273", "432", "508", "514", "519")
ps <- subset_samples(ps, sample_data(ps)$CT %in% CTs)
ps <- prune_taxa(taxa_sums(ps)>0, ps) #24 ITS2 profiles

ps  = transform_sample_counts(ps, function(x) x / sum(x) ) # only use for relative abundance
psmelted <- psmelt(ps) #melt the giant phyloseq object into a dataframe to use with ggplot2

ps <- psmelted

ps$Time <- factor(ps$Time, levels= c("Sept.17", "Mar.18", "Aug.18", "Mar.19", "Aug.19", "Oct.20"))
