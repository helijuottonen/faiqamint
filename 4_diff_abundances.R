# Heli Juottonen 2023 heli.juottonen@alumni.helsinki.fi
# R scripts for bacterial 16S rRNA gene sequence data analysis of Faiqa Atique's PhD manuscript on aquaponics with mint
# https://github.com/helijuottonen/faiqamint

# For the article: xxx
# Atique F 12, Juottonen H 1, Kytöviita MM 1
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Institute of Bioeconomy, JAMK University of Applied Sciences, Finland

# Disclaimer: I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices.
# Use at your own risk!

# 4: Differential abundance analysis

# RUN FIRST:  
# 0_basic_processing.R (gives fameta, faotu.df2)
# 3_taxonomy_plots.R (gives phyloseq object faps)

# library(phyloseq)
# library(tidyverse)
# library(cowplot)
library(readr)

##### filtering OTU data

# function for selecting and splitting data by fish organ

da_filtering <- function(x) {

  # subsetting metadata based on on fish organ & end samples
  meta_subset <- subset(fameta, Organ == x & Treatments != "Start_sample")
                          
  meta_subset <- droplevels(meta_subset)

  # selecting the subsetted samples from the OTU table
  otu_subset = faotu.df2[row.names(faotu.df2) %in% row.names(meta_subset),]

  # selecting only OTUs with >0 reads in at least 10% of samples
  otu_subset10 <- otu_subset[,colSums(otu_subset > 0) > (0.1*(dim(otu_subset)[1]))]

  # transposing OTU data frame (all the DA methods want OTUs in rows, samples in columns)
  otu_subset10t <- as.data.frame(t(otu_subset10))
  
  x_otu_filt <- otu_subset10t  
  x_meta_filt <- meta_subset
  x_meta_filt$sampleID <- rownames(x_meta_filt) # important!
  
  return(list(x_otu_filt, x_meta_filt))

}

muc <- da_filtering(x = "Mucous")
pgut <- da_filtering(x = "Pgut")
agut <- da_filtering(x = "Agut")

# access dataframes in lists: dataframe[[1]] = otu table, dataframe[[2]] = metadata

##### differential abundance analyses

# ANCOM (not ANCOM-BC because there sample size 10 or more recommended)
source("diff_abund_scripts/1_ancom.R")

DAancom(muc[[2]], muc[[1]], "sampleID", "Treatments", "muc")
DAancom(pgut[[2]], pgut[[1]], "sampleID", "Treatments", "pgut")
DAancom(agut[[2]], agut[[1]], "sampleID", "Treatments", "agut")

# DESeq2
source("diff_abund_scripts/2_deseq2.R")

DAdeseq(muc[[1]], muc[[2]], "muc")
DAdeseq(pgut[[1]], pgut[[2]], "pgut")
DAdeseq(agut[[1]], agut[[2]], "agut")

# ALDEx2

source("diff_abund_scripts/3_aldex2.R")

DAaldex(muc[[1]], muc[[2]], "muc")
DAaldex(pgut[[1]], pgut[[2]], "pgut")
DAaldex(agut[[1]], agut[[2]], "agut")

##### filtering and collecting the results

# reading in results and setting 1 = detected as differentially abundance, 0 = not
# criteria different for each DA method

da_res <- function(treat) {

# ancom
treat_ac <- read_csv2(file = paste("results/ancom_", treat, ".csv", sep = "")) %>% 
  dplyr::select(taxa_id, detected_0.7) %>% 
  mutate(ac_sel = case_when(detected_0.7 == TRUE ~ 1, detected_0.7 == FALSE ~ 0))
 
# deseq2
treat_ds <- read_csv2(file = paste("results/deseq_", treat, ".csv", sep = "")) %>% 
  dplyr::rename(taxa_id = 1) %>% 
  dplyr::select(taxa_id, log2FoldChange, padj) %>% 
  mutate(ds_sel = case_when(padj < 0.05 & abs(log2FoldChange) > 1 ~ 1, 
      padj > 0.05 & abs(log2FoldChange) < 1 ~ 0 ))

# aldex2
treat_ax <- read_csv2(file = paste("results/aldex2_", treat, ".csv", sep = "")) %>% 
  dplyr::rename(taxa_id = 1) %>% 
  dplyr::select(taxa_id, effect) %>% 
  mutate(ax_sel = case_when(effect > 1 ~ 1, effect < 1 ~ 0))

  return(list(treat_ac, treat_ds, treat_ax))

}

muc_da <- da_res("muc")
pgut_da <- da_res("pgut")
agut_da <- da_res("agut")

# 1 = ancom, 2 = deseq2, 3 = ancom

# combining results of different DA methods

muc_comb <- as.data.frame(cbind(taxa_id = muc_da[[1]]$taxa_id, ac_sel = muc_da[[1]]$ac_sel, ds_sel = muc_da[[2]]$ds_sel, ax_sel = muc_da[[3]]$ax_sel))
pgut_comb <- as.data.frame(cbind(taxa_id = pgut_da[[1]]$taxa_id, ac_sel = pgut_da[[1]]$ac_sel, ds_sel = pgut_da[[2]]$ds_sel, ax_sel = pgut_da[[3]]$ax_sel))
agut_comb <- as.data.frame(cbind(taxa_id = agut_da[[1]]$taxa_id, ac_sel = agut_da[[1]]$ac_sel, ds_sel = agut_da[[2]]$ds_sel, ax_sel = agut_da[[3]]$ax_sel))

# selecting OTUs detected as differentially abundant by 2 or more DA methods
  # NAs replaced by 0

muc_comb_sel <- muc_comb %>% replace(is.na(.), 0) %>% 
  mutate(da_sum = as.numeric(ac_sel) + as.numeric(ds_sel) + as.numeric(ax_sel)) %>%   
  filter(da_sum > 1)

pgut_comb_sel <- pgut_comb %>% replace(is.na(.), 0) %>% 
  mutate(da_sum = as.numeric(ac_sel) + as.numeric(ds_sel) + as.numeric(ax_sel)) %>%     
  filter(da_sum > 1)

agut_comb_sel <- agut_comb %>% replace(is.na(.), 0) %>% 
  mutate(da_sum = as.numeric(ac_sel) + as.numeric(ds_sel) + as.numeric(ax_sel)) %>%     
  filter(da_sum > 1)

# adding taxonomy to the list of differentially abundance OTUs

# faps = phyloseq object from 2_taxonomy_plots.R
fps.tax = as(tax_table(faps), "matrix")
fps.tax.df = as.data.frame(fps.tax)

muc_tax <- fps.tax.df[row.names(fps.tax.df) %in% muc_comb_sel$taxa_id, ]
pgut_tax <- fps.tax.df[row.names(fps.tax.df) %in% pgut_comb_sel$taxa_id, ]
agut_tax <- fps.tax.df[row.names(fps.tax.df) %in% agut_comb_sel$taxa_id, ]

write.csv2(file = "results/muc_tax.csv", muc_tax)
write.csv2(file = "results/pgut_tax.csv", pgut_tax)
write.csv2(file = "results/agut_tax.csv", agut_tax)

muc_tax <- read.csv2(file = "results/muc_tax.csv", row.names = 1)
pgut_tax <- read.csv2(file = "results/pgut_tax.csv", row.names = 1)
agut_tax <- read.csv2(file = "results/agut_tax.csv", row.names = 1)


muc_ds <- as.data.frame(muc_da[[2]]) 
muc_ds_sel <- muc_ds[muc_ds$taxa_id %in% row.names(muc_tax), ]
muc_tax$log2FoldChange <- muc_ds_sel$log2FoldChange
muc_tax$otugenus <- paste(muc_tax$Genus, row.names(muc_tax), sep=" ")
muc_tax$taxa_id <- row.names(muc_tax)

pgut_ds <- as.data.frame(pgut_da[[2]]) 
pgut_ds_sel <- pgut_ds[pgut_ds$taxa_id %in% row.names(pgut_tax), ]
pgut_tax$log2FoldChange <- pgut_ds_sel$log2FoldChange
pgut_tax$otugenus <- paste(pgut_tax$Genus, row.names(pgut_tax), sep=" ")
pgut_tax$taxa_id <- row.names(pgut_tax)

muc_tax$otugenus <- gsub("_unclassified", " uncl.", muc_tax$otugenus)
muc_tax$otugenus <- gsub("_ge", "", muc_tax$otugenus)
muc_tax$otugenus <- gsub("Candidatus_", "Ca. ", muc_tax$otugenus)
muc_tax$otugenus <- gsub("_marine_group", " marine group", muc_tax$otugenus)
muc_tax$Class <- gsub("Bacteria_unclassified", "Unclassified Bacteria", muc_tax$Class)

muc_tax$otugenus = with(muc_tax, reorder(otugenus, log2FoldChange, median))

discrete_rainbow <- colour("discrete rainbow")
bpal7 <- as.vector(discrete_rainbow(13))  # from: khroma

muc_tax.plot <- ggplot(muc_tax, aes(otugenus, log2FoldChange, color=Class)) +
  geom_hline(yintercept = 0) +
  geom_point(size=4) + 
  theme_cowplot(12) +
  theme(panel.grid.major = element_line(colour = "grey95")) +
  scale_colour_manual(values = bpal7) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

muc_tax.plot

ggsave("plots/muc_da.pdf", muc_tax.plot, width=7, height=6) 
