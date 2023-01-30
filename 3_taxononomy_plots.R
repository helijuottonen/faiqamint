# Heli Juottonen 2023 heli.juottonen@alumni.helsinki.fi
# R scripts for bacterial 16S rRNA gene sequenec data analysis of Faiqa Atique's PhD manuscript on aquaponics with mint
# https://github.com/helijuottonen/faiqamint

# For the article: xxx
# Atique F 12, Juottonen H 1, Kytöviita MM 1
# 1 Department of Biological and Environmental Sciences, University of Jyväskylä, Finland
# 2 Institute of Bioeconomy, JAMK University of Applied Sciences, Finland

# Disclaimer: I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices.
# Use at your own risk!

# 2: Taxonomy plots

# RUN FIRST: 
# 0_basic_processing.R

# library(vegan) #loaded in 0_basic_processing.R
# library(tidyverse) #loaded in 0_basic_processing.R
library(phyloseq) 
library(forcats)
library(cowplot) 
library(khroma)

# To install phyloseq (R v.4.2 or newer), run these lines
# (from https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html):

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
# BiocManager::install("phyloseq")

################ creating a phyloseq object

# defining taxonomy file
taxfile = "data/fams2_all.unique.good.filter.unique.precluster.denovo.uchime.pick.opti_mcc.0.03.cons.taxonomy"

#importing the taxonomy table into phyloseq format
mothur_tax <- import_mothur(mothur_constaxonomy_file = taxfile)

# defining OTU table object
otutable <- otu_table(as.matrix(faotu.r.all), taxa_are_rows=FALSE)
# defining metadata object
metadata <- sample_data(fameta)

# creating the phyloseq object
faps <- merge_phyloseq(otutable, metadata, mothur_tax)

# giving proper names to the columns of the taxonomy tables
colnames(tax_table(faps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# check how the phyloseq object looks 
faps

############ taxonomy plot at order level (change the taxonomic level if needed)

# removing unwanted sample types and splitting data to groups

# water + mucous
faps_wm <- subset_samples(faps, Organ == "Mucous" | Organ == "Water")
# gut
faps_gut <- subset_samples(faps, Organ == "Agut" | Organ == "Pgut")
# mint root
faps_root <- subset_samples(faps, Organ == "Root")


# merging OTUs by order

# water + mucous
faps_wm_order <- faps_wm %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)  

# gut
faps_gut_order <- faps_gut %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)

# mint root
faps_root_order <- faps_root %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)

### Combining orders with relative abundance <1% into 'other'

# WATER + MUCOUS
# getting rid of extra metadata columns because they seem to mess up pivot_wider &
# converting to wide for selecting orders with rel.abund. >1%
faps_wm_order_wide <- faps_wm_order %>% dplyr::select(OTU, Sample, Abundance, name, Treatments, Organ, Order) %>% 
  pivot_wider(names_from = Order, values_from = Abundance)

wm_bcolsums <- colSums(faps_wm_order_wide[6:253], na.rm=TRUE)
wm_brel <- as.data.frame(wm_bcolsums/sum(wm_bcolsums))
names(wm_brel)[1] <- "rel.abundance"
# selecting orders with overall relative abundance over 1%
wm_brel_sub <- subset(wm_brel, rel.abundance > 0.01)
# orders > 1% as a vector
wm_ord_over1 <- row.names(wm_brel_sub)

wm_faps_order_other <- faps_wm_order %>%
  # renaming orders <1% as 'other'
  mutate(Order_mod = ifelse(Order %in% wm_ord_over1, Order, "other")) %>%
  # converting read counts to relative abundances
  mutate(Rel.Abundance = Abundance/sum(Abundance)) %>%
  # removing orders/rows with no reads
  filter(Abundance !=0)

# check the number of orders left after combining orders <1% relative abundance to 'other'
# the aim is to make sure there are not too many for selecting colours for each order
length(unique(wm_faps_order_other[,"Order_mod"]))

# GUT
# getting rid of extra metadata columns because they seem to mess up pivot_wider &
# converting to wide for selecting orders with rel.abund. >1%
faps_gut_order_wide <- faps_gut_order %>% dplyr::select(OTU, Sample, Abundance, name, Treatments, Organ, Order) %>% 
  pivot_wider(names_from = Order, values_from = Abundance)

gut_bcolsums <- colSums(faps_gut_order_wide[6:253], na.rm=TRUE)
gut_brel <- as.data.frame(gut_bcolsums/sum(gut_bcolsums))
names(gut_brel)[1] <- "rel.abundance"
# selecting orders with overall relative abundance over 1%
gut_brel_sub <- subset(gut_brel, rel.abundance > 0.01)
# orders > 1% as a vector
gut_ord_over1 <- row.names(gut_brel_sub)

gut_faps_order_other <- faps_gut_order %>%
  # renaming orders <1% as 'other'
  mutate(Order_mod = ifelse(Order %in% gut_ord_over1, Order, "other")) %>%
  # converting read counts to relative abundances
  mutate(Rel.Abundance = Abundance/sum(Abundance)) %>%
  # removing orders/rows with no reads
  filter(Abundance !=0)

# renaming and reordering for plotting
gut_faps_order_other$Organ <- fct_recode(gut_faps_order_other$Organ, Anterior_gut = "Agut", Posterior_gut = "Pgut")

# checking the number of colours needed
length(unique(gut_faps_order_other[,"Order_mod"]))

# MINT ROOT
faps_root_order_wide <- faps_root_order %>% dplyr::select(OTU, Sample, Abundance, name, Treatments, Organ, Order) %>% 
  pivot_wider(names_from = Order, values_from = Abundance)

root_bcolsums <- colSums(faps_root_order_wide[6:253], na.rm=TRUE)
root_brel <- as.data.frame(root_bcolsums/sum(root_bcolsums))
names(root_brel)[1] <- "rel.abundance"
root_brel_sub <- subset(root_brel, rel.abundance > 0.01)
root_ord_over1 <- row.names(root_brel_sub)

root_faps_order_other <- faps_root_order %>%
  mutate(Order_mod = ifelse(Order %in% root_ord_over1, Order, "other")) %>%
  mutate(Rel.Abundance = Abundance/sum(Abundance)) %>%
  filter(Abundance !=0)

root_faps_order_other$Treatments <- fct_recode(root_faps_order_other$Treatments, Start_sample = "Start_sample", End_sample = "aquaponics")

length(unique(root_faps_order_other[,"Order_mod"]))

##### PLOTTING

# defining colours for plotting
discrete_rainbow <- colour("discrete rainbow")
bpal <- as.vector(discrete_rainbow())  # from: khroma

# water+mucous
wm_order_plot <- ggplot(wm_faps_order_other, aes(x = name, y = Rel.Abundance, fill = Order_mod, color=Order_mod)) + 
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = bpal, limits=force, name="Order") +
  scale_colour_manual(values = bpal, limits=force, guide="none") +
  ylab("Relative abundance") +
  theme_cowplot(10) + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(. ~ Organ + Treatments, scales="free_x", space="free")

wm_order_plot

# gut
gut_order_plot <- ggplot(gut_faps_order_other, aes(x = name, y = Rel.Abundance, fill = Order_mod, color=Order_mod)) + 
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = bpal, limits=force, name="Order") +
  scale_colour_manual(values = bpal, limits=force, guide="none") +
  ylab("Relative abundance") +
  theme_cowplot(10) + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(. ~ Organ + Treatments, scales="free_x", space="free")

gut_order_plot

# root
root_order_plot <- ggplot(root_faps_order_other, aes(x = name, y = Rel.Abundance, fill = Order_mod, color=Order_mod)) + 
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = bpal, limits=force, name="Order") +
  scale_colour_manual(values = bpal, limits=force, guide="none") +
  ylab("Relative abundance") +
  theme_cowplot(10) + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(~Treatments, scales="free_x", space="free")

root_order_plot

#save the plotS as a pdf file in folder 'plots' on the R project folder
ggsave("plots/wm_order_plot.pdf", wm_order_plot, width=12, height=5) 
ggsave("plots/gut_order_plot.pdf", gut_order_plot, width=12, height=5) 
ggsave("plots/root_order_plot.pdf",root_order_plot, width=8, height=5) 


