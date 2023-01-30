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

# 6: Alpha diversity

# RUN FIRST:
# 0_basic_processing.R

library(readr)
library(ggplot2)
library(gridExtra)

# converting otu table into dataframe
faotu_div.df <- as.data.frame(faotu)

# removing shoot and control samples
faotu_div.df <- faotu_div.df[faotu_div.df$Group %in% row.names(fameta),]

min(rowSums(faotu_div.df[,4:ncol(faotu_div.df)]))
# 256

# removing start samples
fameta_div <- subset(fameta, Treatments != "Start_sample")
faotu_div.df <- faotu_div.df[faotu_div.df$Group %in% row.names(fameta_div),]

min(rowSums(faotu_div.df[,4:ncol(faotu_div.df)]))
# 1332

# saving otu table as shared file for diversity calculations in mothur
write_delim(file = "data/faotu_div.shared", faotu_div.df)

# calculating diversity indices in mothur
# in terminal (navigate to the mothur folder first):

# ./mothur "#summary.single(shared=/Users/helij/Documents/Työt/Faiqa/ms2/R_ms2/data/faotu_div.shared, calc=nseqs-sobs-coverage-shannon, subsample=T)"
 
fa_div <- read.table("data/faotu_div.groups.ave-std.summary", header = T)

fa_div2 <- subset(fa_div, method == "ave")
row.names(fa_div2) <- fa_div2$group
fa_div2 <- fa_div2[order(row.names(fa_div2)),]

div_meta <- merge(fa_div2, fameta_div, by = "row.names")
div_meta <- subset(div_meta, Organ != "Root")

# setting  labels
div_meta$Organ <- fct_recode(div_meta$Organ, Mucous = "Mucous", Water = "Water", Anterior_gut = "Agut", Posterior_gut = "Pgut")
div_meta$Treatments <- fct_recode(div_meta$Treatments, Aquaponics = "aquaponics", RAS = "RAS")
div_meta$Organ <- fct_relevel(div_meta$Organ, "Mucous", "Water", "Anterior_gut", "Posterior_gut")

shplot <- ggplot(div_meta, aes(Treatments, shannon, fill = Treatments)) +
  geom_boxplot() +
  facet_wrap(~ Organ) + 
  ylab("Shannon") +
  theme_cowplot(11) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_blank())
  
shplot

sobsplot <- ggplot(div_meta, aes(Treatments, sobs, fill = Treatments)) +
  geom_boxplot() +
  facet_wrap(~ Organ) +
  ylab("Number of OTUs") +
  theme_cowplot(11) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_blank())

sobsplot

sobs_shannon <- grid.arrange(sobsplot, shplot)

ggsave(file = "plots/sobs_shannon.pdf", sobs_shannon, width = 6, height = 8)


t.test(shannon ~ Treatments, data=div_meta)

# Welch Two Sample t-test
# 
# data:  shannon by Treatments
# t = -1.9397, df = 62.311, p-value = 0.05694
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.94057752  0.01409788
# sample estimates:
#   mean in group Aquaponics        mean in group RAS 
# 2.996475                 3.459715 
