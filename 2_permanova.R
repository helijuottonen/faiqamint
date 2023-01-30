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

# 2: PERMANOVA

# RUN FIRST:
# 0_basic_processing.R
# 1_NMDS.R

# library(vegan)

# remove start samples from metadata
w_no_start <- subset(wr, Treatments != "Start_sample" & Organ != "Root")
muc_no_start <- subset(muc, Treatments != "Start_sample")
gut_no_start <- subset(gut, Treatments != "Start_sample")

# same for gut sections separately
agut <- subset(gut_no_start, Organ == "Agut")
pgut <- subset(gut_no_start, Organ == "Pgut")


# choose correct rows from the otu table
w_otu_ns <- faotu.st[row.names(faotu.st) %in% row.names(w_no_start),]
muc_otu_ns <- faotu.st[row.names(faotu.st) %in% row.names(muc_no_start),]
gut_otu_ns <- faotu.st[row.names(faotu.st) %in% row.names(gut_no_start),]

agut_otu <- faotu.st[row.names(faotu.st) %in% row.names(agut),]
pgut_otu <- faotu.st[row.names(faotu.st) %in% row.names(pgut),]

# run PERMANOVA:
adonis2(w_otu_ns ~ Treatments, distance="bray", permutations=999, data=w_no_start)
adonis2(muc_otu_ns ~ Treatments, distance="bray", permutations=999, data=muc_no_start)
adonis2(gut_otu_ns ~ Treatments + Organ, distance="bray", permutations=999, data=gut_no_start, by="margin")

adonis2(agut_otu ~ Treatments, distance="bray", permutations=999, data=agut)
adonis2(pgut_otu ~ Treatments, distance="bray", permutations=999, data=pgut)




