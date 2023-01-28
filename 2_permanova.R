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

# PERMANOVA results: 
#### water
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = w_otu_ns ~ Treatments, data = w_no_start, permutations = 999, distance = "bray")
# Df          SumOfSqs     R2      F Pr(>F)  
# Treatments  1   0.5480 0.1624 1.9389  0.081 .
# Residual   10   2.8263 0.8376                
# Total      11   3.3743 1.0000                  

#### mucous
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = muc_otu_ns ~ Treatments, data = muc_no_start, permutations = 999, distance = "bray")
#             Df SumOfSqs      R2      F Pr(>F)    
# Treatments  1   0.6464 0.17131 3.3075  0.001 ***
# Residual   16   3.1268 0.82869                  
# Total      17   3.7731 1.00000     

#### gut
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = gut_otu_ns ~ Treatments + Organ, data = gut_no_start, permutations = 999, by = "margin", distance = "bray")
#             Df SumOfSqs      R2      F Pr(>F)   
# Treatments  1   0.7901 0.08125 3.0772  0.005 **
# Organ       1   0.4606 0.04737 1.7940  0.098 . 
# Residual   33   8.4728 0.87137                 
# Total      35   9.7235 1.00000               

### agut
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = agut_otu ~ Treatments, data = agut, permutations = 999, distance = "bray")
#           Df SumOfSqs      R2      F Pr(>F)  
# Treatments  1   0.6874 0.14872 2.7952  0.036 *
# Residual   16   3.9348 0.85128                
# Total      17   4.6222 1.00000   

### pgut
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = pgut_otu ~ Treatments, data = pgut, permutations = 999, distance = "bray")
#           Df SumOfSqs      R2      F Pr(>F)
# Treatments  1   0.4069 0.08768 1.5377  0.153
# Residual   16   4.2338 0.91232              
# Total      17   4.6407 1.00000 


