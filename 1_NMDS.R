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

# 1: Running NMDS

# RUN FIRST: 0_basic_processing.R

# packages needed:
# library(vegan) # loaded in 0_basic_processing.R
library(ggplot2)
library(gridExtra)

# separating fish samples vs. water & root

wr <- subset(fameta, Organ == "Water" | Organ == "Root")
muc <- subset(fameta, Organ == "Mucous")
gut <- subset(fameta, Organ == "Agut" | Organ == "Pgut")


# defining a function to run NMDS:

subset_nmds_ms2 <- function(x) {

# choose correct rows from the otu table
otu <- faotu.st[row.names(faotu.st) %in% row.names(x),]

set.seed(1732)

# run NMDS on the standardized OTU table with the vegan package function metaMDS and Bray-Curtis distances
mds <- metaMDS(otu, distance = "bray", trymax=200, trace=FALSE)

# extracting the sample scores from the NMDS result for plotting
mds.sc <- scores(mds, display="sites", choices=c(1,2))

# making a dataframe of the sample scores
mds.sc.df <- as.data.frame(mds.sc)

  print(mds$stress)
  return(mds.sc.df)
}

# running NMDS for water+roots, mucous, gut
wr.mds.sc.df <- subset_nmds_ms2(wr)
# stress 0.1013588
muc.mds.sc.df <- subset_nmds_ms2(muc)
# stress 0.143949
gut.mds.sc.df <- subset_nmds_ms2(gut)
# stress 0.1919543

# plotting the NMDS with ggplot2

# defining factors for plotting
organ_wr = factor(wr$Organ,c("Water", "Root"))
organ_muc = factor(muc$Organ,c("Mucous"))
organ_gut = factor(gut$Organ,c("Agut", "Pgut"))
treat_wr = factor(wr$Treatments,c("Start_sample", "Without_mint", "With_mint"))
treat_muc = factor(muc$Treatments,c("Start_sample", "Without_mint", "With_mint"))
treat_gut = factor(gut$Treatments,c("Start_sample", "Without_mint", "With_mint"))

# defining colours for plotting
wr_pal <- c("#7BAFDE", "#90C087") #"#1965B0"
muc_pal <- c("#F1932D")
gut_pal <- c("#D1BBD7", "#882E72")

# making the plot
faotuplot_wr <- ggplot(wr.mds.sc.df, aes(NMDS1, NMDS2, color=organ_wr, shape=treat_wr)) + geom_point(size=3, stroke=1)
faotuplot_muc <- ggplot(muc.mds.sc.df, aes(NMDS1, NMDS2, color=organ_muc, shape=treat_muc)) + geom_point(size=3, stroke=1)
faotuplot_gut <- ggplot(gut.mds.sc.df, aes(NMDS1, NMDS2, color=organ_gut, shape=treat_gut)) + geom_point(size=3, stroke=1)


# modifying the plot further (adding a theme, colours, symbols)
# water and roots:
faotuplot_wr2 <- faotuplot_wr + theme_classic() + 
  scale_colour_manual(values=wr_pal, name="Sample type") + 
  scale_shape_manual(values=c(3, 1, 17), name="Treatment", labels = c("Start sample", "RAS", "Aquaponics")) +
  ggtitle("Water and mint root") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=11),
        plot.title = element_text(size=12, face="bold")) + 
  guides(colour = guide_legend(order = 1), 
        shape = guide_legend(order = 2))

# mucous: 
faotuplot_muc2 <- faotuplot_muc + theme_classic() + 
  scale_colour_manual(values=muc_pal, name="Sample type") + 
  scale_shape_manual(values=c(3, 1, 17), name="Treatment", labels = c("Start sample", "RAS", "Aquaponics")) +
  ggtitle("Fish mucous") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=11),
        plot.title = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

# gut:
faotuplot_gut2 <- faotuplot_gut + theme_classic() + 
  scale_colour_manual(values=gut_pal, name="Sample type", labels = c("Anterior gut", "Posterior gut")) + 
  scale_shape_manual(values=c(3, 1, 17), name="Treatment", labels = c("Start sample", "RAS", "Aquaponics")) +
  ggtitle("Fish gut") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=11),
        plot.title = element_text(size=12, face="bold")) +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

# checking how the plots look
faotuplot_wr2
faotuplot_muc2
faotuplot_gut2

# defining a layout for a combined plot
lay = rbind(c(1,1,1,2,2,2),
            c(3,3,3,3,3,NA))

# using gridExtra to combine the plots
nmds_all <- grid.arrange(arrangeGrob(faotuplot_wr2), arrangeGrob(faotuplot_muc2), arrangeGrob(faotuplot_gut2), layout_matrix = lay)

#save the plot as a pdf file in folder 'plots' on the R project folder
ggsave(file = "plots/ms2_nmds_2023.pdf", nmds_all, width = 8.5, height = 6) 

