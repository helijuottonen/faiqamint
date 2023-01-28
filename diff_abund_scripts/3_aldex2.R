
DAaldex <- function(filt_otu, filt_meta, tr_name) {

library(ALDEx2)
library(cowplot)

filtdata.ax <- aldex.clr(filt_otu, filt_meta$Treatments, mc.samples=256, denom="all", verbose=F)
filtdata.ax.tt <- aldex.ttest(filtdata.ax, paired.test=FALSE, verbose=FALSE)
filtdata.ax.effect <- aldex.effect(filtdata.ax, CI=T, verbose=FALSE)
filtdata.ax.all <- data.frame(filtdata.ax.tt, filtdata.ax.effect)
write.csv2(filtdata.ax.all, paste("results/aldex2_", tr_name, ".csv", sep = ""))

#plotting PLOT EMPTY ?!
pdf(paste("plots/aldex2_", tr_name, ".pdf", sep = ""), width=9.5, height=4)
par(mfrow = c(1,2))
aldex.plot(filtdata.ax.all, type="MA", test="welch")
aldex.plot(filtdata.ax.all, type="MW", test="welch")
dev.off()

}

