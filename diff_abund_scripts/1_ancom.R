
# https://github.com/FrederickHuangLin/ANCOM

# not ANCOM-BC because samples size per group >10 (5 not recommended)

# run first: 0_DA_fitlering.R

# ANCOM preprocessing

# data = list of dataframes resulting from DA filtering in 4_diff_abundances.R
# svar = variable with sample names
# gvar = treatment / factor being tested

DAancom <- function(filt_meta, filt_otu, svar, gvar, tr_name) {

source("diff_abund_scripts/ancom.R")

meta_data = filt_meta
feature_table = filt_otu; sample_var = svar; group_var = gvar
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE

prepro = feature_table_pre_process(feature_table, meta_data, sample_var = svar, group_var, 
                                   out_cut=0.05, zero_cut=0.90, lib_cut=1000, neg_lb=FALSE)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = gvar; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL; lme_control = NULL

res_data = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula, lme_control)

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa - 1), 0.8 * (n_taxa - 1), 0.7 * (n_taxa - 1), 0.6 * (n_taxa - 1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

cut_off
#detected_0.9 detected_0.8 detected_0.7 detected_0.6 

# Annotation data
dat_ann = data.frame(x = min(res_data$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res_data$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig
pdf(paste("plots/ancom_",tr_name,".pdf", sep = ""))
print(fig)
dev.off()

write_csv2(res_data$out, paste("results/ancom_",tr_name,".csv", sep = ""))

}









