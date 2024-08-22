# Script to make figures of the correlations between population growth rates:
## correlation matrix plot
## histogram of the correlation btwn all population growth rates over all years

# libraries ----

library(here)
library(dplyr)
library(tidyr)
library(mvgam)
library(tidybayes)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggdist)
library(gganimate)
theme_set(theme_pubr() +
            theme(panel.grid.major.x = element_line()))

set.seed(12)

# read the previous data objects and the model object --------------------------

biomass_scaling = readRDS("outputs/biomass_scaling.rds")
mod1 = readRDS(paste0("outputs/gam_hierarchical_gp.rds"))
data_train = mod1$obs_data

# calculate species correlations -----------------------------------------------
sp_correlations = lv_correlations(mod1)
saveRDS(sp_correlations, "outputs/gam_hierarchical_species_correlations.rds")

# clean up the species names 

colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# Plot as heatmap --------------------------------------------------------------

png(height=1800, width=1800, file="figures/fig3_species_associations.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

# Plot as a histogram ----------------------------------------------------------

corr_df = data.frame("value" = sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))],
                     "group" = rep("A", length(sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))])))
(correlations_histogram = ggplot(data = corr_df, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x))) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", x = "Correlation between species (R)", fill = "R") +
    geom_vline(xintercept = mean(corr_df$value)) +
    geom_vline(xintercept = mean(corr_df$value) - sd(corr_df$value), lty = 2) + 
    geom_vline(xintercept = mean(corr_df$value) + sd(corr_df$value), lty = 2)
)

# summary statistics -----------------------------------------------------------
corrs = sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))]
corrs |> quantile(probs = c(0.05, 0.5, 0.95))
length(corrs[which(corrs > .75)]) / length(corrs)
length(corrs[which(corrs < -.75)]) / length(corrs)
