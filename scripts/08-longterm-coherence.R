# Script to calculate and plot correlations (short and long-term)

library(ggplot2)
library(tidyverse)

## DATA ------------------------------------------------------------------------

# read in short-term correlations
sp_correlations = readRDS("outputs/gam_hierarchical_species_correlations.rds")
# clean up the species names 
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# read in the data to calculate long-term changes
dat = readRDS("outputs/data_train.rds")


## LONG-TERM SLOPES ------------------------------------------------------------

# fit linear model to each series
sp = unique(dat$series)
m <- list()
for(i in 1:length(sp)){
  m[[i]] = lm(y ~ time + 0, 
              data = filter(dat, series == sp[i]))
}

# extract slopes and intercepts
coefs = lapply(m, coef) |> 
  bind_rows() |>
  mutate("species" = sp)

# plot to just visualize what we're about to do here
dat_plot = dat
coefs_plot = coefs
colnames(coefs_plot)[1] = "slope"
dat_plot = left_join(dat_plot, coefs_plot, by = c("series" = "species"))
dat_plot |>
  ggplot(aes(
    x = time+1981,
    y = y,
    group = series,
    fill = slope
  )) +
  labs(x = "", y = "Biomass (kg per trawl)", fill = "Slope") +
  geom_smooth(method = "lm", alpha = .3, col = "black", linewidth = .3) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(legend.position = "right") 
ggsave("figures/longterm-coherence-lineartrends.png", width = 6, height = 4.5)

# make a pairwise dataframe to multiply the slopes
mat = cbind(expand_grid(coefs$species, coefs$species), expand_grid(coefs$time, coefs$time))
colnames(mat)[1:4] = c("sp1", "sp2", "slope1", "slope2")
mat$multiplied = mat$slope1*mat$slope2

# convert to a wide df
mat_w = mat |>
  select(c(sp1, sp2, multiplied)) |>
  pivot_wider(
    names_from = "sp1", 
    id_cols = "sp2", 
    values_from = "multiplied"
  ) |> as.data.frame()
rownames(mat_w) = mat_w$sp2
mat_w = mat_w[,-1]

# convert to matrix
mat_m = as.matrix(mat_w)
# clean up the species names
colnames(mat_m) = gsub("_", " ", colnames(mat_w)) |> stringr::str_to_sentence()
rownames(mat_m) = gsub("_", " ", rownames(mat_w)) |> stringr::str_to_sentence()

# remove same-species pairs
mat = mat |> filter(sp1 != sp2)
mat$sppair = NA
for(i in 1:nrow(mat)){
  mat$sppair[i] = sort(as.character(c(mat$sp1[i], mat$sp2[i]))) |> paste0(collapse = "_")
}
mat = distinct(mat, sppair,.keep_all = TRUE)

# plot species pairs as a grid
ggplot(data = mat) +
  geom_point(aes(x = slope1, y = slope2, fill = multiplied), 
             alpha = 1, size = 3, pch = 21) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  scale_fill_distiller(palette="RdBu", 
                       direction = 1, 
                       limits = c(-max(mat$multiplied), 
                                  max(mat$multiplied))) +
  ggpubr::theme_pubclean() +
  labs(x = "Slope of Species A", 
       y = "Slope of Species B",
       fill = "Correlation") +
  coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2))
ggsave("figures/longterm-coherence-gridplot.png", width = 7.66, height = 8)



# plot short-term as heatmap ---------------------------------------------------

png(height=1800, width=1800, file="figures/fig3b_species_associations_shortterm_latents.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", is.corr = FALSE,diag = FALSE,
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

# Plot as a histogram ----------------------------------------------------------

corr_df = data.frame("value" = sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))],
                     "group" = rep("A", length(sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))])))
(pC = ggplot(data = corr_df, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x))) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", x = "Correlation between species (R)", fill = "R") +
    geom_vline(xintercept = mean(corr_df$value)) +
    geom_vline(xintercept = mean(corr_df$value) - sd(corr_df$value), lty = 2) + 
    geom_vline(xintercept = mean(corr_df$value) + sd(corr_df$value), lty = 2) +
    theme(legend.position = "none")
)
ggsave("figures/fig3a_species_associations_shortterm_latent_histogram.png", width = 5, height = 2.5)

# plot long-term, B and D histogram --------------------------------------------

# plot as heatmap
png(height=1800, width=1800, file="figures/fig3b_species_associations_longterm_linearslope.png", type = "cairo")
corrplot::corrplot(mat_m, 
                   type = "lower",
                   method = "color", is.corr = FALSE,diag = FALSE,
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

# Plot as a histogram ---------------------------------------------------------

corr_df = data.frame("value" = mat_m[which(lower.tri(mat_m))],
                     "group" = rep("A", length(mat_m[which(lower.tri(mat_m))])))
(pD = ggplot(data = corr_df, 
                                 aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x)), col = "grey50", linewidth = .03) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-.04, .04)) +
    labs(y = "Frequency", x = "Correlation between species (R)", fill = "R") +
    geom_vline(xintercept = mean(corr_df$value)) +
    geom_vline(xintercept = mean(corr_df$value) - sd(corr_df$value), lty = 2) + 
    geom_vline(xintercept = mean(corr_df$value) + sd(corr_df$value), lty = 2) 
  + coord_cartesian(xlim = c(-.04, .04)) +
    theme(legend.position = "none")
)
ggsave("figures/fig3b_species_associations_longterm_linearslope_histogram.png", width = 5, height = 2.5)

# summary statistics -----------------------------------------------------------
corrs = mat_m[which(lower.tri(mat_m))]
mean(corrs)
sd(corrs)
corrs |> quantile(probs = c(0.05, 0.5, 0.95))
length(corrs[which(corrs > .75)]) / length(corrs)
length(corrs[which(corrs < -.75)]) / length(corrs)




## ALTERNATIVE - correlations between annual rates of change

# # read in the data
# derivs_pops_df = readRDS("outputs/derivs_pops_df.rds")
# 
# ggplot(data = derivs_pops_df) +
#   geom_line(aes(x = year, y = value, col = name))
# 
# # select just the slope information
# slopes = derivs_pops_df |>
#   select(c(year,name,value))
# 
# # pivot to wider
# slopes_m = slopes |>
#   pivot_wider(names_from = name, values_from = value) |>
#   slice(-1) # remove baseline 0 year
# rownames(slopes_m) = slopes_m$year
# slopes_m = select(slopes_m, -c(year))
# 
# # species pairs
# pairs = expand.grid(colnames(slopes_m), colnames(slopes_m))
# pairs$cor = NA
# for(i in 1:nrow(pairs)){
#   pairs$cor[i] = cor(select(slopes_m, pairs[i,1]),
#                      select(slopes_m, pairs[i,2]))
# }
# pairs_nodupes = pairs[-which(pairs$Var1 == pairs$Var2),]
# hist(pairs_nodupes$cor)
# 
# # make a matrix
# pairs_w = pivot_wider(
#   pairs,
#   names_from = Var1, values_from = cor, id_cols = Var2
# )
# 
# # clean up the species names 
# colnames(pairs_w) = gsub("_", " ", colnames(pairs_w)) |> 
#   stringr::str_to_sentence()
# rownames(pairs_w) = gsub("_", " ", pairs_w$Var2) |> 
#   stringr::str_to_sentence()
# pairs_w = select(pairs_w, -c(Var2))
# pairs_m = as.matrix(pairs_w)
# rownames(pairs_m) = pairs_w$Var2
# # plot as heatmap
# png(height=1800, width=1800, file="figures/fig3b_species_associations_longterm.png", type = "cairo")
# corrplot::corrplot(pairs_m, 
#                    type = "lower",
#                    method = "color", 
#                    tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
# dev.off()
# 
# # Plot as a histogram ----------------------------------------------------------
# 
# corr_df = data.frame("value" = pairs_m[which(lower.tri(pairs_m))],
#                      "group" = rep("A", length(pairs_m[which(lower.tri(pairs_m))])))
# (correlations_histogram = ggplot(data = corr_df, aes(x = value)) +
#     geom_histogram(aes(fill = after_stat(x))) +
#     scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
#     labs(y = "Frequency", x = "Correlation between species (R)", fill = "R") +
#     geom_vline(xintercept = mean(corr_df$value)) +
#     geom_vline(xintercept = mean(corr_df$value) - sd(corr_df$value), lty = 2) + 
#     geom_vline(xintercept = mean(corr_df$value) + sd(corr_df$value), lty = 2) +
#     ylim(c(0,45))
# )
# ggsave("figures/fig3b_species_associations_longterm_histogram.png", width = 6.4, height = 7)
# 
# 
# # take correlations
# cor(slopes_m$`Anarhichas denticulatus`, slopes_m$`Anarhichas lupus`)
