# Script to take the derivatives of each species trend (to get their annual growth rates)

# libraries --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(mvgam)
library(tidybayes)
library(ggplot2)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr() +
            theme(panel.grid.major.x = element_line()))

set.seed(12)

# load modified function that extracts derivatives from a plot function
source("scripts/plot_mvgam_trend_custom.R")

# read the previous data objects and the model object ----

biomass_scaling = readRDS("outputs/biomass_scaling.rds")
mod1 = readRDS(paste0("outputs/gam_hierarchical_gp.rds"))
npops = biomass_scaling$npops
time = biomass_scaling$time

# Take the derivative of each population trend ---------------------------------

trend_vals = list()
derivs_ls = list()
preds_ls = list()
for(i in 1:npops){
  trend_vals[[i]] = plot_mvgam_trend_custom(mod1, derivatives = TRUE, series = i)
  preds_ls[[i]] = trend_vals[[i]]$preds
  derivs_ls[[i]] = trend_vals[[i]]$derivs
}
preds = do.call(rbind, preds_ls)
derivs = do.call(rbind, derivs_ls)
saveRDS(preds, "outputs/predictions.rds")
saveRDS(derivs, "outputs/derivatives.rds")

# Take each species' median predicted trend in biomass
preds_pops = do.call(cbind, lapply(preds_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(preds_pops, x = time+1981, type = "l")

# each species' median derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE)))
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))
matplot(derivs_pops, x = time+1981, type = "l")

# Calculate each species' mean derivative over the whole time series -----------

(temporal_trend = data.frame(
  species = colnames(biomass_scaling$biomass),
  mu_deriv = derivs_pops |> apply(2, median, na.rm = TRUE),
  lower = derivs_pops |> apply(2, quantile, probs = .05, na.rm = TRUE),
  upper = derivs_pops |> apply(2, quantile, probs = .95, na.rm = TRUE),
  sd = derivs_pops |> apply(2, sd, na.rm = TRUE),
  mu = derivs_pops |> apply(2, mean, na.rm = TRUE)
))
temporal_trend$species = gsub("_", " ", temporal_trend$species) |> stringr::str_to_sentence()
temporal_trend$species =  factor(temporal_trend$species,
                                 levels = temporal_trend$species[order(temporal_trend$mu_deriv)])
saveRDS(temporal_trend, "outputs/temporal_trend.rds")
temporal_trend = readRDS("outputs/temporal_trend.rds")

#### Average and variance of growth rates of the whole assemblage per year -----

(temporal_trend_yearly = data.frame(
  year = rownames(Year_Geom_Means_all),
  mu_deriv = derivs_pops |> apply(1, median, na.rm = TRUE),
  lower = derivs_pops |> apply(1, quantile, probs = .05, na.rm = TRUE),
  upper = derivs_pops |> apply(1, quantile, probs = .95, na.rm = TRUE),
  sd = derivs_pops |> apply(1, sd, na.rm = TRUE),
  mu = derivs_pops |> apply(1, mean, na.rm = TRUE)
))
temporal_trend_yearly$year = as.character(temporal_trend_yearly$year)
saveRDS(temporal_trend_yearly, "outputs/temporal_trend_yearly.rds")

# Convert the derivatives to long format

colnames(derivs_pops) = colnames(biomass_scaling$biomass)
colnames(derivs_pops_lower) = colnames(biomass_scaling$biomass)
colnames(derivs_pops_upper) = colnames(biomass_scaling$biomass)
derivs_pops_df = derivs_pops |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year)
derivs_pops_df_lower = derivs_pops_lower |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "lower")
derivs_pops_df_upper = derivs_pops_upper |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "upper")
derivs_pops_df = full_join(derivs_pops_df, derivs_pops_df_lower) |>
  full_join(derivs_pops_df_upper)
# set first time step's derivative to 0
derivs_pops_df$value[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$lower[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$upper[which(derivs_pops_df$year == 1981)] <- 0
# clean up the species names
derivs_pops_df$name = gsub("_", " ", derivs_pops_df$name) |> stringr::str_to_sentence()
derivs_pops_df$name =  factor(derivs_pops_df$name,
                              levels = temporal_trend$species[order(temporal_trend$mu_deriv)])
saveRDS(derivs_pops_df, "outputs/derivs_pops_df.rds")

# Average derivative across all species per year -------------------------------

avg_deriv_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs, 2, mean, na.rm = T),
  "cilo" = apply(derivs, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs, 2, quantile, prob = .95, na.rm = T)
)
avg_deriv_trend[1,2:4] = 0
saveRDS(avg_deriv_trend, "outputs/gam_hierarchical_df_overall_deriv.rds")

# Ridgeplot of derivatives per species -----------------------------------------

(plot_deriv_densities = ggplot(data = derivs_pops_df) +
   ggridges::geom_density_ridges_gradient(aes(x = value, y = name, fill = after_stat(x)),
                                          size = .2, scale = 3) +
   scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-.5,.5)) +
   geom_vline(xintercept = 0, lwd = .3, lty = 2) +
   theme(panel.grid.major.x = element_line(),
   ) +
   labs(x = "Annual rate of change (1981-2013)",
        y = "",
        fill = "")) +
  coord_cartesian(xlim = c(-.5,.5))
ggsave("figures/distribution_poptrends_ridgeplot.png", width = 8.23, height = 7)