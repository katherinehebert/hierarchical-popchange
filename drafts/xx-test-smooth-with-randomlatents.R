# Script to test the effect of the order of entry of the predictors and the
# latent variables in the mvgam

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

# data -------------------------------------------------------------------------

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)
time_m <- as.matrix(time)

# remove Notacanthidae, which was 0 in year 1
Year_Geom_Means_all = Year_Geom_Means_all[,-which(colnames(Year_Geom_Means_all) %in% "NOTACANTHIDAE")]

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

YData <- Year_Geom_Means_all 
# center on the baseline biomass and standardise by the mean
for(i in 1:ncol(YData)){
  YData[,i] = (YData[,i] - YData[1,i])/(mean(YData[,i], na.rm = TRUE))
}

biomass <- YData |> as.data.frame()
YMeans = apply(Year_Geom_Means_all, 2, mean, na.rm = TRUE)
baselines = Year_Geom_Means_all[1,] |> as.matrix() |> as.vector()
weightings = YMeans/sum(YMeans)*100

# format into long
biomass$time = time
dat = pivot_longer(biomass, cols = -c(time), names_to = "series", values_to = "y")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time)
data_train = dat

knots = ceiling(tsl/4)

################################################################################
# 1. only random effect and latent variables -----------------------------------
# random walk only

# train the model
mod1 <- mvgam(data = data_train,
              formula =  y ~ s(time, bs = "fs", k = knots) + s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'RW',
              use_stan = TRUE,
              chains = 3, 
              burnin = 10,
              samples = 500
)
#saveRDS(mod1, paste0("outputs/testing/gam_04_smoothwithrandomlatents.rds")) 


plot_mvgam_smooth(mod1, residuals = TRUE)


# Plot latent factors ----------------------------------------------------------

png("outputs/testing/fig_gam_04_smoothwithrandomlatents_latentfactors.png", height=1000, width=1800, type = "cairo",pointsize = 20)
mvgam::plot_mvgam_factors(mod1)
dev.off()

# Plot the random effect on species --------------------------------------------

png("outputs/testing/fig_gam_04_smoothwithrandomlatents_randomeffect.png", height=1000, width=1800, type = "cairo",pointsize = 20)
mvgam::plot_mvgam_randomeffects(mod1)
dev.off()

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
saveRDS(preds, "outputs/testing/gam_04_smoothwithrandomlatents_predictions.rds")
saveRDS(derivs, "outputs/testing/gam_04_smoothwithrandomlatents_derivatives.rds")

# plot a histogram of the average derivatives across all species ---------------

# each species' median derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE)))
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))

# format into a long df
colnames(derivs_pops) = colnames(biomass)[-ncol(biomass)]
colnames(derivs_pops_lower) = colnames(biomass)[-ncol(biomass)]
colnames(derivs_pops_upper) = colnames(biomass)[-ncol(biomass)]
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
saveRDS(derivs_pops_df, "outputs/testing/gam_04_smoothwithrandomlatents_derivs_pops_df.rds")

derivs_without1981 = dplyr::filter(derivs_pops_df, year != 1981)

# Average derivative across all species per year -------------------------------

avg_deriv_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs, 2, mean, na.rm = T),
  "cilo" = apply(derivs, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs, 2, quantile, prob = .95, na.rm = T)
)
avg_deriv_trend[1,2:4] = 0
saveRDS(avg_deriv_trend, "outputs/testing/gam_04_smoothwithrandomlatents_df_overall_deriv.rds")

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
saveRDS(temporal_trend_yearly, "outputs/testing/gam_04_smoothwithrandomlatents_temporal_trend_yearly.rds")

################################################################################
# 3-panel plot: histogram, mean trend, coherence -------------------------------
################################################################################

## Panel A - Histogram of all population growth rates of the full time series ----
(plot_trenddensity = 
   ggplot(data = derivs_without1981) +
   geom_histogram(aes(x = value, fill = after_stat(x)), 
                  col = "black", linewidth = .2, bins = 19) + 
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE)) +
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE) - sd(derivs_without1981$value, na.rm = TRUE), lty = 2) +
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE) + sd(derivs_without1981$value, na.rm = TRUE), lty = 2) +
   theme(panel.grid.major.x = element_line()) +
   scale_y_sqrt() +
   labs(x = "Annual rate of change (α)", 
        y = "Frequency", 
        fill = "α") +
   scale_fill_distiller(palette = "RdYlGn", 
                        direction = 1, 
                        limits = c(-1,1)) +
   coord_cartesian(xlim = c(-1, 1))) 
# variance
mean(derivs_without1981$value, na.rm = TRUE)
var(derivs_without1981$value, na.rm = TRUE)

## Panel B - Trend of the distribution of rates of change through time ---------

(avg_derivative_pointplot = 
   ggplot(data = avg_deriv_trend, aes(x = year)) +
   ggpattern::geom_area_pattern(aes(y = cihi),
                                pattern = "gradient", 
                                fill = "#ffffbf",
                                pattern_fill  = "#ffffbf80",
                                pattern_fill2 = "#77c15c") + #prepended w/ 50% transparency hex code (80)
   ggpattern::geom_area_pattern(aes(y = cilo),
                                pattern = "gradient", 
                                fill = "#ffffbf",
                                pattern_fill  = "#f98c59",
                                pattern_fill2 = "#ffffbf80") + #prepended w/ 50% transparency hex code (80)
   geom_line(aes(y = avg_trend, col = avg_trend), lwd = .1) +
   geom_line(aes(y = avg_trend), lwd = .8, col = "black") +
   geom_line(aes(y = cilo), lwd = 1, col = "white") +
   geom_line(aes(y = cihi), lwd = 1, col = "white") +
   geom_hline(yintercept = 0, lwd = .3) +
   scale_color_gradientn(colours = c("#d73027", "#f98c59", "#fee08b", "#ffffbf", "#d9ef8b", "#77c15c", "#1a9850"), limits = c(-0.4, 0.4)) +
   labs(x = "",
        col = "α",
        y = "Rate of change (α)") +
   coord_cartesian(ylim = c(-1, 1)) +
   theme(panel.grid.major = element_line(),
         legend.position = "right") +
   scale_y_continuous(labels = scales::percent))
ggsave("outputs/testing/fig_04_smooth-withrandomlatents_assemblagevariability.png", width = 7, height = 3)

## Panel C - Trend of the distribution's variance through time -----------------

(plot_coherence = ggplot(data = temporal_trend_yearly) +
   geom_line(aes(y = sd, x = as.numeric(year), col = sd), lwd = 1.5) +
   geom_point(aes(y = sd, x = as.numeric(year), col = sd), 
              size = .9) +
   labs(y = "Community variability (σ)", 
        x = "",
        col = "σ") +
   scale_color_distiller(palette = "YlGnBu", direction = 1, limits = c(0, .4)) +
   coord_cartesian(ylim = c(0, .4)) +
   theme_pubr() +
   theme(panel.grid.major = element_line())
)

# assemble all the panels
# ((plot_trenddensity + 
#     theme(legend.position = "right")) /
#     (avg_derivative_pointplot) /
#     (plot_coherence + 
#        theme(legend.position = "right"))) +
#   plot_annotation(tag_levels = "a")
#ggsave("outputs/testing/fig_03_onlyrandomlatents_assemblagevariability.png", width = 7.5, height = 7)
