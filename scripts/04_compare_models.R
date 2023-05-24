# Script to compare the models

library(brms)
library(HMSC)
library(dplyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())

# loading the data ----

biomass_l <- readRDS("outputs/biomass_l.rds")


## comparing the variation explained ----
# (i.e. which one has less residual variance)

## loading the models

linear_single <- lapply(list.files("outputs/linear_singlepop/")[-c(1,32)], function(x) readRDS(paste0("outputs/linear_singlepop/", x)))
linear_multi <- readRDS("outputs/linear_multipop.RDS")
linear_hier <- readRDS("outputs/linear_hierarchical.RDS")

## get the variance explained by the models

# each pop biomass ~ time
R2_single = lapply(linear_single, loo_R2) %>% lapply(as.data.frame) %>% bind_rows()
R2_single$pop = gsub(".rds", "", list.files("outputs/linear_singlepop/")[-c(1,32)])
# biomass ~ time
R2_multi = Rsquared(linear_multi)
# + a latent variable for year
R2_hier = Rsquared(linear_hier) # major improvement!!


## comparing the population trends ----

## loading the slope summaries
poptrends_single <- readRDS("outputs/linear_singlepop_population_trends.rds") #---- to fix
poptrends_multi = readRDS("outputs/linear_multipop_population_trends.rds")
poptrends_hier = readRDS("outputs/linear_hierarchical_population_trends.rds")

# join together into a long data frame
poptrends_single$model = "univariate"
poptrends_multi$model = "multivariate"
poptrends_hier$model = "hierarchical"
population_trends = rbind(poptrends_single, poptrends_multi, poptrends_hier) |> as.data.frame()
# order the species names so the points show up nicely
population_trends$pop = factor(population_trends$pop,
                               levels = poptrends_hier$pop[order(poptrends_hier$x1)])
(plot_poptrends = 
    ggplot(data = population_trends) +
    geom_segment(aes(x = x1 - x1_SD,
                     xend = x1 + x1_SD,
                     y = pop, yend = pop, col = model)) +
    geom_point(aes(y = pop, x = x1, col = model), size = 3) +
    labs(y = "", x = "temporal slope") +
    geom_vline(xintercept = 0, lwd = .3) +
    coord_cartesian(xlim = c(-1, 1))) +
  scale_color_manual(values = c("#cc5c76", "#1d457f", "#f9ad2a"))


## comparing the model's overall trend prediction ----

## loading the slope summaries
pred_single = readRDS("outputs/linear_singlepop_pred_l.rds")
pred_multi = readRDS("outputs/linear_multipop_pred_l.rds")
pred_hier = readRDS("outputs/linear_hierarchical_pred_l.rds")
df_overall_single <- readRDS("outputs/linear_singlepop/avg_trend.rds")
df_overall_multi = readRDS("outputs/linear_multipop_df_overall.RDS")
df_overall_hier = readRDS("outputs/linear_hierarchical_df_overall.RDS")

ggplot() +
  # raw data
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = pop), 
            lwd = .3, col = "grey85") +
  # average of single-population models
  geom_ribbon(data = df_overall_single, aes(x = year, ymin = cilo, ymax = cihi), 
              alpha = .3, fill = "goldenrod1") +
  geom_line(data = df_overall_single, aes(x = year, y = biomass), 
            linewidth = 1, col = "goldenrod1") +
  # multipopulation model without year effect
  geom_ribbon(data = df_overall_multi, aes(x = year, ymin = q05, ymax = q95), 
              alpha = .6, fill = "#cc5c76") +
  geom_line(data = df_overall_multi, aes(x = year, y = mu), 
            linewidth = 1, col = "#cc5c76") +
  # hierarchical model
  geom_ribbon(data = df_overall_hier, aes(x = year, ymin = q05, ymax = q95), 
              alpha = .6, fill = "#1d457f") +
  geom_line(data = df_overall_hier, aes(x = year, y = mu), 
            linewidth = 1, col = "#1d457f") +
  labs(x = "Year", y = "Biomass")
# note: the univariate and multivariate slopes overlap completely

## compare distribution of slope estimates ====

linear_multi$results$estimation$paramX[,2,] |> 
  density() |> 
  plot(col = "#cc5c76", frame = F, main = "Distribution of temporal slopes", ylim = c(0, 1.8))
m = linear_multi
# average slope with CI
paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.05, .5, .95))
abline(v = paramX_overall_mean[2], col = "#cc5c76", lwd = 2)
abline(v = paramX_overall_quantile[1,2], col = "#cc5c76", lty = 2)
abline(v = paramX_overall_quantile[3,2], col = "#cc5c76", lty = 2)

linear_hier$results$estimation$paramX[,2,] |> density() |> lines(col =  "#1d457f")
m = linear_hier
# average slope with CI
paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.05, .5, .95))
abline(v = paramX_overall_mean[2], col = "#1d457f", lwd = 2)
abline(v = paramX_overall_quantile[1,2], col = "#1d457f", lty = 2)
abline(v = paramX_overall_quantile[3,2], col = "#1d457f", lty = 2)

singlepop_all_draws <- readRDS("~/Documents/GitHub/hierarchical-lpi/outputs/linear_singlepop/all_draws.rds")
# singlepop_all_draws[,2] |> density() |> lines(col =  "goldenrod1") # the density is super high, so comment out to see the other distributions better
# average slope with CI
paramX_overall_mean = singlepop_all_draws[,2] |> mean()
paramX_overall_sd = singlepop_all_draws[,2] |> sd()
paramX_overall_quantile = singlepop_all_draws[,2] |> quantile(prob = c(.05, .5, .95))
abline(v = paramX_overall_mean, col = "goldenrod1", lwd = 2)
abline(v = paramX_overall_quantile[1], col = "goldenrod1", lty = 2)
abline(v = paramX_overall_quantile[3], col = "goldenrod1", lty = 2)
