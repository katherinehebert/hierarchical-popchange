# Script to compare the models

library(brms)
library(bayesGAM)
library(mvgam)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# loading the data ----

biomass_l <- readRDS("outputs/biomass_l.rds")


## comparing the variation explained ----
# (i.e. which one has less residual variance)

## loading the models

gam_single <- lapply(list.files("outputs/gam_singlepop/")[-c(1,32,33)], function(x) readRDS(paste0("outputs/gam_singlepop/", x)))
gam_multi <- readRDS("outputs/gam_multipop.RDS")
gam_multi_latent <- readRDS("outputs/gam_multipop_latent.RDS")
gam_hier <- readRDS("outputs/gam_hierarchical.RDS")
gam_hier_latent <- readRDS("outputs/gam_hierarchical_latent.RDS")

## comparing the population trends ----

## loading the slope summaries
poptrends_single <- readRDS("outputs/gam_singlepop_population_trends.rds") #---- to fix
poptrends_multi = readRDS("outputs/gam_multipop_population_trends.rds")
poptrends_multi_latent = readRDS("outputs/gam_multipop_latent_population_trends.rds")
poptrends_hier = readRDS("outputs/gam_hierarchical_population_trends.rds")
poptrends_hier_latent = readRDS("outputs/gam_hierarchical_latent_population_trends.rds")

# join together into a long data frame
poptrends_single$model = "univariate"
poptrends_multi$model = "multivariate"
poptrends_multi$model = "multivariate-latent"
poptrends_hier$model = "hierarchical"
poptrends_hier_latent$model = "hierarchical-latent"

population_trends = rbind(poptrends_single, 
                          poptrends_multi,
                          poptrends_hier,
                          poptrends_hier_latent) |> as.data.frame()
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
    #coord_cartesian(xlim = c(-1, 1))) +
  scale_color_manual(values = c("#cc5c76", "#1d457f", "#f9ad2a", "#625a94", "seagreen")))


## comparing the model's overall trend prediction ----

## loading the slope summaries
pred_single = readRDS("outputs/gam_singlepop_pred_l.rds")
pred_multi = readRDS("outputs/gam_multipop_pred_l.rds")
pred_multi_latent = readRDS("outputs/gam_multipop_latent_pred_l.rds")
pred_hier = readRDS("outputs/gam_hierarchical_pred_l.rds")
pred_hier_latent = readRDS("outputs/gam_hierarchical_latent_pred_l.rds")

df_overall_single <- readRDS("outputs/gam_singlepop/avg_trend.rds")
df_overall_multi = readRDS("outputs/gam_multipop_df_overall.rds")
df_overall_multi_latent = readRDS("outputs/gam_multipop_latent_df_overall.rds")
df_overall_hier = readRDS("outputs/gam_hierarchical_df_overall.rds")
df_overall_hier_latent = readRDS("outputs/gam_hierarchical_latent_df_overall.rds")

df_overall_multi_d = readRDS("outputs/gam_multipop_df_overall_deriv.rds")
df_overall_multi_latent_d = readRDS("outputs/gam_multipop_latent_df_overall_deriv.rds")
df_overall_hier_d = readRDS("outputs/gam_hierarchical_latent_df_overall_deriv.rds")
df_overall_hier_latent_d = readRDS("outputs/gam_hierarchical_latent_df_overall_deriv.rds")


A = ggplot() +
  # raw data
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = pop), 
            lwd = .3, col = "grey85") +
  # average of single-population models
  geom_ribbon(data = df_overall_single, aes(x = year, ymin = cilo, ymax = cihi), 
              alpha = .3, fill = "goldenrod1") +
  geom_line(data = df_overall_single, aes(x = year, y = biomass), 
            linewidth = 1, col = "goldenrod1") +
  labs(x = "Year", y = "Biomass")


B = ggplot() +
  # raw data
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = pop), 
            lwd = .3, col = "grey85") +
  # multipopulation model without year effect
  geom_ribbon(data = df_overall_multi, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "#cc5c76") +
  geom_line(data = df_overall_multi, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "#cc5c76") +
  # multipopulation-latent model without year effect
  geom_ribbon(data = df_overall_multi_latent, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "seagreen") +
  geom_line(data = df_overall_multi_latent, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "seagreen") +
  labs(x = "Year", y = "Biomass")

  
C = ggplot() +
  # raw data
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = pop), 
            lwd = .3, col = "grey85") +
  # hierarchical model
  geom_ribbon(data = df_overall_hier, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "#1d457f") +
  geom_line(data = df_overall_hier, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "#1d457f") +
  # hierarchical-latent model
  geom_ribbon(data = df_overall_hier_latent, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "#625a94") +
  geom_line(data = df_overall_hier_latent, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "#625a94") +
  labs(x = "Year", y = "Biomass")

A + B + C

## derivatives

B = ggplot() +
  # multipopulation model without year effect
  geom_ribbon(data = df_overall_multi_d, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "#cc5c76") +
  geom_line(data = df_overall_multi_d, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "#cc5c76") +
  # multipopulation-latent model without year effect
  geom_ribbon(data = df_overall_multi_latent_d, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "seagreen") +
  geom_line(data = df_overall_multi_latent_d, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "seagreen") +
  labs(x = "Year", y = "Biomass")


C = ggplot() +
  # hierarchical model
  geom_ribbon(data = df_overall_hier_d, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "#1d457f") +
  geom_line(data = df_overall_hier_d, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "#1d457f") +
  # hierarchical-latent model
  geom_ribbon(data = df_overall_hier_latent_d, aes(x = year, min = cilo, ymax = cihi), 
              alpha = .6, fill = "#625a94") +
  geom_line(data = df_overall_hier_latent_d, aes(x = year, y = avg_trend), 
            linewidth = 1, col = "#625a94") +
  labs(x = "Year", y = "Biomass")

B + C


## compare distribution of slope estimates ====
coef(gam_multi, summarise = FALSE)[,-1] |> 
   density() |> 
   plot(col = "#cc5c76", frame = F, main = "Distribution of temporal splines", ylim = c(0, 1.8))

coef(gam_multi_latent, summarise = FALSE)[,-1] |> 
  density() |> 
  lines(col =  "seagreen")
# # average slope with CI
# paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
# paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
# paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.05, .5, .95))
# abline(v = paramX_overall_mean[2], col = "#cc5c76", lwd = 2)
# abline(v = paramX_overall_quantile[1,2], col = "#cc5c76", lty = 2)
# abline(v = paramX_overall_quantile[3,2], col = "#cc5c76", lty = 2)
 
coef(gam_hier, summarise = FALSE)[,-1] |> 
  density() |> 
  lines(col =  "#1d457f")

coef(gam_hier_latent, summarise = FALSE)[,-1] |> 
  density() |> 
  lines(col =  "#625a94")
# m = gam_hier
# # average slope with CI
# paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
# paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
# paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.05, .5, .95))
# abline(v = paramX_overall_mean[2], col = "#1d457f", lwd = 2)
# abline(v = paramX_overall_quantile[1,2], col = "#1d457f", lty = 2)
# abline(v = paramX_overall_quantile[3,2], col = "#1d457f", lty = 2)


# gam_single[[1]]@model@names_beta
# gam_single[[1]]@model@names_u
# gam_single[[1]]@model@names_y
# gam_single[[1]]@results
#singlepop_all_draws <- readRDS("~/Documents/GitHub/hierarchical-lpi/outputs/gam_singlepop_alldraws.rds")


# singlepop_all_draws[,2] |> density() |> lines(col =  "goldenrod1") # the density is super high, so comment out to see the other distributions better
# average slope with CI
# paramX_overall_mean = singlepop_all_draws[,2] |> mean()
# paramX_overall_sd = singlepop_all_draws[,2] |> sd()
# paramX_overall_quantile = singlepop_all_draws[,2] |> quantile(prob = c(.05, .5, .95))
# abline(v = paramX_overall_mean, col = "goldenrod1", lwd = 2)
# abline(v = paramX_overall_quantile[1], col = "goldenrod1", lty = 2)
# abline(v = paramX_overall_quantile[3], col = "goldenrod1", lty = 2)

