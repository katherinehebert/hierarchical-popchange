# Inspect the latent factors

library(mvgam)
library(rstan)
library(patchwork)
library(tidyverse)

biomass = readRDS("outputs/biomass_scaling.rds")
biomasses = as.data.frame(biomass$YMeans)
colnames(biomasses) = "biomass_mean"

derivs_pops_df = readRDS("outputs/derivs_pops_df.rds") |> 
  #dplyr::filter(year > 1998) |>
  dplyr::group_by(name) |>
  dplyr::summarise(mean_rate = mean(value),
                   sd_rate = sd(value))

mod1 = readRDS(paste0("outputs/gam_hierarchical_gp.rds"))

variables(mod1) |> names()

mvgam::plot_mvgam_factors(mod1)
mod1$n_lv
preds <- mvgam:::mcmc_chains(mod1$model_output, "LV")
loadings <- mvgam:::mcmc_chains(mod1$model_output, "trend")
posterior_df = readRDS("outputs/gam_hierarchical_posterior.rds")

loadings = posterior_df[grep("lv", posterior_df$variable),] 
lv = posterior_df[grep("LV", posterior_df$variable),] 

loadings = loadings |> tidyr::separate(variable, sep = ",", into = c("species", "LV"))
loadings$species = colnames(biomass_scaling$biomass)[readr::parse_number(loadings$species)]
loadings$species = gsub("_", " ", loadings$species) |> stringr::str_to_sentence()
loadings$LV = readr::parse_number(loadings$LV)

loadings_plot =
  data.frame(
    "name" = loadings$species[which(loadings$LV == 1)],
    "LV1" = loadings$mean[which(loadings$LV == 1)],
    "LV2" = loadings$mean[which(loadings$LV == 2)]
  ) 
loadings_plot = dplyr::full_join(loadings_plot, derivs_pops_df)

biomasses$name <- loadings_plot$name
loadings_plot = dplyr::full_join(loadings_plot, biomasses)

library(ggplot2)
ggplot(data = loadings_plot) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_label(aes(x = LV1, y = LV2, 
                 fill = mean_rate, label = gsub(" ", "\n", name)), 
             col = "black") +
  theme_minimal() +
  scale_fill_distiller(palette = "Spectral", direction = 1)

ggplot(data = loadings_plot) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(x = LV1, y = LV2, fill = mean_rate), pch = 21, size = 4) +
  ggpubr::theme_pubr()  +
  scale_fill_distiller(palette = "Spectral", direction = 1, limits = c(-0.09, 0.09)) +
  labs(fill = "Average\ngrowth rate", x = "Latent variable 1", y = "Latent variable 2") +
  theme(legend.position  = "right") +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-2.5, 2.5))
ggsave("figures/ordinationplot_latents.png", height = 5, width = 6)


ggplot(data = loadings_plot) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(x = LV1, y = LV2, fill = sd_rate), pch = 21, size = 4) +
  ggpubr::theme_pubr()  +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(fill = "Average\ngrowth rate", x = "Latent variable 1", y = "Latent variable 2") +
  theme(legend.position  = "right") +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-2.5, 2.5))
ggsave("figures/ordinationplot_latents.png", height = 5, width = 6)

plot(loadings_plot$sd_rate/loadings_plot$mean_rate ~ loadings_plot$LV2)
m1 = lm(loadings_plot$sd_rate ~ loadings_plot$LV2)
summary(m1)
abline(a = 0.076, b = 0.01045)

# try to correlate latents with fisheries and ????

fishing_summary <- readr::read_csv("~/Documents/GitHub/groundfish-data-analysis/data/fishing_effort_data.csv")
fishing_summary <- fishing_summary[which(fishing_summary$year > 1981) & which(fishing_summary$year < 2015),]
fishing_summary <- dplyr::filter(fishing_summary, type == "benthic")
fishing_effort_data = fishing_summary$catch[fishing_summary$type=="benthic"]/1e6
fishing_summary$catch = fishing_effort_data

lv = tidyr::separate(lv, col = variable, into = c("time", "LV"), sep = ",", remove = TRUE)
lv$time = readr::parse_number(lv$time)
lv$LV = readr::parse_number(lv$LV)

A = ggplot(data = fishing_summary) +
  geom_line(aes(x = year, y = catch)) +
  labs(x = "Year", y = "Catch", title = "Potential drivers") +
  ggpubr::theme_pubr()


B = ggplot(data = dplyr::filter(lv, LV == 1)) +
  geom_ribbon(aes(x = time+1981, ymin = q5, ymax = q95), alpha = .3, fill = "navyblue") +
  geom_line(aes(x = time+1981, y = median)) +
  labs(x = "Year", y = "Latent variable 1", title = "Latent variation") +
  ggpubr::theme_pubr()

A + B

climate_index <- readr::read_csv("~/Documents/GitHub/groundfish-data-analysis/data/composite_index.csv")

C = ggplot(data = dplyr::filter(climate_index, year > 1981 & year <2015)) +
  geom_line(aes(x = year, y = climate_index)) +
  labs(x = "Year", y = "Composite climate index") +
  ggpubr::theme_pubr()

D = ggplot(data = dplyr::filter(lv, LV == 2)) +
  geom_ribbon(aes(x = time+1981, ymin = q5, ymax = q95), alpha = .3, fill = "navyblue") +
  geom_line(aes(x = time+1981, y = median)) +
  labs(x = "Year", y = "Latent variable 2") +
  ggpubr::theme_pubr()

(A + B) / (C + D) + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave("figures/latent_potentialdrivers.png", width = 7, height = 7)
