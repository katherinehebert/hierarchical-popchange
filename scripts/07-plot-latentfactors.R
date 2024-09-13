# Inspect the latent factors

library(mvgam)

biomass_scaling = readRDS("outputs/biomass_scaling.rds")
derivs_pops_df = readRDS("outputs/derivs_pops_df.rds") |> 
  dplyr::filter(year > 1998) |>
  dplyr::group_by(name) |>
  dplyr::summarise(mean_rate = mean(value))

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
  geom_point(aes(x = LV1, y = LV2, fill = mean_rate), pch = 21, size = 3) +
  theme_minimal()  +
  scale_fill_distiller(palette = "Spectral", direction = 1)


mod1$model_output |> View()

draws_fit = mod1$model_output |> posterior::as_draws_matrix() # mvgam_draws() now exists too
posterior_df = posterior::summarize_draws(draws_fit)


matplot(preds)

mvgam::plot_mvgam_trend(mod1, series = 1)

conditional_effects(mod1)


mod1$model_output |> dplyr::filter()

# this is the "LPI", should take the derivative of just this. this is the 
# common trajectory of the whole community.
plot(mod1, type = "smooth")

plot(mod1, type = "factors")

plot(mod1, type = "residuals")

