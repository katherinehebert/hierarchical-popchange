# Script to inspect and reformat the outputs from the multivariate 
# hierarchical generalized additive model 

# libraries ----

library(here)
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

# read the previous data objects and the model object ----
biomass_scaling = readRDS("outputs/biomass_scaling.rds")
mod1 = readRDS(paste0("outputs/gam_hierarchical_gp_2025.rds"))
biomass = biomass_scaling$biomass

# to view the Stan model file:
m <- mod1
# save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(m)

# convergence
png("figures/traceplots_timebases.png", height=1800, width=1800, type = "cairo")
rstan::stan_trace(mod1$model_output, 'b')
dev.off()

## forecast ----
newdata = mod1$obs_data
newdata$time = newdata$time + 32
fore = mvgam::forecast(mod1, newdata = newdata)
saveRDS(fore, "outputs/forecast.rds")
plot(fore,series = 10)
forescore = score(fore)

## Plot the temporal trend (time smooth) without population random effect----

png("figures/timesmooth.png", height=1000, width=1000, type = "cairo")
plot_mvgam_smooth(mod1, smooth = "time") 
dev.off()

# smooth for all species is the same (global smooth)
par(mfrow = c(3,3))
for(i in 1:npops){
  plot_mvgam_smooth(mod1, smooth = "time", series = i)
}

# Plot the random effect on species
png("figures/mvgam_randomeffect.png", height=1000, width=1800, type = "cairo",pointsize = 20)
mvgam::plot_mvgam_randomeffects(mod1)
dev.off()

# Plot the latent factors
png("figures/mvgam_latentfactorplot.png", height=1000, width=1800, type = "cairo",pointsize = 20)
plot_mvgam_factors(mod1)
dev.off()

# Plot the ppc
pdf("figures/mvgam_pp_check.pdf",onefile = TRUE)
  pp_check(mod1, type = "ecdf_overlay",
           newdata = list(
    time = mod1$obs_data$time,
    series = factor(mod1$obs_data$series),
    y = mod1$obs_data$y,
    knots = 9,
    npops = 29))
  
  pp_check(mod1, type = "resid_ribbon",
           newdata = list(
             time = mod1$obs_data$time,
             series = factor(mod1$obs_data$series),
             y = mod1$obs_data$y,
             knots = 9,
             npops = 29))
  
dev.off()

# Plot the ppc
pdf("figures/mvgam_ppc.pdf",onefile = TRUE)
for(i in 1:ncol(biomass)){
ppc(mod1,series = i)
}
dev.off()


# Plot the ppc
pdf("figures/mvgam_ppc_cdf.pdf",onefile = TRUE)
for(i in 1:ncol(biomass)){
  ppc(mod1,series = i, type = "cdf")
}
dev.off()

# Plot the ppc
pdf("figures/mvgam_ppc_mean.pdf",onefile = TRUE)
for(i in 1:ncol(biomass)){
  ppc(mod1,series = i, type = "mean")
}
dev.off()

ppc(mod1,type = "mean")

# plot(mod1)

# extract posterior draws in an array format
draws_fit = m$model_output |> posterior::as_draws_matrix() # mvgam_draws() now exists too
posterior_df = posterior::summarize_draws(draws_fit)
saveRDS(posterior_df, "outputs/gam_hierarchical_posterior_2025.rds")

coefs_summ = coef(m, summarise = TRUE) |> as.data.frame()
saveRDS(coefs_summ, "outputs/gam_hierarchical_coefs_summary.rds")

# Summarise the coefficients
coefs = coef(m, summarise = FALSE) 
coef_df = data.frame(
  "variable" = colnames(coefs),
  "mean" = apply(coefs, 2, mean),
  "q5" = apply(coefs, 2, quantile, prob = .05),
  "q95" = apply(coefs, 2, quantile, prob = .95)
)
saveRDS(coef_df, "outputs/gam_hierarchical_coefs.rds")

# extract the model summary
summary_df = posterior_df[-c(grep("ypred", posterior_df$variable),
                             grep("mu_raw", posterior_df$variable),
                             grep("sigma_raw", posterior_df$variable)),]

# assign the population name to the numbers in the parameter names
pops = data.frame(
  "pop" = colnames(biomass)[-ncol(biomass)],
  "number" = 1:(ncol(biomass)-1)
)

# predict from the model
predictions = posterior_df[c(grep("ypred", posterior_df$variable)),] |>
  tidyr::separate(variable, sep = "\\[", into = c("variable", "pop_number")) |>
  tidyr::separate(pop_number, sep = ",", into = c("time", "pop")) |>
  mutate(pop = readr::parse_number(pop)) |>
  left_join(pops, by = c("pop" = "number"), keep = FALSE) |>
  dplyr::select(c(pop.y, time, mean, q5, q95)) |>
  rename("sp" = "pop.y",
         "biomass" = "mean",
         "cilo" = "q5",
         "cihi" = "q95") |>
  mutate(time = as.integer(time)+1981)
saveRDS(predictions, "outputs/gam_hierarchical_pred_l.rds")
