# Script to make a multivariate hierarchical generalized additive model 
# (i.e. with a year effect) and with latent temporal processs

# libraries ----

library(here)
library(dplyr)
library(tidyr)
library(mvgam)
library(tidybayes)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())

set.seed(12)

# data ----

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)
time_m <- as.matrix(time)

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

biomass <- Year_Geom_Means_all %>% apply(2, scale, center = TRUE) |> as.data.frame()
matplot(biomass, type = "l")

# format into long
biomass$time = time
dat = pivot_longer(biomass, cols = -c(time), names_to = "series", values_to = "y")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time)
data_train = dat

################################################################################

# hierarchical gam on all populations ----

# prepare the priors ----
knots = ceiling(tsl/4)
mvgam_prior <- mvgam(data = data_train,
                     formula = y ~ s(time, bs = "tp", k = knots) + 
                       s(series, bs = 're', k = npops),
                     family = "gaussian",
                     trend_model = 'AR1',
                     use_lv =  TRUE, n_lv = 2,
                     chains = 3,
                     use_stan = TRUE,
                     prior_simulation = TRUE)
mvgam_prior$model_file

# record the priors
test_priors <- get_mvgam_priors(y ~ s(time, bs = "tp", k = knots) + 
                                  s(series, bs = 're', k = npops),
                                family = "gaussian",
                                data = data_train,
                                use_lv = TRUE, n_lv = 2,
                                trend_model = 'AR1',
                                use_stan = TRUE)
test_priors

# look at the priors
plot(mvgam_prior, type = 'smooths', realisations = TRUE)
plot(mvgam_prior, type = 'trend')


# train the model on data ----
mod1 <- mvgam(data = data_train,
              formula = y ~ s(time, bs = "tp", k = knots) + 
                s(series, bs = 're', k = npops), 
              use_lv = TRUE,
              n_lv = 2,
              family = "gaussian",
              trend_model = 'AR1', # latent temporal effect
              use_stan = TRUE,
              chains = 3,
              burnin = 100,
              samples = 1000
)
saveRDS(mod1, paste0("outputs/gam_hierarchical_latent.rds")) # a looot of divergence
# to view the Stan model file:
code(mod1)
m <- mod1

# convergence?
rstan::stan_trace(mod1$model_output, 'rho')
rstan::stan_trace(mod1$model_output, 'b')

# extract posterior draws in an array format
draws_fit = m$model_output |> posterior::as_draws_matrix()
posterior_df = posterior::summarize_draws(draws_fit)
saveRDS(posterior_df, "outputs/gam_hierarchical_latent_posterior.rds")

coefs = coef(m, summarise = FALSE)
coef_df = data.frame(
  "variable" = colnames(coefs),
  "mean" = apply(coefs, 2, mean),
  "q5" = apply(coefs, 2, quantile, prob = .05),
  "q95" = apply(coefs, 2, quantile, prob = .95)
)
saveRDS(coef_df, "outputs/gam_hierarchical_latent_coefs.rds")

# assign the population name to the numbers in the parameter names
# extract the model summary
summary_df = posterior_df[-c(grep("ypred", posterior_df$variable),
                             grep("mu_raw", posterior_df$variable),
                             grep("sigma_raw", posterior_df$variable)),]
pops = data.frame(
  "pop" = colnames(biomass)[-31],
  "number" = 1:(ncol(biomass)-1)
)
# get coefficients for each series (species)
coef_df$sd = apply(coefs, 2, sd)
coef_pop_df = coef_df[grep("series", coef_df$variable),] |>
  tidyr::separate(variable, sep = "\\.", into = c("variable", "number")) |>
  mutate(number = as.integer(number)) |>
  left_join(pops)
temp = coef_df[grep("Intercept", coef_df$variable),]
temp$number = NA
temp$pop = NA
temp = subset(temp, select = c(variable, number, mean, q5, q95, sd, pop))
coef_pop_df = rbind(coef_pop_df, temp)
coef_pop_df = coef_pop_df |> 
  subset(select = c(mean, sd, pop)) |>
  mutate("Intercept" = temp$mean,
         "Intercept_SD" = temp$sd) |>
  subset(select = c(Intercept, mean, Intercept_SD, sd, pop)) |>
  rename("x1" = "mean",
         "x1_SD" = "sd")
saveRDS(coef_pop_df, "outputs/gam_hierarchical_latent_population_trends.rds") 

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
  mutate(time = as.integer(time)+1980)
saveRDS(predictions, "outputs/gam_hierarchical_latent_pred_l.rds")


## average trend

source("~/Documents/GitHub/hierarchical-lpi/scripts/plot_mvgam_trend_custom.R")
trend_vals = list()
preds_ls = list()
derivs_ls = list()
for(i in 1:npops){
  trend_vals[[i]] = plot_mvgam_trend_custom(mod1, derivatives = TRUE, series = i)
  preds_ls[[i]] = trend_vals[[i]]$preds
  derivs_ls[[i]] = trend_vals[[i]]$derivs
}
preds = do.call(rbind, preds_ls)
derivs = do.call(rbind, derivs_ls)

avg_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds, 2, mean),
  "cilo" = apply(preds, 2, quantile, prob = .05),
  "cihi" = apply(preds, 2, quantile, prob = .95),
  "sd" = apply(preds, 2, sd)
)
saveRDS(avg_trend, "outputs/gam_hierarchical_latent_df_overall.rds")

avg_deriv_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs, 2, mean, na.rm = T),
  "cilo" = apply(derivs, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs, 2, quantile, prob = .95, na.rm = T),
  "sd" = apply(derivs, 2, sd)
)
saveRDS(avg_deriv_trend, "outputs/gam_hierarchical_latent_df_overall_deriv.rds")

plot(avg_trend ~ time, data = avg_trend, 
     type = "l", lty = 1, ylim = c(-3, 3),
     main = "average biomass trend")
lines(cilo ~ time, data = avg_trend, lty = 2)
lines(cihi ~ time, data = avg_trend, lty = 2)

plot(avg_trend ~ time, data = avg_deriv_trend, 
     type = "l", lty = 1, ylim = c(-1, 1), col = NULL,
     main = "average of the derivative trend")
abline(h = 0, col = "grey")
lines(avg_trend ~ time, data = avg_deriv_trend)
lines(cilo ~ time, data = avg_deriv_trend, lty = 2)
lines(cihi ~ time, data = avg_deriv_trend, lty = 2)

# get species correlations
sp_correlations = lv_correlations(mod1)
saveRDS(sp_correlations, "outputs/gam_hierarchical_latent_species_correlations.rds")

# Plot as heatmap
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = .4, tl.col = "black")

