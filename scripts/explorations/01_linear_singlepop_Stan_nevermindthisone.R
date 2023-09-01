# Script to make a linear model per population, then average them
# This is the current model we will compare to hierarchical models

# libraries ----

library(here)
library(dplyr)
library(cmdstanr)
library(rstan)
library(tidybayes)
library(bayesplot)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

set.seed(12)

# data ----

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

biomass <- Year_Geom_Means_all %>% apply(2, scale, center = FALSE)
matplot(biomass, type = "l")

################################################################################

# simulate some priors ----

# need better priors...
prior_biomass = list()
for(i in 1:tsl){
  prior_alpha = rlnorm(1, 0, .5)
  prior_beta_year = rnorm(1, 0, .1)
  prior_sigma = rexp(1, 1)
  (prior_mus = prior_alpha + prior_beta_year * time)
  prior_biomass[[i]] = rnorm(100, mean = prior_mus[i], sd = prior_sigma)
}
par(mfrow=c(1,2))
boxplot(prior_biomass)
boxplot(t(biomass))


################################################################################

# linear model per population ----

# compile the model coded in Stan
m = cmdstan_model(stan_file = "scripts/linear_singlepop.stan", pedantic = TRUE)

# fit the model with the data
m_fit = vector("list", npops)
for(i in 1:npops){
  # make a data list for the model
  m_data = list(
    n = tsl,
    biomass = biomass[,i],
    year = time,
    npost = tsl
  )
  # sample from the model
  m_fit[[i]] = m$sample(data = m_data)
  
  # save
  m_fit[[i]]$save_object(file = paste0("outputs/linear_singlepop/model_", i, ".rds"))
}


## posterior predictive checks ----

par(mfrow = c(3,2))
draws_matrix <- list()
for(i in 1:npops){
  draws <- m_fit[[i]]$draws(variables = c("y_rep"))
  draws_matrix[[i]] <- posterior::as_draws_matrix(draws)
  bayesplot::ppc_dens_overlay(y = biomass[,i],
                              yrep = head(draws_matrix[[i]], 50)) |> print()
}


## plot per population ----

posterior_ls = list()
par(mfrow = c(3,2))
for(i in 1:npops){
  m_fit_temp <- read_stan_csv(m_fit[[i]]$output_files())
  posterior <- as.data.frame(m_fit_temp)
  
  # plot the fitted Line
  p = ggplot() +
    geom_point(aes(x = m_data$year, y = biomass[,i])) +
    stat_function(fun = function(x) mean(posterior$alpha) + mean(posterior$beta_year) * x)
  print(p)
  
  posterior_ls[[i]] = posterior
}
names(posterior_ls) = colnames(biomass)

# prep per-population posterior summary as one long data frame for the next plot

posterior_pops  = data.frame(
  "sp" = names(posterior_ls),
  "slope" = unlist(lapply(posterior_ls, function(x) mean(x$beta_year))),
  "intercept" = unlist(lapply(posterior_ls, function(x) mean(x$alpha))),
  "slope_lo" = unlist(lapply(posterior_ls, function(x) quantile(x$beta_year, prob = .025))),
  "slope_hi" = unlist(lapply(posterior_ls, function(x) quantile(x$beta_year, prob = .975))),
  "intercept_lo" = unlist(lapply(posterior_ls, function(x) quantile(x$alpha, prob = .025))),
  "intercept_hi" = unlist(lapply(posterior_ls, function(x) quantile(x$beta_year, prob = .975)))
)

################################################################################

# calculate the average trend across populations ----

posterior <- bind_rows(posterior_ls, .id = "sp")
saveRDS(posterior, "outputs/linear_singlepop/posterior.rds")

slope = posterior$beta_year |> quantile(probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
intercept = posterior$alpha |> quantile(probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
sigma = posterior$sigma |> mean(na.rm = TRUE)

avg_trend = data.frame(
  "year" = time+1981,
  "biomass" = slope[2]*time + intercept[2],
  "cilo" = slope[1]*time + intercept[1],
  "cihi" = slope[3]*time + intercept[3]
)
saveRDS(avg_trend, "outputs/linear_singlepop/avg_trend.rds")

# predict
pop_trends = list()
for(i in 1:npops){
  pop_trends[[i]] = data.frame(
    "year" = time+1981,
    "biomass" = posterior_pops$slope[i]*time + posterior_pops$intercept[i],
    "cilo" = posterior_pops$slope_lo[i]*time + posterior_pops$intercept_lo[i],
    "cihi" = posterior_pops$slope_hi[i]*time + posterior_pops$intercept_hi[i]
  )
}
names(pop_trends) = colnames(biomass)
pop_trends = bind_rows(pop_trends, .id = "sp")

## plot the average trend  ----

# convert biomass to long
biomass_l = as.data.frame(biomass) |> 
  mutate("year" = time+1981) |>
  tidyr::pivot_longer(cols = -year, names_to = "sp", values_to = "biomass")

ggplot() +
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = sp), size = .3, alpha = .4) +
  geom_ribbon(data = avg_trend,
            aes(x = year, ymin = cilo, ymax = cihi), alpha = .3) + # credible interval?
  geom_line(data = avg_trend,
              aes(x = year, y = biomass), linewidth = 1, col = "red") +
  # geom_ribbon(data = pop_trends,
  #           aes(x = year, ymin = cilo, ymax = cihi, group = sp), alpha = .05) +
  geom_line(data = pop_trends,
            aes(x = year, y = biomass, group = sp), linewidth = .1, col = "blue")
