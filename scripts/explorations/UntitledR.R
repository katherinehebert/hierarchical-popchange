# Script to make a multivariate hierarchical generalized additive model 
# with an AR1 latent factor model

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

# data ----

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)
time_m <- as.matrix(time)

# remove taxa that are not IDed to species level
Year_Geom_Means_all = Year_Geom_Means_all[,-c(colnames(Year_Geom_Means_all) %in% c("NOTACANTHIDAE", "EUMICROTREMUS SP"))]

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

YData <- Year_Geom_Means_all 
# center on the baseline biomass and standardise by the mean
for(i in 1:ncol(YData)){
  YData[,i] = (YData[,i] - YData[1,i])/(mean(YData[,i], na.rm = TRUE))
}

biomass <- YData |> as.data.frame()
matplot(biomass, type = "l")

YMeans = apply(Year_Geom_Means_all, 2, mean, na.rm = TRUE)

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
                     formula = y ~ 
                       # global smoother for all pops over time
                       s(time, bs = "tp", k = knots) + 
                       # random intercept per group
                       s(series, bs = 're', k = npops),
                     family = "gaussian",
                     trend_model = 'AR1',
                     chains = 3,
                     use_stan = TRUE,
                     prior_simulation = TRUE)

mvgam_prior$model_file

# record the priors
test_priors <- get_mvgam_priors(y ~ 
                                  # global smoother for all pops over time
                                  s(time, bs = "tp", k = knots) + 
                                  # random intercept per group
                                  s(series, bs = 're', k = npops),
                                family = "gaussian",
                                data = data_train,
                                trend_model = 'AR1',
                                use_stan = TRUE)

# look at the priors
plot(mvgam_prior, type = 'smooths', realisations = TRUE)
plot(mvgam_prior, type = 'trend')
plot(mvgam_prior, type = 're')

# train the model on data ----
mod1 <- mvgam(data = data_train,
              formula =  y ~ s(time, bs = "tp", k = knots) + 
                s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'AR1',
              use_stan = TRUE,
              chains = 3, 
              burnin = 5000,
              samples = 10000
)
saveRDS(mod1, paste0("outputs/gam_hierarchical_ar1.rds")) 
m <- mod1
# execution time: 2100.7 seconds

# to save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(mod1)

# plot the model residuals
plot(mod1)

# visually check for convergence
png("figures/traceplots_timebases.png", height=1800, width=1800, type = "cairo")
rstan::stan_trace(mod1$model_output, 'b')
dev.off()

# plot the time smooth (without the species random effect) ----

png("figures/timesmooth.png", height=1000, width=1000, type = "cairo")
plot_mvgam_smooth(mod1, smooth = "time") # this is essentially the LPI curve
dev.off()

# plot the random intercept effect
mvgam::plot_mvgam_randomeffects(mod1)

png("figures/mvgam_latentfactorplot.png", height=1000, width=1800, type = "cairo",pointsize = 20)
plot_mvgam_factors(mod1)
dev.off()


# extract posterior draws ----

draws_fit = m$model_output |> posterior::as_draws_matrix() # mvgam_draws() now exists too
posterior_df = posterior::summarize_draws(draws_fit)
saveRDS(posterior_df, "outputs/gam_hierarchical_posterior.rds")

# summarise model coefficients ----

coefs_summ = coef(m, summarise = TRUE) |> as.data.frame()
saveRDS(coefs_summ, "outputs/gam_hierarchical_coefs_summary.rds")

# # assign the population name to the numbers in the parameter names
pops = data.frame(
  "pop" = colnames(biomass)[-30],
  "number" = 1:(ncol(biomass)-1)
)

################################################################################

# PREDICTIONS ----

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


# plot the predictions next to the data ----

# get the species names ready for plotting
data_train = rename(data_train, "sp" = "series")
data_train$sp = gsub("_", " ", data_train$sp) |> stringr::str_to_sentence()
data_train$sp = gsub(" ", "\n", data_train$sp)
predictions$sp = gsub("_", " ", predictions$sp) |> stringr::str_to_sentence()
predictions$sp = gsub(" ", "\n", predictions$sp)

# plot the predictions on the data
ggplot(data = predictions,
       aes(x = time, group = sp)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = sp), alpha = .2) +
  geom_line(aes(y = biomass)) +
  geom_line(data = data_train,
            aes(x = time+1981, y = y), lty = 2) +
  facet_wrap(~sp, scales = "free_y") +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic")) +
  labs(x = "", y = "Biomass (kg per tow)")
ggsave("figures/mvgam_predictions_perspecies.png", width = 12.8, height = 7.63)


# rescale the predictions to the response scale ----

predictions_respscale = predictions |> group_by(sp) |> group_split()
for(i in 1:length(predictions_respscale)){
  predictions_respscale[[i]]$biomass = predictions_respscale[[i]]$biomass*YMeans[i]
  predictions_respscale[[i]]$cilo = predictions_respscale[[i]]$cilo*YMeans[i]
  predictions_respscale[[i]]$cihi = predictions_respscale[[i]]$cihi*YMeans[i]
}
predictions_respscale = bind_rows(predictions_respscale)
saveRDS(predictions_respscale, "outputs/gam_hierarchical_pred_responsescale_l.rds")

# plot the rescaled predictions
(plot_respscale_predictions = 
    ggplot(data = predictions_respscale) +
    geom_ribbon(aes(x = time,
                    ymin = cilo,
                    ymax = cihi,
                    group = sp), alpha = .3, fill = "#6497b1") +
    geom_line(aes(x = time,
                  y = biomass,
                  group = sp)))


################################################################################

# RATE OF CHANGE (Derivatives) ----

# take the derivative of each population trend
source("~/Documents/GitHub/hierarchical-lpi/scripts/plot_mvgam_trend_custom.R")
trend_vals = list()
preds_ls = list()
preds_respscale_ls = list()
derivs_ls = list()
for(i in 1:npops){
  trend_vals[[i]] = plot_mvgam_trend_custom(mod1, derivatives = TRUE, series = i)
  preds_ls[[i]] = trend_vals[[i]]$preds
  preds_respscale_ls[[i]] = trend_vals[[i]]$preds*YMeans[i]
  derivs_ls[[i]] = trend_vals[[i]]$derivs
}
preds = do.call(rbind, preds_ls)
preds_respscale = do.call(rbind, preds_respscale_ls)
derivs = do.call(rbind, derivs_ls)

# each species' median predicted trend in biomass
preds_pops = do.call(cbind, 
                     lapply(preds_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(preds_pops, x = time+1981, type = "l")

preds_respscale_pops = do.call(cbind, 
                               lapply(preds_respscale_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(preds_respscale_pops, x = time+1981, type = "l")

# each species' median derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE)))
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))
matplot(derivs_pops, x = time+1981, type = "l")

# plot each species' mean derivative over the whole time series ----

(temporal_trend = data.frame(
  species = colnames(Year_Geom_Means_all),
  mu_deriv = derivs_pops |> apply(2, mean, na.rm = TRUE),
  lower = derivs_pops |> apply(2, quantile, probs = .05, na.rm = TRUE),
  upper = derivs_pops |> apply(2, quantile, probs = .95, na.rm = TRUE)
))
# make species names cleaner for plotting
temporal_trend$species = gsub("_", " ", temporal_trend$species) |> stringr::str_to_sentence()
temporal_trend$species =  factor(temporal_trend$species,
                                 levels = temporal_trend$species[order(temporal_trend$mu_deriv)])
colnames(derivs_pops) = colnames(Year_Geom_Means_all)
colnames(derivs_pops_lower) = colnames(Year_Geom_Means_all)
colnames(derivs_pops_upper) = colnames(Year_Geom_Means_all)

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

derivs_pops_df$name = gsub("_", " ", derivs_pops_df$name) |> stringr::str_to_sentence()
derivs_pops_df$name =  factor(derivs_pops_df$name,
                              levels = temporal_trend$species[order(temporal_trend$mu_deriv)])


# plot a histogram of the average derivatives across all species ----
(plot_trenddensity = 
   ggplot(derivs_pops_df) +
   geom_histogram(aes(x = value, fill = after_stat(x)), 
                  col = "black", linewidth = .2, bins = 19) + 
   geom_vline(xintercept = mean(derivs_pops, na.rm = TRUE)) +
   geom_vline(xintercept = mean(derivs_pops, na.rm = TRUE) - sd(derivs_pops, na.rm = TRUE), lty = 2) +
   geom_vline(xintercept = mean(derivs_pops, na.rm = TRUE) + sd(derivs_pops, na.rm = TRUE), lty = 2) +
   theme(panel.grid.major.x = element_line()) +
   scale_y_sqrt() +
   labs(x = "Annual rate of change", 
        y = "Frequency", 
        fill = "") +
   scale_fill_distiller(palette = "RdYlGn", 
                        direction = 1, 
                        limits = c(-1.1,1.1)) +
   coord_cartesian(xlim = c(-1.1, 1.1))) 

(plot_poptrends = ggplot(data = temporal_trend[order(temporal_trend$mu_deriv),]) +
    geom_segment(aes(x = lower, xend = upper,
                     y = species, yend = species), lwd = .3) +
    geom_point(aes(x = mu_deriv, y = species, fill = mu_deriv), 
               size = 4, pch = 21) +
    geom_vline(xintercept = 0, lwd = .3, lty = 2) +
    labs(x = "Average rate of change", 
         y = "",
         fill = "") +
    scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-.1,.1)) +
    coord_cartesian(xlim = c(-.5, .5)) +
    theme_pubr() +
    theme(axis.text.y = element_text(face = "italic"),
          panel.grid.major.x = element_line()) 
)

################################################################################

# Fig 2 - distribution of species' rates of change ----

# arrange the plot panels and save
(plot_trenddensity + 
   coord_cartesian(xlim = c(-1.1, 1.1)) +
   theme(axis.title.y = element_text(vjust = -60)) 
) / (plot_poptrends + 
       coord_cartesian(xlim = c(-.5, .5)) +
       theme(axis.text.y = element_text(face = "italic"),
             panel.grid.major.x = element_line())) + 
  plot_layout(heights = (c(1,3))) +
  plot_annotation(tag_levels = "a")
ggsave("figures/distribution_poptrends.png", width = 8.23, height = 9)

# plot estimated derivative trend ----
(plot_deriv = ggplot(data = derivs_pops_df) +
   geom_hline(yintercept = 0, col = "grey") +
   geom_ribbon(aes(ymin = lower, ymax = upper,
                   x = year, group = name, group = name), alpha = .1, fill = "#6497b1") +
   geom_line(aes(y = value, x = year, group = name), lwd = .2) +
   labs(x = "", 
        y = "Average rate of change in biomass") +
   theme(legend.position = "none") +
   coord_cartesian(ylim = c(-1.1, 1.1)))



# calculative cumulative sum of the derivative ----

derivs_pops_df$cumulative_deriv = NA
for(n in unique(derivs_pops_df$name)){
  temp = derivs_pops_df$value[which(derivs_pops_df$name == n)]
  cumulative_sum_derivatives = c(0)
  for(i in 2:length(temp)){
    cumulative_sum_derivatives[i] = sum(temp[1:i], na.rm = TRUE)
  }
  derivs_pops_df$cumulative_deriv[which(derivs_pops_df$name == n)] <- cumulative_sum_derivatives
}
