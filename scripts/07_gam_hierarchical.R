# Script to make a multivariate hierarchical generalized additive model 

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
matplot(biomass, type = "l")

YMeans = apply(Year_Geom_Means_all, 2, mean, na.rm = TRUE)
baselines = Year_Geom_Means_all[1,] |> as.matrix() |> as.vector()
weightings = YMeans/sum(YMeans)*100

# plot the contribution of each species to the community's abundance:
weightings_toplot = data.frame("sp" = stringr::str_to_sentence(gsub("_", " ", names(weightings))),
           "contribution" = baselines/sum(baselines))
weightings_toplot$sp = factor(weightings_toplot$sp,
                              levels = weightings_toplot$sp[order(weightings_toplot$contribution)])

ggplot(data = weightings_toplot) +
  geom_segment(aes(y = sp, yend = sp,
                   x = 0, xend = contribution), lwd = .5, col = "#6497b1") +
  geom_point(aes(y = sp, x = contribution), size = 3, col = "#03396c") +
  labs(x = "Contribution to the community's baseline biomass",
       y = "Species") +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_x_continuous(labels = scales::percent)
ggsave("figures/contributions_to_communitybiomass.png", width = 8, height = 5)

contributions_to_baseline = Year_Geom_Means_all[1,]/sum(Year_Geom_Means_all[1,])*100

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
                     trend_model = 'GP',
                     chains = 3,
                     use_stan = TRUE,
                     prior_simulation = TRUE)

# record the priors
test_priors <- get_mvgam_priors(y ~ 
                                  # global smoother for all pops over time
                                  s(time, bs = "tp", k = knots) + 
                                  # random intercept per group
                                  s(series, bs = 're', k = npops),
                                family = "gaussian",
                                data = data_train,
                                trend_model = 'GP',
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
              trend_model = 'GP',
              use_stan = TRUE,
              chains = 3, 
              burnin = 5000,
              samples = 10000
)
saveRDS(mod1, paste0("outputs/gam_hierarchical_gp.rds")) 
mod1 = readRDS(paste0("outputs/gam_hierarchical_gp.rds"))

# to view the Stan model file:
m <- mod1
# save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(mod1)

# convergence
png("figures/traceplots_timebases.png", height=1800, width=1800, type = "cairo")
rstan::stan_trace(mod1$model_output, 'b')
dev.off()

## lfo cross-validation just added to mvgam ----
#lfo_cv.mvgam()

## plot the temporal trend (time smooth) without population random effect----

png("figures/timesmooth.png", height=1000, width=1000, type = "cairo")
plot_mvgam_smooth(mod1, smooth = "time") # this is essentially the LPI curve
dev.off()

par(mfrow = c(3,3))
for(i in 1:npops){
  plot_mvgam_smooth(mod1, smooth = "time", series = i)
}

mvgam::plot_mvgam_randomeffects(mod1)

png("figures/mvgam_latentfactorplot.png", height=1000, width=1800, type = "cairo",pointsize = 20)
plot_mvgam_factors(mod1)
dev.off()

# plot(mod1)

# extract posterior draws in an array format
draws_fit = m$model_output |> posterior::as_draws_matrix() # mvgam_draws() now exists too
posterior_df = posterior::summarize_draws(draws_fit)
saveRDS(posterior_df, "outputs/gam_hierarchical_posterior.rds")

coefs_summ = coef(m, summarise = TRUE) |> as.data.frame()
saveRDS(coefs_summ, "outputs/gam_hierarchical_coefs_summary.rds")

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

# plot the predictions next to the data ----

data_train = rename(data_train, "sp" = "series")
data_train$sp = gsub("_", " ", data_train$sp) |> stringr::str_to_sentence()
data_train$sp = gsub(" ", "\n", data_train$sp)
predictions$sp = gsub("_", " ", predictions$sp) |> stringr::str_to_sentence()
predictions$sp = gsub(" ", "\n", predictions$sp)
ggplot(data = predictions,
       aes(x = time, group = sp)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = sp), alpha = .2) +
  geom_line(aes(y = biomass)) +
  geom_line(data = data_train,
            aes(x = time+1981, y = y), lty = 2) +
  facet_wrap(~sp, scales = "free_y") +
  theme(legend.position = "none",strip.text = element_text(face = "italic")) +
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

predictions_respscale_summary = predictions_respscale |>
  group_by(sp) |>
  summarise("mu" = mean(biomass),
            "sd" = sd(biomass))
predictions_respscale_summary$sp = factor(predictions_respscale_summary$sp,
                                          levels = predictions_respscale_summary$sp[order(predictions_respscale_summary$mu)])

predictions_respscale_perc = predictions |> group_by(sp) |> group_split()
for(i in 1:length(predictions_respscale_perc)){
# convert back to proportion of the baseline biomass
  predictions_respscale_perc[[i]]$biomass = ((predictions_respscale_perc[[i]]$biomass*YMeans[i]) + baselines[i])/baselines[i]
  predictions_respscale_perc[[i]]$cilo = ((predictions_respscale_perc[[i]]$cilo*YMeans[i]) + baselines[i])/baselines[i]
  predictions_respscale_perc[[i]]$cihi = ((predictions_respscale_perc[[i]]$cihi*YMeans[i]) + baselines[i])/baselines[i]
  
}
predictions_respscale_perc = bind_rows(predictions_respscale_perc)


(plot_biomassdiff = ggplot(data = predictions_respscale_perc) +
    geom_line(aes(y = biomass, x = time, group = sp)) +
    geom_hline(yintercept = 0, lwd = .3, lty = 2) +
    labs(x = "Population size relative to the baseline", 
         y = "",
         fill = "") +
    #scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-20,1)) +
    #coord_cartesian(xlim = c(-.5, .5)) +
    theme_pubr() +
    theme(axis.text.y = element_text(face = "italic"),
          panel.grid.major.x = element_line()) 
)


# take the derivative of each population trend ----
source("~/Documents/GitHub/hierarchical-lpi/scripts/plot_mvgam_trend_custom.R")
trend_vals = list()
preds_ls = list()
preds_respscale_ls = list()
preds_vsbaseline_ls = list() # proportion of baseline abundance
derivs_ls = list()
derivs_respscale_ls = list()
for(i in 1:npops){
  trend_vals[[i]] = plot_mvgam_trend_custom(mod1, derivatives = TRUE, series = i)
  preds_ls[[i]] = trend_vals[[i]]$preds
  preds_respscale_ls[[i]] = trend_vals[[i]]$preds*YMeans[i]
  preds_vsbaseline_ls[[i]] = (trend_vals[[i]]$preds*YMeans[i] + baselines[i])/baselines[i]
  derivs_ls[[i]] = trend_vals[[i]]$derivs
}
preds = do.call(rbind, preds_ls)
preds_respscale = do.call(rbind, preds_respscale_ls)
preds_vsbaseline = do.call(rbind, preds_vsbaseline_ls)
derivs = do.call(rbind, derivs_ls)

# each species' median predicted trend in biomass
preds_pops = do.call(cbind, lapply(preds_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(preds_pops, x = time+1981, type = "l")

# each species' median derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE)))
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))
matplot(derivs_pops, x = time+1981, type = "l")

## weighted mean???? ----

# plot each species' mean derivative over the whole time series ----

(temporal_trend = data.frame(
  species = colnames(Year_Geom_Means_all),
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

derivs_pops_df$name = gsub("_", " ", derivs_pops_df$name) |> 
  stringr::str_to_sentence()
derivs_pops_df$name =  factor(derivs_pops_df$name,
                                 levels = temporal_trend$species[order(temporal_trend$mu_deriv)])

# ridgeplot of derivatives per species ----
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

# ridgeplot of derivatives per year ----
derivs_pops_df$year_factor = factor(derivs_pops_df$year,
                                    levels = (as.character(time+1981)))
(plot_deriv_densities_annual = ggplot(data = derivs_pops_df) +
    ggridges::geom_density_ridges_gradient(aes(x = value, y = year_factor, 
                                               fill = after_stat(x)),
                                           quantile_lines = TRUE,
                                           quantile_fun = function(x,...) median(x),
                                           size = .4, scale = 2.5) +
    scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-.5,.5)) +
    geom_vline(xintercept = 0, lwd = .5, lty = 2) +
    theme(panel.grid.major.x = element_line(), legend.position = "none") +
    labs(x = "Annual rate of change (1981-2013)",
         y = "",
         fill = "")) +
  coord_cartesian(xlim = c(-.5,.5)) +
  coord_flip()
ggsave("figures/distribution_yearlytrends_ridgeplot.png", width = 4.03, height = 7.58)


derivs_pops_df$cumulative_deriv = NA
for(n in unique(derivs_pops_df$name)){
  temp = derivs_pops_df$value[which(derivs_pops_df$name == n)]
  cumulative_sum_derivatives = c(0)
  for(i in 2:length(temp)){
    cumulative_sum_derivatives[i] = sum(temp[1:i], na.rm = TRUE)
  }
  derivs_pops_df$cumulative_deriv[which(derivs_pops_df$name == n)] <- cumulative_sum_derivatives
}

# plot a histogram of the average derivatives across all species ----
derivs_without1981 = dplyr::filter(derivs_pops_df, year != 1981)
(plot_trenddensity = 
   ggplot(data = derivs_without1981) +
   geom_histogram(aes(x = value, fill = after_stat(x)), 
                  col = "black", linewidth = .2, bins = 19) + 
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE)) +
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE) - sd(derivs_without1981$value, na.rm = TRUE), lty = 2) +
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE) + sd(derivs_without1981$value, na.rm = TRUE), lty = 2) +
   theme(panel.grid.major.x = element_line()) +
   scale_y_sqrt() +
   labs(x = "Annual rate of change", 
        y = "Frequency", 
        fill = "") +
   scale_fill_distiller(palette = "RdYlGn", 
                        direction = 1, 
                        limits = c(-.4,.4)) +
   coord_cartesian(xlim = c(-.4, .4))) 
# variance
var(derivs_without1981$value, na.rm = TRUE)

## animated per year ----
library(gganimate)
derivs_without1981$frame = as.character(derivs_without1981$year)
derivs_without1981_summary = derivs_without1981 |>
  group_by(year) |>
  summarise(mu = mean(value),
            lower = mean(value) - sd(value),
            upper = mean(value) + sd(value))
derivs_without1981_summary$frame = as.character(derivs_without1981_summary$year)

(plot_trenddensity_anim = 
   ggplot(data = derivs_without1981) +
   geom_histogram(aes(x = value, fill = after_stat(x)), 
                  col = "black", linewidth = .2, bins = 19) + 
    geom_vline(data = derivs_without1981_summary, aes(xintercept = mu)) +
    geom_vline(data = derivs_without1981_summary, aes(xintercept = lower), lty = 2) +
    geom_vline(data = derivs_without1981_summary, aes(xintercept = upper), lty = 2) +
   theme(panel.grid.major.x = element_line(),
         axis.text = element_text(size = 16),
         legend.text = element_text(size = 16),
         title = element_text(size = 20),
         legend.position = "none") +
   scale_y_sqrt() +
   labs(x = "Taux de croissance annuel", 
        y = "Fréquence", 
        fill = "",
        title = "Année: {closest_state}") +
   scale_fill_distiller(palette = "RdYlGn", 
                        direction = 1, 
                        limits = c(-.4,.4)) +
   coord_cartesian(xlim = c(-.4, .4)) +
   transition_states(frame, state_length = 50)) 
anim_save("figures/trend_density.gif", plot_trenddensity_anim)

#----

# plot the average derivative with 90% quantiles per species ----

(plot_poptrends = ggplot(data = temporal_trend[order(temporal_trend$mu_deriv),]) +
   geom_segment(aes(x = lower, xend = upper,
                    y = species, yend = species), lwd = .3) +
   geom_point(aes(x = mu_deriv, y = species, fill = mu_deriv), 
              size = 4, pch = 21) +
   geom_vline(xintercept = 0, lwd = .3, lty = 2) +
   labs(x = "Median rate of change", 
        y = "",
        fill = "") +
   scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-.1,.1)) +
   coord_cartesian(xlim = c(-.4, .4)) +
   theme_pubr() +
   theme(axis.text.y = element_text(face = "italic"),
         panel.grid.major.x = element_line()) 
)

# avg biomass trend ----

avg_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds, 2, mean),
  "cilo" = apply(preds, 2, quantile, prob = .05),
  "cihi" = apply(preds, 2, quantile, prob = .95)
)
saveRDS(avg_trend, "outputs/gam_hierarchical_df_overall.rds")

# weighted avg biomass trend ----

avg_weighted_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds, 2, Hmisc::wtd.mean, weights = rep(weightings, each = 30000), normwt = TRUE),
  "cilo" = apply(preds, 2, Hmisc::wtd.quantile, prob = .05, weights = rep(weightings, each = 30000), normwt = TRUE),
  "cihi" = apply(preds, 2, Hmisc::wtd.quantile, prob = .95, weights = rep(weightings, each = 30000), normwt = TRUE)
)
saveRDS(avg_weighted_trend, "outputs/gam_hierarchical_df_overall_weighted.rds")

# avg biomass trend on response scale ----

avg_trend_respscale = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds_respscale, 2, mean),
  "cilo" = apply(preds_respscale, 2, quantile, prob = .05),
  "cihi" = apply(preds_respscale, 2, quantile, prob = .95)
)
saveRDS(avg_trend_respscale, "outputs/gam_hierarchical_df_overall_respscale.rds")

# weighted avg biomass trend ----

avg_weighted_trend_respscale = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds_respscale, 2, Hmisc::wtd.mean, weights = rep(weightings, each = 30000), normwt = TRUE),
  "cilo" = apply(preds_respscale, 2, Hmisc::wtd.quantile, prob = .05, weights = rep(weightings, each = 30000), normwt = TRUE),
  "cihi" = apply(preds_respscale, 2, Hmisc::wtd.quantile, prob = .95, weights = rep(weightings, each = 30000), normwt = TRUE)
)
saveRDS(avg_weighted_trend_respscale, "outputs/gam_hierarchical_df_overall_weighted_respscale.rds")


# avg biomass trend as proportion of baseline biomass ====
avg_trend_vsbaseline = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds_vsbaseline, 2, mean, na.rm = TRUE),
  "cilo" = apply(preds_vsbaseline, 2, quantile, prob = .05),
  "cihi" = apply(preds_vsbaseline, 2, quantile, prob = .95)
)
saveRDS(avg_trend_respscale, "outputs/gam_hierarchical_df_overall_preds_vsbaseline.rds")


# avg derivative ----

avg_deriv_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs, 2, mean, na.rm = T),
  "cilo" = apply(derivs, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs, 2, quantile, prob = .95, na.rm = T)
  )
avg_deriv_trend[1,2:4] = 0
saveRDS(avg_deriv_trend, "outputs/gam_hierarchical_df_overall_deriv.rds")

# weighted avg derivative ----
# matrix_weightings = rep(weightings, each = 30000)
# matrix_weightings[which(is.na(derivs))] <- 0
# avg_weighted_deriv_trend = data.frame(
#   "year" = time+1981,
#   "avg_trend" = apply(derivs, 2, Hmisc::wtd.mean, weights = rep(weightings, each = 30000), normwt = TRUE, na.rm = TRUE)#,
#  # "cilo" = apply(derivs, 2, Hmisc::wtd.quantile, prob = .05, weights = matrix_weightings, normwt = TRUE, na.rm = TRUE)#,
#  # "cihi" = apply(derivs, 2, Hmisc::wtd.quantile, prob = .95, weights = rep(weightings, each = 30000), normwt = TRUE, na.rm = TRUE)
# )
# avg_weighted_deriv_trend[1,2:4] = 0
# saveRDS(avg_weighted_deriv_trend, "outputs/gam_hierarchical_df_overall_deriv_weighted.rds")

# average cumulative change trend ----
# avg_cumu_deriv_trend = derivs_pops_df |>
#   group_by(year) |>
#   summarise(avg_trend = median(cumulative_deriv),
#             cilo = quantile(cumulative_deriv, probs = .05),
#             cihi = quantile(cumulative_deriv, probs = .95)) 

# cumsum_deriv <- function(dt){
#   cumulative_dt = 0
#   for(i in 2:tsl){
#     cumulative_dt[i] = sum(dt[1:i], na.rm = TRUE)
#   }
#   return(cumulative_dt)
# }
# 
# # calculate cumulative derivative for all permutations
# 
# avg_cumu_deriv_trend = data.frame(
#   "year" = avg_deriv_trend$year,
#   "avg_trend" = cumsum_deriv(avg_deriv_trend$avg_trend),
#   "lower" = cumsum_deriv(avg_deriv_trend$cilo),
#   "upper" = cumsum_deriv(avg_deriv_trend$cihi)
# )
# saveRDS(avg_cumu_deriv_trend, "outputs/gam_hierarchical_df_overall_cumuderiv.rds")

# plot average biomass trend
(A = ggplot(data = avg_trend, aes(x = year)) +
    geom_hline(yintercept = 0, lwd = .3) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = .8) +
    labs(x = "",
         y = "Average difference \nfrom baseline biomass \n(kg per tow)") +
    coord_cartesian(ylim = c(-5,5)) +
    theme(panel.grid.major.x = element_line()))

# plot avg derivative trend
(B = ggplot(data = avg_deriv_trend, aes(x = year)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = .8) +
    geom_hline(yintercept = 0, lwd = .3) +
    labs(x = "",
         y = "Average rate of change \n(kg per tow / year)") +
    coord_cartesian(ylim = c(-.5, .5)) +
    theme(panel.grid.major.x = element_line()))

(C = ggplot(data = avg_weighted_trend, aes(x = year)) +
    geom_hline(yintercept = 0, lwd = .3) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = .8) +
    labs(x = "",
         y = "Average difference \nfrom baseline biomass \n(kg per tow)") +
    coord_cartesian(ylim = c(-5,5)) +
    theme(panel.grid.major.x = element_line()))


A + B + plot_annotation(tag_levels = "a")
ggsave("figures/mvgam_prediction.png", width = 8.46, height = 4.06)


# FIGURES ######################################################################


# fig 1. data and model predictions of biomass ----

(fig1a = ggplot(data = dat) +
   geom_hline(yintercept = 0, lwd = .3, col = "grey") +
  geom_line(aes(x = time+1981,
            y = y*YMeans,
            group = series), lwd = .4, col = "#03396c") +
   theme(legend.position = "none") +
  labs(x = "Year", y = "Difference in biomass (kg per tow)", col = "Species") +
   coord_cartesian(ylim = c(-30,30)))

(fig1b = ggplot(data = predictions_respscale,
                aes(x = time, group = sp)) +
    geom_hline(yintercept = 0, lwd = .3, col = "grey") +
    geom_ribbon(aes(ymin = cilo, ymax = cihi, group = sp), alpha = .2, fill = "#6497b1") +
    geom_line(aes(y = biomass), lwd = .4, col = "#03396c") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(-30,10)) +
    labs(x = "Year", y = "Predicted difference \nfrom baseline biomass \n(kg per tow)", fill = "Species"))

(fig1c <- ggplot(data = avg_trend_respscale, aes(x = year)) +
    geom_hline(yintercept = 0, lwd = .3, col = "grey") +
    geom_ribbon(aes(ymin = cilo, ymax = cihi),
              alpha = .3, fill = "#6497b1") +
  geom_line(aes(y = avg_trend), col = "#03396c", lwd = .4) +
  labs(x = "Year",
       y = "Average predicted difference \nfrom baseline biomass \n(kg per tow)") +
  coord_cartesian(ylim = c(-10,10)) +
  theme(panel.grid.major.x = element_line()))

(avg_derivative_plot = ggplot(data = avg_deriv_trend, aes(x = year)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c") +
    geom_hline(yintercept = 0, lwd = .3) +
    labs(x = "",
         y = "Average annual rate of change") +
    coord_cartesian(ylim = c(-.5, .5)) +
    theme(panel.grid.major.x = element_line()) +
    scale_y_continuous(labels = scales::percent))


#(fig1a + fig1b + fig1c) + plot_annotation(tag_levels = "a")
fig1a + fig1c + avg_derivative_plot + plot_annotation(tag_levels = "a")
ggsave("figures/fig1.png", width = 11, height = 4.21)

# fig 2 - distribution of species' rates of change ----

# arrange the plot panels and save
(plot_trenddensity + 
    coord_cartesian(xlim = c(-.5, .5)) +
    theme(axis.title.y = element_text(vjust = -60)) #+
) / (plot_poptrends + 
       coord_cartesian(xlim = c(-.5, .5)) +
       theme(axis.text.y = element_text(face = "italic"),
             panel.grid.major.x = element_line())) + 
  plot_layout(heights = (c(1,3))) +
  plot_annotation(tag_levels = "a")
ggsave("figures/distribution_poptrends.png", width = 8, height = 9)


# fig 3 - species correlations ----

# get species correlations
data_train = dat
sp_correlations = lv_correlations(mod1)
saveRDS(sp_correlations, "outputs/gam_hierarchical_species_correlations.rds")

# sp_correlations = readRDS("outputs/gam_hierarchical_species_correlations.rds")

# edit the species names to be prettier
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# Plot as heatmap
png(height=1800, width=1800, file="figures/fig3_species_associations.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

png(height=500, width=500, file="figures/fig3b_species_associations_histogram.png", type = "cairo")
sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))] |>
  hist(col = "grey20", border = "white", lwd = .2,
       xlab = "Corrélation entre espèces",
       ylab = "Fréquence",
       main = "", cex = 3)
dev.off()



corrs = sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))]

corrs |> quantile(probs = c(0.05, 0.5, 0.95))


length(corrs[which(corrs > .75)]) / length(corrs)
length(corrs[which(corrs < -.75)]) / length(corrs)
temporal_trend$mu_deriv

# fig 4 - summary of population change with indices ----

## Plot rate of change trend versus the LPI for comparison 

rlpi_results <- readRDS("outputs/rlpi_results.rds")
rlpi_results$time <- as.integer(rownames(rlpi_results))
rlpi_results <- rlpi_results[-nrow(rlpi_results),] # remove last value

# import the population trends
rlpi_poptrends <- read.csv("default_infile_pops_lambda.csv", row.names = 1) |>
  subset(select = -c(Freq)) |>
  rename("species" = "SpeciesSSet")
rlpi_poptrends$X1981 = 0
# create long version for plotting
rlpi_poptrends <- pivot_longer(rlpi_poptrends, cols = c(-species)) |>
  rename("time" = "name")
rlpi_poptrends$time = readr::parse_number(rlpi_poptrends$time)

rlpi_avgtrend <- read.csv("default_infile_pops_dtemp.csv", row.names = 1) |>
  pivot_longer(cols = everything()) |>
  rename("time" = "name")
rlpi_avgtrend$time = readr::parse_number(rlpi_avgtrend$time)
rlpi_avgtrend = full_join(data.frame(time = 1980, value = 0), rlpi_avgtrend)

(rlpi_index = ggplot() +
    geom_ribbon(data = rlpi_results, 
                aes(ymin = CI_low, ymax = CI_high, x = time), alpha = .3, fill = "#6497b1") + 
    geom_line(data = rlpi_results, aes(x = time, y = LPI_final), col = "#03396c") +
    geom_hline(yintercept = 1, lwd = .3) +
    labs(y = "Living Planet Index", x = "") +
    coord_cartesian(ylim = c(0, 2))+
    theme(panel.grid.major.x = element_line()))

(avg_derivative_plot = ggplot(data = avg_deriv_trend, aes(x = year)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c") +
    geom_hline(yintercept = 0, lwd = .3) +
    labs(x = "",
         y = "Average annual rate of change") +
    coord_cartesian(ylim = c(-.5, .5)) +
    theme(panel.grid.major.x = element_line()) +
    scale_y_continuous(labels = scales::percent))

avg_derivative_plot + rlpi_index + plot_annotation(tag_levels = "a")
rlpi_index
ggsave("figures/compare_lpi_mvgam.png", width = 6, height = 4)

## calculate an index of change ----

# # make a function to calculate the change in a population from a baseline value
# calc_index = function(derivative, baseline = 1){
#   index = baseline
#   for(i in 2:tsl){
#     index[i] = index[i-1]*(1 + derivative[i])
#   }
#   return(index)
# }
# 
# #### draw 100 random samples from each species' posterior to make a credible interval
# # for the average derivative trend
# index_allperms_perseries = lapply(derivs_ls, function(x) apply(x, 1, FUN = calc_index))
# # # function(x) apply(x[runif(100, 1, nrow(x)),], 1, FUN = calc_index))
# # index_perspecies = lapply(index_allperms_perseries,
# # function(x) apply(x, 1, quantile, probs = c(.05, .5, .95)))
# # names(index_perspecies) = colnames(YData)
# # index_perspecies = lapply(index_perspecies, function(x) as.data.frame(t(x)))
# # index_perspecies = bind_rows(index_perspecies, .id = "species")
# # colnames(index_perspecies) = c("species", "q05", "q50", "q95")
# # index_perspecies$year = rep(time+1981, npops)
# 
# # temp = index_perspecies |>
# #   group_by(year) |>
# #   summarise(mean_index = mean(q50),
# #             sd_index = sd(q50))
# 
# # draw 100 random samples per species
# randoms = lapply(index_allperms_perseries, function(x) x[,runif(100, 1, nrow(x))])
# index_allperms = bind_cols(randoms)
# index_allperms_summary = data.frame(
#   #"index" = temp$mean_index,
#   # "lower_sd" = temp$mean_index - temp$sd_index,
#   # "upper_sd" = temp$mean_index + temp$sd_index,
#   "index" = apply(index_allperms, 1, quantile, probs = c(.5), na.rm = TRUE),
#   "lower" = apply(index_allperms, 1, quantile, probs = c(.05), na.rm = TRUE),
#   "upper" = apply(index_allperms, 1, quantile, probs = c(.95), na.rm = TRUE)
# )
# index_allperms_summary$year = time+1981
# 
# 
# calc_index(avg_deriv_trend$avg_trend)
# calc_index(avg_deriv_trend$cilo)
# calc_index(avg_deriv_trend$cihi)
# 
# (model_index = ggplot(data = index_allperms_summary, aes(x = year)) +
#   geom_hline(yintercept = 1, lty = 2) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3, fill = "#6497b1") +
#   geom_line(aes(y = index)) +
#   labs(y = "Index of change in population size", x = "")) + coord_cartesian(ylim = c(0,2))
# ggsave("figures/index_popsize.png", width = 5.72, height = 4.82)

# plot the derivatives per species with highlighted main 4 commercially fished species ----

ggplot(data = derivs_pops_df, 
       aes(x = year, group = name)) +
  geom_hline(yintercept = 0, col= "grey90") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = name), 
              alpha = .3, lwd = 0) +
  geom_line(aes(y = value, col = name), lwd = .6) +
  labs(x = "",
       fill = "Species", color = "Species",
       y = "Rate of change in biomass") +
  gghighlight::gghighlight(name %in% c("Gadus morhua",
                                       "Hippoglossoides platessoides",
                                       "Reinhardtius hippoglossoides",
                                       "Sebastes mentella"),
                           unhighlighted_colour = "grey85",
                           unhighlighted_params = list(fill = "grey90"),
                           use_direct_label = FALSE) +
  theme(legend.position = "right") +
  #viridis::scale_fill_viridis(discrete = TRUE) +
  #viridis::scale_color_viridis(discrete = TRUE, end = .95)
  ggsci::scale_fill_nejm() + ggsci::scale_color_nejm() +
  scale_y_continuous(labels = scales::percent)
ggsave("figures/species_derivative_trends.png", width = 8.33, height = 4.25)
