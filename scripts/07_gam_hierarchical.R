# Script to make a multivariate hierarchical generalized additive model 
# (i.e. with a year effect)

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
# this model is like GS
# assumption: there is a global trend, but populations can vary in their response to time.
mvgam_prior <- mvgam(data = data_train,
                     formula = y ~ 
                       # global smoother for all pops over time
                       s(time, bs = "tp", k = knots) + 
                       # independent smoothers for each group
                       #s(time, by = series, m = 1, bs = "tp") +
                       # random intercept per group
                       s(series, bs = 're', k = npops),
                     family = "gaussian",
                     trend_model = 'GP',
                     chains = 3,
                     use_stan = TRUE,
                     prior_simulation = TRUE)

mvgam_prior$model_file

# record the priors
test_priors <- get_mvgam_priors(y ~ 
                                  # global smoother for all pops over time
                                  s(time, bs = "tp", k = knots) + 
                                  # # independent smoothers for each group
                                  # s(time, by = series, bs = "tp", k = knots) +
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

# nvm. the number of knots seems to be important to choose here, bc too high = the model has too many things to estimate?
# train the model on data ----
mod1 <- mvgam(data = data_train,
              #formula = y ~ s(time, bs = "tp", k = knots) + s(series, bs = "re"),
              formula =  y ~ s(time, bs = "tp", k = knots) + 
                # s(time, series, bs = "fs", k = knots, m = 2) + 
                s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'GP',
              use_stan = TRUE,
              chains = 3, 
              burnin = 5000,
              samples = 10000
)
saveRDS(mod1, paste0("outputs/gam_hierarchical.rds")) 
mod1 = readRDS(paste0("outputs/gam_hierarchical.rds"))

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

# assign the population name to the numbers in the parameter names
# extract the model summary
summary_df = posterior_df[-c(grep("ypred", posterior_df$variable),
                             grep("mu_raw", posterior_df$variable),
                             grep("sigma_raw", posterior_df$variable)),]
pops = data.frame(
  "pop" = colnames(biomass)[-31],
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


# # mean trend from the model (without obs error)
# posterior_mus = posterior_df[grep("mus", posterior_df$variable),] |>
#   separate(col = variable, into = c("year", "sp"), sep = ",") |>
#   mutate(year = readr::parse_number(year),
#          sp = readr::parse_number(sp))

# take the derivative of each population trend
source("~/Documents/GitHub/hierarchical-lpi/scripts/plot_mvgam_trend_custom.R")
trend_vals = list()
preds_ls = list()
preds_respscale_ls = list()
derivs_ls = list()
derivs_respscale_ls = list()
for(i in 1:npops){
  trend_vals[[i]] = plot_mvgam_trend_custom(mod1, derivatives = TRUE, series = i)
  preds_ls[[i]] = trend_vals[[i]]$preds
  preds_respscale_ls[[i]] = trend_vals[[i]]$preds*YMeans[i]
  derivs_ls[[i]] = trend_vals[[i]]$derivs
  derivs_respscale_ls[[i]] = trend_vals[[i]]$derivs*YMeans[i]
}
preds = do.call(rbind, preds_ls)
preds_respscale = do.call(rbind, preds_respscale_ls)
derivs = do.call(rbind, derivs_ls)
derivs_respscale = do.call(rbind, derivs_respscale_ls)

# each species' median predicted trend in biomass
preds_pops = do.call(cbind, lapply(preds_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(preds_pops, x = time+1981, type = "l")

preds_respscale_pops = do.call(cbind, lapply(preds_respscale_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(preds_respscale_pops, x = time+1981, type = "l")

derivs_respscale_pops = do.call(cbind, lapply(derivs_respscale_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
matplot(derivs_respscale_pops, x = time+1981, type = "l")

# each species' median derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE)))
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))

matplot(derivs_pops, x = time+1981, type = "l")

# plot each species' mean derivative over the whole time series ----

(temporal_trend = data.frame(
  species = colnames(Year_Geom_Means_all),
  mu_deriv = derivs_pops |> apply(2, mean, na.rm = TRUE),
  #sd_deriv = derivs_pops |> apply(2, sd, na.rm = TRUE)
  lower = derivs_pops |> apply(2, quantile, probs = .05, na.rm = TRUE),
  upper = derivs_pops |> apply(2, quantile, probs = .95, na.rm = TRUE)
))
# temporal_trend$lower = temporal_trend$mu_deriv - temporal_trend$sd_deriv
# temporal_trend$upper = temporal_trend$mu_deriv + temporal_trend$sd_deriv
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
# temp_sd = derivs_pops_sd |> 
#   as.data.frame() |> 
#   mutate(year = time + 1981) |>
#   pivot_longer(cols = -year, values_to = "sd")
derivs_pops_df = full_join(derivs_pops_df, derivs_pops_df_lower) |>
  full_join(derivs_pops_df_upper)
# set first time step's derivative to 0
derivs_pops_df$value[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$lower[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$upper[which(derivs_pops_df$year == 1981)] <- 0

#derivs_pops_df$sd[which(derivs_pops_df$year == 1981)] <- 0

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
                                    levels = rev(as.character(time+1981)))
(plot_deriv_densities_annual = ggplot(data = derivs_pops_df) +
    ggridges::geom_density_ridges_gradient(aes(x = value, y = year_factor, 
                                               fill = after_stat(x)),
                                           quantile_lines = TRUE,
                                           quantile_fun = function(x,...) median(x),
                                           size = .4, scale = 2.5) +
    # scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", limits = c(-.5,.5)) +
    scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-.5,.5)) +
    geom_vline(xintercept = 0, lwd = .5, lty = 2) +
    theme(panel.grid.major.x = element_line(), legend.position = "none") +
    labs(x = "Annual rate of change (1981-2013)",
         y = "",
         fill = "")) +
  coord_cartesian(xlim = c(-.5,.5))
ggsave("figures/distribution_yearlytrends_ridgeplot.png", width = 4.03, height = 7.58)


derivs_pops_df$cumulative_deriv = NA
#derivs_pops_df$cumulative_deriv_sd = NA
for(n in unique(derivs_pops_df$name)){
  temp = derivs_pops_df$value[which(derivs_pops_df$name == n)]
  #temp2 = derivs_pops_df$sd[which(derivs_pops_df$name == n)]
  
  cumulative_sum_derivatives = c(0)
  #cumulative_sum_derivatives_sd = c(0)
  for(i in 2:length(temp)){
    cumulative_sum_derivatives[i] = sum(temp[1:i], na.rm = TRUE)
    #cumulative_sum_derivatives_sd[i] = sum(temp2[1:i], na.rm = TRUE)
  }
  
  derivs_pops_df$cumulative_deriv[which(derivs_pops_df$name == n)] <- cumulative_sum_derivatives
  #derivs_pops_df$cumulative_deriv_sd[which(derivs_pops_df$name == n)] <- cumulative_sum_derivatives_sd
}


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
                        limits = c(-.5,.5)) +
   coord_cartesian(xlim = c(-.5, .5))) 

# plot the average derivative with 90% quantiles per species ----

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

# plot the average derivative with 90% quantiles per year ----

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

# plot estimated derivative trend ----
ggplot(data = derivs_respscale) +
  geom_line(aes(y = value, x = year, group = name), lwd = .2) +
  coord_cartesian(ylim = c(-3,1)) +
  labs(x = "",
       y = "Estimated rate of change in biomass") +
  theme(legend.position = "right") +
  gghighlight::gghighlight(name %in% c("Gadus morhua",
                                       "Hippoglossoides platessoides",
                                       "Reinhardtius hippoglossoides",
                                       "Sebastes mentella"),
                           use_direct_label = FALSE
  )

(plot_deriv = ggplot(data = derivs_pops_df) +
   geom_hline(yintercept = 0, col = "grey") +
    geom_ribbon(aes(ymin = lower, ymax = upper,
                    x = year, group = name, group = name), alpha = .1, fill = "#6497b1") +
    geom_line(aes(y = value, x = year, group = name), lwd = .2) +
    labs(x = "", 
         y = "Average rate of change in biomass") +
    theme(legend.position = "none") +
    # gghighlight::gghighlight(name %in% c("Gadus morhua", 
    #                                      "Hippoglossoides platessoides",
    #                                      "Reinhardtius hippoglossoides",
    #                                      "Sebastes mentella"),
    #                          use_direct_label = FALSE,
    #                          unhighlighted_params = list(colour = NULL)) +
    coord_cartesian(ylim = c(-.6, .6)))

# # plot cumulative derivative trend ----
# (plot_cumu_deriv = ggplot(data = derivs_pops_df) +
#    geom_hline(yintercept = 0, col = "grey") +
#     # geom_ribbon(aes(ymin = cumulative_deriv - cumulative_deriv_sd,
#     #                 ymax = cumulative_deriv + cumulative_deriv_sd,
#     #                 x = year, group = name), alpha = .1, fill = "#6497b1") +
#   geom_line(aes(y = cumulative_deriv, x = year, group = name), lwd = .2) +
#     labs(x = "", y = "Cumulative change in biomass", fill = "Species") +
#     theme(legend.position = "right",
#           legend.text = element_text(face = "italic")) +
#   # gghighlight::gghighlight(name %in% c("Gadus morhua", 
#   #                                      "Hippoglossoides platessoides",
#   #                                      "Reinhardtius hippoglossoides",
#   #                                      "Sebastes mentella"),
#   #                          use_direct_label = FALSE,
#   #                          unhighlighted_params = list(colour = NULL)) +
#   coord_cartesian(ylim = c(-6, 6)))

#|> plotly::ggplotly() 

# plot_deriv + plot_cumu_deriv + plot_annotation(tag_levels = "a")
# ggsave("figures/mvgam_rateofchange.png", width = 11, height = 4.06)

# avg biomass trend ----

avg_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds, 2, mean),
  "cilo" = apply(preds, 2, quantile, prob = .05),
  "cihi" = apply(preds, 2, quantile, prob = .95)
)
saveRDS(avg_trend, "outputs/gam_hierarchical_df_overall.rds")

# avg biomass trend on response scale ----

avg_trend_respscale = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds_respscale, 2, mean),
  "cilo" = apply(preds_respscale, 2, quantile, prob = .05),
  "cihi" = apply(preds_respscale, 2, quantile, prob = .95)
)
saveRDS(avg_trend_respscale, "outputs/gam_hierarchical_df_overall_respscale.rds")


# avg derivative ----

avg_deriv_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs, 2, mean, na.rm = T),
  "cilo" = apply(derivs, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs, 2, quantile, prob = .95, na.rm = T)
  )
avg_deriv_trend[1,2:4] = 0
saveRDS(avg_deriv_trend, "outputs/gam_hierarchical_df_overall_deriv.rds")

avg_deriv_respscale_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs_respscale, 2, mean, na.rm = T),
  "cilo" = apply(derivs_respscale, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs_respscale, 2, quantile, prob = .95, na.rm = T)
)
avg_deriv_respscale_trend[1,2:4] = 0
#saveRDS(avg_deriv_trend, "outputs/gam_hierarchical_df_overall_deriv.rds")


# average cumulative change trend ----
avg_cumu_deriv_trend = derivs_pops_df |>
  group_by(year) |>
  summarise(avg_trend = mean(cumulative_deriv),
            cilo = quantile(cumulative_deriv, probs = .05),
            cihi = quantile(cumulative_deriv, probs = .95)) 
saveRDS(avg_cumu_deriv_trend, "outputs/gam_hierarchical_df_overall_cumuderiv.rds")

# plot average biomass trend
(A = ggplot(data = avg_trend, aes(x = year)) +
    # geom_line(data = predictions, aes(x = time, y = biomass, group = sp),
    #           linewidth = .2, col = "grey") +
    geom_hline(yintercept = 0, lwd = .3) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = .8) +
    labs(x = "",
         y = "Average difference \nfrom baseline biomass \n(kg per tow)") +
    coord_cartesian(ylim = c(-5,5)) +
    theme(panel.grid.major.x = element_line()))

# plot average biomass trend
(A2 = ggplot(data = avg_trend_respscale, aes(x = year)) +
    # geom_line(data = predictions, aes(x = time, y = biomass, group = sp),
    #           linewidth = .2, col = "grey") +
    geom_hline(yintercept = 0, lwd = .3) +
    # geom_line(data = predictions_respscale,
    #           aes(x = time, y = biomass, group = sp), lwd = .1, col = "#6497b1") +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = .8) +
    coord_cartesian(ylim = c(-10,10)) +
    labs(x = "",
         y = "Average difference \nfrom baseline biomass \n(kg per tow)") +
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

(B2 = ggplot(data = avg_deriv_respscale_trend, aes(x = year)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = .8) +
    geom_hline(yintercept = 0, lwd = .3) +
    labs(x = "",
         y = "Average rate of change \n(kg per tow / year)") +
    #coord_cartesian(ylim = c(-.5, .5)) +
    theme(panel.grid.major.x = element_line()))

A + B + plot_annotation(tag_levels = "a")
A2 + B2 + plot_annotation(tag_levels = "a") #---- key figure
ggsave("figures/mvgam_prediction.png", width = 8.46, height = 4.06)


# FIGURES ######################################################################


# fig 1. data and model predictions of biomass ----

(fig1a = ggplot(data = data_train) +
   geom_hline(yintercept = 0, lwd = .3, col = "grey") +
  geom_line(aes(x = time+1981,
            y = y,
            group = sp), lwd = .4, col = "#03396c") +
   theme(legend.position = "none") +
  labs(x = "Year", y = "Difference in biomass (kg per tow)", col = "Species") +
   coord_cartesian(ylim = c(-30,5)))

(fig1b = ggplot(data = predictions_respscale,
                aes(x = time, group = sp)) +
    geom_hline(yintercept = 0, lwd = .3, col = "grey") +
    geom_ribbon(aes(ymin = cilo, ymax = cihi, group = sp), alpha = .2, fill = "#6497b1") +
    geom_line(aes(y = biomass), lwd = .4, col = "#03396c") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(-30,5)) +
    labs(x = "Year", y = "Predicted difference \nfrom baseline biomass \n(kg per tow)", fill = "Species"))

(fig1c <- ggplot(data = avg_trend_respscale, aes(x = year)) +
    geom_hline(yintercept = 0, lwd = .3, col = "grey") +
    geom_ribbon(aes(ymin = cilo, ymax = cihi),
              alpha = .3, fill = "#6497b1") +
  geom_line(aes(y = avg_trend), col = "#03396c", lwd = .4) +
  labs(x = "Year",
       y = "Average predicted difference \nfrom baseline biomass \n(kg per tow)") +
  coord_cartesian(ylim = c(-10,2)) +
  theme(panel.grid.major.x = element_line()))

# plot the trend for an average population ----

newdf <- data.frame(time = time, series = rep(0, tsl)) # zero out ranef terms
prediction_typicalpop = predict(mod1, newdata = newdf)
prediction_typicalpop_df = data.frame(
  time = time+1981,
  mu = apply(prediction_typicalpop, 2, mean),
  cilo = apply(prediction_typicalpop, 2, quantile, probs = .05),
  cihi = apply(prediction_typicalpop, 2, quantile, probs = .95)
)
#(fig1c <- 
  (ggplot(prediction_typicalpop_df) +
  geom_hline(yintercept = 0, lwd = .3, col = "grey80") +
  geom_ribbon(aes(x = time, ymin = cilo, ymax = cihi), alpha = .3, fill = "#6497b1") +
  geom_line(aes(x = time, y = mu), col = "#03396c") +
  labs(y = "Mean predicted biomass difference", x = "") +
  coord_cartesian(ylim = c(-5,5)))

(fig1a + fig1b + fig1c) + plot_annotation(tag_levels = "a")
ggsave("figures/fig1.png", width = 11, height = 4.21)

# fig 2 - distribution of species' rates of change ----

# These plots are coded above
# arrange the plot panels and save
(plot_trenddensity + 
    coord_cartesian(xlim = c(-.5, .5)) +
    theme(axis.title.y = element_text(vjust = -60)) #+
  # annotate("text", x = coefs$meansParamX[2]-0.4, y = 1.8, 
  #          label = paste0("\u03bc = ", round(coefs$meansParamX[2], digits = 3),
  #                         "\n \u03c3\u00b2 = ", round(coefs$varX[2,2], digits = 3)))
) / (plot_poptrends + 
       coord_cartesian(xlim = c(-.25, .25)) +
       theme(axis.text.y = element_text(face = "italic"),
             panel.grid.major.x = element_line())) + 
       #scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1,1))) +
  plot_layout(heights = (c(1,3))) +
  plot_annotation(tag_levels = "a")
ggsave("figures/distribution_poptrends.png", width = 8.23, height = 9)


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


# fig 4 - summary of population change with indices ----

## Plot cumulative change trend versus the LPI for comparison 

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

(avg_derivative_plot = ggplot(data = avg_deriv_respscale_trend, aes(x = year)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c") +
    geom_hline(yintercept = 0, lwd = .3) +
    labs(x = "",
         y = "Average rate of change \n(kg per tow / year)") +
    coord_cartesian(ylim = c(-.8, .8)) +
    theme(panel.grid.major.x = element_line()))

avg_derivative_plot + rlpi_index + plot_annotation(tag_levels = "a")
ggsave("figures/compare_lpi_mvgam.png", width = 8.68, height = 4.16)


# (B = ggplot(data = avg_cumu_deriv_trend, aes(x = year)) +
#     geom_ribbon(aes(ymin = avg_trend - sd, ymax = avg_trend + sd), 
#                 alpha = .3, fill = "#6497b1") +
#     geom_line(aes(y = avg_trend), col = "#03396c") +
#     geom_hline(yintercept = 0, lwd = .3) +
#     labs(x = "",
#          y = "Cumulative change in biomass\n(kg per tow)") +
#     coord_cartesian(ylim = c(-3, 3)) +
#     theme(panel.grid.major.x = element_line()))

# (C = ggplot(data = avg_deriv_trend, aes(x = year)) +
#     geom_ribbon(aes(ymin = cilo, ymax = cihi), 
#                 alpha = .3, fill = "#6497b1") +
#     geom_line(aes(y = avg_trend), col = "#03396c") +
#     geom_hline(yintercept = 0, lwd = .3) +
#     labs(x = "",
#          y = "Average rate of change\n(kg per tow / year)") +
#     coord_cartesian(ylim = c(-.32, .32)) +
#     theme(panel.grid.major.x = element_line()))

#C / B / A + plot_annotation(tag_levels = "a")

#B + A


# avg_deriv_trend$lpi_proxy = 0
# avg_deriv_trend$lpi_proxy_sd = 0
# for(i in 2:nrow(avg_deriv_trend)){
#     
#   avg_deriv_trend$lpi_proxy[i] = avg_deriv_trend$lpi_proxy[i-1] + avg_deriv_trend$avg_trend[i]
#   
#   # cumulative_sum_derivatives[i] = sum(avg_deriv_trend$avg_trend[1:i], na.rm = TRUE)
#   # cumulative_sum_derivatives_sd[i] = sum(temp2[1:i], na.rm = TRUE)
# }
# 
# avg_index = data.frame(
#   "time" = time+1981,
#   "index" = calc_index(avg_deriv_trend$avg_trend),
#   "lower" = calc_index(avg_deriv_trend$cilo),
#   "upper" = calc_index(avg_deriv_trend$cihi)
# )

# ggplot(data = avg_index) +
#   #geom_ribbon(aes(x = time, ymin = lower, ymax = upper)) +
#   geom_line(aes(x = time, y = index))

# 
# 
# 
# 
# ## calculate an index of change ----
# 
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
# index_allperms_perseries = lapply(derivs_respscale_ls, function(x) apply(x, 1, FUN = calc_index))
# # function(x) apply(x[runif(100, 1, nrow(x)),], 1, FUN = calc_index))
# index_perspecies = lapply(index_allperms_perseries,
#                           function(x) apply(x, 1, quantile, probs = c(.05, .5, .95)))
# names(index_perspecies) = colnames(YData)
# index_perspecies = lapply(index_perspecies, function(x) as.data.frame(t(x))) 
# index_perspecies = bind_rows(index_perspecies, .id = "species")
# colnames(index_perspecies) = c("species", "q05", "q50", "q95")
# index_perspecies$year = rep(time+1981, npops)
# 
# temp = index_perspecies |>
#   group_by(year) |>
#   summarise(mean_index = mean(q50))
# 
# # draw 100 random samples per species
# randoms = lapply(index_allperms_perseries, function(x) x[,runif(100, 1, nrow(x))])
# index_allperms = bind_cols(randoms)
# index_allperms_summary = data.frame(
#   "index" = temp$mean_index,
#   "lower" = apply(index_allperms, 1, quantile, probs = c(.05), na.rm = TRUE),
#   "upper" = apply(index_allperms, 1, quantile, probs = c(.95), na.rm = TRUE)
# )
# index_allperms_summary$year = time+1981
# 
# (model_index = ggplot(data = index_allperms_summary, aes(x = year)) +
#   geom_hline(yintercept = 1, lty = 2) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3, fill = "#6497b1") +
#   geom_line(aes(y = index)) +
#   labs(y = "Index of change in population size", x = ""))


# plot the derivatives per species with highlighted main 4 commercially fished species ----

# each species' median derivative 
derivs_respscale_pops = do.call(cbind, lapply(derivs_respscale_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE))) 
derivs_respscale_pops_lower = do.call(cbind, lapply(derivs_respscale_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_respscale_pops_upper = do.call(cbind, lapply(derivs_respscale_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))
colnames(derivs_respscale_pops) = colnames(Year_Geom_Means_all)
colnames(derivs_respscale_pops_lower) = colnames(Year_Geom_Means_all)
colnames(derivs_respscale_pops_upper) = colnames(Year_Geom_Means_all)
derivs_respscale_pops_df = derivs_respscale_pops |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year)
derivs_respscale_pops_df_lower = derivs_respscale_pops_lower |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "lower")
derivs_respscale_pops_df_upper = derivs_respscale_pops_upper |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "upper")
derivs_respscale_pops_df = full_join(derivs_respscale_pops_df, derivs_respscale_pops_df_lower) |>
  full_join(derivs_respscale_pops_df_upper)
# set first time step's derivative to 0
derivs_respscale_pops_df$value[which(derivs_respscale_pops_df$year == 1981)] <- 0
derivs_respscale_pops_df$lower[which(derivs_respscale_pops_df$year == 1981)] <- 0
derivs_respscale_pops_df$upper[which(derivs_respscale_pops_df$year == 1981)] <- 0

derivs_respscale_pops_df$name = gsub("_", " ", derivs_respscale_pops_df$name) |> 
  stringr::str_to_sentence()
derivs_respscale_pops_df$name =  factor(derivs_respscale_pops_df$name,
                              levels = temporal_trend$species[order(temporal_trend$mu_deriv)])

ggplot(data = derivs_respscale_pops_df, 
       aes(x = year, group = name)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = name), alpha = .3) +
  geom_line(aes(y = value), lwd = .2) +
  labs(x = "",
       fill = "Species",
       y = "Rate of change in biomass \n(kg per tow / year)") +
  gghighlight::gghighlight(name %in% c("Gadus morhua",
                                       "Hippoglossoides platessoides",
                                       "Reinhardtius hippoglossoides",
                                       "Sebastes mentella"),
                           use_direct_label = FALSE,
                           unhighlighted_params = list(colour = NULL)) +
    theme(legend.position = "right") 
ggsave("figures/species_derivative_trends.png", width = 8.33, height = 4.25)
