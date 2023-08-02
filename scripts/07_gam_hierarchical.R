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

YData <- Year_Geom_Means_all # |> apply(2, scale, center = FALSE)
# instead of centering on the mean, center on the baseline biomass
for(i in 1:ncol(YData)){
  YData[,i] = YData[,i] - YData[1,i]
}

biomass <- YData |>
  #apply(2, scale, center = FALSE) |> # try without scaling? the non-centered modelled looks crazy. 
  as.data.frame()
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
                                  #s(time, bs = "tp", k = knots) + 
                                  # independent smoothers for each group
                                  s(time, by = series, bs = "tp", k = knots) +
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
              formula = y ~ s(time, bs = "tp", k = knots) + s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'GP',
              use_stan = TRUE,
              chains = 3, 
              burnin = 5000,
              samples = 15000
)
saveRDS(mod1, paste0("outputs/gam_hierarchical.rds")) # 4% ended in divergence

# to view the Stan model file:
m <- mod1
# save the model file
write_stan_file(m$model_file, "stan")
summary(mod1)

# convergence?
rstan::stan_trace(mod1$model_output, 'b')

## plot the temporal trend (time smooth) ----

plot_mvgam_smooth(mod1, smooth = "time")
mvgam::plot_mvgam_randomeffects(mod1)
plot_mvgam_factors(mod1)

# extract posterior draws in an array format
draws_fit = m$model_output |> posterior::as_draws_matrix()
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
  mutate(time = as.integer(time)+1980)
saveRDS(predictions, "outputs/gam_hierarchical_pred_l.rds")

# plot the predictions next to the data ----
data_train = rename(data_train, "sp" = "series")
data_train$sp = gsub("_", "\n", data_train$sp) |> stringr::str_to_sentence()
predictions$sp = gsub("_", "\n", predictions$sp) |> stringr::str_to_sentence()
ggplot(data = predictions,
       aes(x = time, group = sp)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = sp), alpha = .2) +
  geom_line(aes(y = biomass)) +
  geom_line(data = data_train,
            aes(x = time+1981, y = y), lty = 2) +
  facet_wrap(~sp, scales = "free_y") +
  theme(legend.position = "none",strip.text = element_text(face = "italic")) +
  labs(x = "", y = "Biomass (observed)")
ggsave("figures/mvgam_predictions_perspecies.png", width = 12.8, height = 7.63)


# mean trend from the model
# posterior_mus = posterior_df[grep("mus", posterior_df$variable),] |> 
#   separate(col = variable, into = c("year", "sp"), sep = ",") |>
#   mutate(year = readr::parse_number(year),
#          sp = readr::parse_number(sp))

# take the derivative of each population trend
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

# each species' mean predicted trend in biomass
preds_pops = do.call(cbind, lapply(preds_ls, FUN = function(x) apply(x, 2, mean, na.rm = TRUE))) 
matplot(preds_pops, x = time+1981, type = "l")
# each species' mean derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, mean, na.rm = TRUE)))
derivs_pops_sd = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, sd, na.rm = TRUE)))
matplot(derivs_pops, x = time+1981, type = "l")

# # plot the temporal basis coefficients ----
# # note: (b2 through b9, when 9 knots are used) but need to update if more/less knots are used
# temporary = draws_fit[,grep("b\\[", colnames(draws_fit))]
# b_time = temporary[,2:9]
# plot(density(b_time))
# abline(v = mean(b_time, na.rm = TRUE))
# b_time_long = b_time |> 
#   as.data.frame() |> 
#   pivot_longer(cols = everything())
# 
# (plot_trenddensity = ggdensity(b_time_long, x = "value",
#                                add = "mean", rug = FALSE, add.params = list(linetype = 1)) +
#     labs(x = "Posterior estimates of basis coefficients (time)", y = "Density") +
#     geom_vline(xintercept = 0, lwd = .3, lty = 2) +
#     theme(legend.position = "none", panel.grid.major.x = element_line()))

# plot each species' mean derivative over the whole time series ----
(temporal_trend = data.frame(
  species = colnames(Year_Geom_Means_all),
  mu_deriv = derivs_pops |> apply(2, mean, na.rm = TRUE),
  sd_deriv = derivs_pops |> apply(2, sd, na.rm = TRUE)
  #lower = derivs_pops |> apply(2, quantile, probs = .025, na.rm = TRUE),
  #upper = derivs_pops |> apply(2, quantile, probs = .975, na.rm = TRUE)
))
temporal_trend$lower = temporal_trend$mu_deriv - temporal_trend$sd_deriv
temporal_trend$upper = temporal_trend$mu_deriv + temporal_trend$sd_deriv
temporal_trend$species = gsub("_", " ", temporal_trend$species) |> stringr::str_to_sentence()
temporal_trend$species =  factor(temporal_trend$species,
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
  labs(x = "Annual rates of change\n of all species", 
       y = "Frequency", 
       fill = "") +
  scale_fill_distiller(palette = "RdBu", 
                       direction = 1, 
                       limits = c(-2.5,2.5)))

# plot the average derivative with sd limits per species ----
(plot_poptrends = ggplot(data = temporal_trend[order(temporal_trend$mu_deriv),]) +
  geom_segment(aes(x = lower, xend = upper,
                   y = species, yend = species)) +
  geom_point(aes(x = mu_deriv, y = species, fill = mu_deriv), size = 4, pch = 21) +
  labs(x = "Average rate of change (1981-2013)", 
       y = "",
       fill = "") +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1,1)) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_pubr() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major.x = element_line()) 
)

# arrange the plot panels and save
(plot_trenddensity + 
    coord_cartesian(xlim = c(-2.5, 2.5)) +
    theme(axis.title.y = element_text(vjust = -60)) #+
    # annotate("text", x = coefs$meansParamX[2]-0.4, y = 1.8, 
    #          label = paste0("\u03bc = ", round(coefs$meansParamX[2], digits = 3),
    #                         "\n \u03c3\u00b2 = ", round(coefs$varX[2,2], digits = 3)))
  ) / (plot_poptrends + 
         coord_cartesian(xlim = c(-2.5, 2.5)) +
     theme(axis.text.y = element_text(face = "italic"),
           panel.grid.major.x = element_line()) + 
     scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1,1))) +
  plot_layout(heights = (c(1,3))) +
  plot_annotation(tag_levels = "a")
ggsave("figures/distribution_poptrends.png", width = 8.23, height = 9)


colnames(derivs_pops) = colnames(Year_Geom_Means_all)
colnames(derivs_pops_sd) = colnames(Year_Geom_Means_all)
derivs_pops_df = derivs_pops |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year)
temp_sd = derivs_pops_sd |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "sd")
derivs_pops_df = full_join(derivs_pops_df, temp_sd)
# set first time step's derivative to 0
derivs_pops_df$value[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$sd[which(derivs_pops_df$year == 1981)] <- 0

derivs_pops_df$name = gsub("_", " ", derivs_pops_df$name) |> 
  stringr::str_to_sentence()
derivs_pops_df$name =  factor(derivs_pops_df$name,
                                 levels = temporal_trend$species[order(temporal_trend$mu_deriv)])

# ridgeplot of derivatives per species
# (plot_deriv_densities = ggplot(data = derivs_pops_df) +
#   geom_density_ridges_gradient(aes(x = value, y = name, fill = after_stat(x)),
#                                size = .2, scale = 5) +
#   scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-3,3)) +
#   theme(panel.grid.major.x = element_line(),
#       ) + 
#   labs(x = "Annual rate of change (1981-2013)",
#        y = "",
#        fill = ""))

derivs_pops_df$cumulative_deriv = NA
derivs_pops_df$cumulative_deriv_sd = NA
for(n in unique(derivs_pops_df$name)){
  temp = derivs_pops_df$value[which(derivs_pops_df$name == n)]
  temp2 = derivs_pops_df$sd[which(derivs_pops_df$name == n)]
  
  cumulative_sum_derivatives = c(0)
  cumulative_sum_derivatives_sd = c(0)
  for(i in 2:length(temp)){
    cumulative_sum_derivatives[i] = sum(temp[1:i], na.rm = TRUE)
    cumulative_sum_derivatives_sd[i] = sum(temp2[1:i], na.rm = TRUE)
  }
  
  derivs_pops_df$cumulative_deriv[which(derivs_pops_df$name == n)] <- cumulative_sum_derivatives
  derivs_pops_df$cumulative_deriv_sd[which(derivs_pops_df$name == n)] <- cumulative_sum_derivatives_sd
}

# plot estimated derivative trend ----
# ggplot(data = derivs_pops_df) +
#   geom_line(aes(y = value, x = year, group = name), lwd = .2) +
#   coord_cartesian(ylim = c(-3,1)) +
#   labs(x = "",
#        y = "Estimated rate of change in biomass") +
#   theme(legend.position = "right") +
#   gghighlight::gghighlight(name %in% c("Gadus morhua", 
#                                        "Hippoglossoides platessoides",
#                                        "Reinhardtius hippoglossoides",
#                                        "Sebastes mentella"),
#                            use_direct_label = FALSE
#   )
(plot_deriv = ggplot(data = derivs_pops_df) +
    geom_ribbon(aes(ymin = value - sd,
                    ymax = value + sd,
                    x = year, group = name, fill = name), alpha = .3) +
    geom_line(aes(y = value, x = year, group = name), lwd = .2) +
    labs(x = "", 
         y = "Annual rate of change in biomass") +
    theme(legend.position = "none") +
    gghighlight::gghighlight(name %in% c("Gadus morhua", 
                                         "Hippoglossoides platessoides",
                                         "Reinhardtius hippoglossoides",
                                         "Sebastes mentella"),
                             use_direct_label = FALSE,
                             unhighlighted_params = list(colour = NULL))
)
# plot cumulative derivative trend ----
(plot_cumu_deriv = ggplot(data = derivs_pops_df) +
    geom_ribbon(aes(ymin = cumulative_deriv - cumulative_deriv_sd,
                    ymax = cumulative_deriv + cumulative_deriv_sd,
                    x = year, group = name, fill = name), alpha = .3) +
  geom_line(aes(y = cumulative_deriv, x = year, group = name), lwd = .2) +
    labs(x = "", y = "Cumulative change in biomass", fill = "Species") +
    theme(legend.position = "right",
          legend.text = element_text(face = "italic")) +
  gghighlight::gghighlight(name %in% c("Gadus morhua", 
                                       "Hippoglossoides platessoides",
                                       "Reinhardtius hippoglossoides",
                                       "Sebastes mentella"),
                           use_direct_label = FALSE,
                           unhighlighted_params = list(colour = NULL))
)
# Anarhichas denticulatus is the species that has declined a lot too, but isn't 
# one of the big commercially fished ones. it is considered threatened by COSEWIC since 2001
#|> plotly::ggplotly() 

plot_deriv + plot_cumu_deriv
ggsave("figures/mvgam_rateofchange.png", width = 11, height = 4.06)

# should this be avg by species, then avg over community? ----
avg_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(preds, 2, mean),
  "cilo" = apply(preds, 2, quantile, prob = .05),
  "cihi" = apply(preds, 2, quantile, prob = .95),
  "sd" = apply(preds, 2, sd)
)
saveRDS(avg_trend, "outputs/gam_hierarchical_df_overall.rds")

avg_deriv_trend = data.frame(
  "year" = time+1981,
  "avg_trend" = apply(derivs, 2, mean, na.rm = T),
  "cilo" = apply(derivs, 2, quantile, prob = .05, na.rm = T),
  "cihi" = apply(derivs, 2, quantile, prob = .95, na.rm = T),
  "sd" = apply(derivs, 2, sd)
)
avg_deriv_trend[1,2:5] = 0
# avg_deriv_trend = data.frame(
#   "year" = time+1981,
#   "avg_trend" = apply(derivs_pops, 1, mean, na.rm = TRUE),
#   "cilo" = apply(derivs_pops, 1, quantile, prob = .025, na.rm = TRUE),
#   "cihi" = apply(derivs_pops, 1, quantile, prob = .975, na.rm = TRUE),
#   "sd" = apply(derivs_pops, 1, sd)
# )
saveRDS(avg_deriv_trend, "outputs/gam_hierarchical_df_overall_deriv.rds")

# plot average biomass trend
(A = ggplot(data = avg_trend, aes(x = year)) +
    # geom_line(data = predictions, aes(x = time, y = biomass, group = sp),
    #           linewidth = .2, col = "grey") +
    geom_hline(yintercept = 0, lwd = .3) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = 1) +
    labs(x = "",
         y = "Biomass") +
    coord_cartesian(ylim = c(-12,12)) +
    theme(panel.grid.major.x = element_line()))
# plot avg derivative trend
(B = ggplot(data = avg_deriv_trend, aes(x = year)) +
    geom_ribbon(aes(ymin = cilo, ymax = cihi), 
                alpha = .3, fill = "#6497b1") +
    geom_line(aes(y = avg_trend), col = "#03396c", lwd = 1) +
    geom_hline(yintercept = 0, lwd = .3) +
    labs(x = "",
         y = "Rate of change") +
    coord_cartesian(ylim = c(-1.3, 1.3)) +
    theme(panel.grid.major.x = element_line()))
A + B + plot_annotation(tag_levels = "a")
ggsave("figures/mvgam_prediction.png", width = 8.46, height = 4.06)


# get species correlations
sp_correlations = lv_correlations(mod1)
saveRDS(sp_correlations, "outputs/gam_hierarchical_species_correlations.rds")

# sp_correlations = readRDS("outputs/gam_hierarchical_species_correlations.rds")

# Plot as heatmap
png(height=1800, width=1800, file="figures/species_associations.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

