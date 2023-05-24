# Script to make a linear model per population, then average them
# This is the current model we will compare to hierarchical models

# libraries ----

library(here)
library(dplyr)
library(brms)
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

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

biomass <- Year_Geom_Means_all %>% apply(2, scale, center = TRUE)
matplot(biomass, type = "l")

################################################################################

# linear model per population ----

for(i in 1:npops){
  df = data.frame(
    "biomass" = biomass[,i],
    "time" = time
  )
  m = brm(biomass ~ time + 1, data = df)
  saveRDS(m, paste0("outputs/linear_singlepop/model_", i, ".rds"))
}


## plot per population ----

posterior_ls = list()
par(mfrow = c(3,2))
for(i in 1:npops){
  # Read in the model 
  m = readRDS(paste0("outputs/linear_singlepop/model_", i, ".rds"))
  m_data = m$data
  
  # extract posterior draws in an array format
  (draws_fit <- as_draws_matrix(m))
  post = posterior::summarize_draws(draws_fit)
  
  # plot the fitted Line
  p = ggplot() +
    geom_point(data = m_data, aes(x = time, y = biomass)) +
    stat_function(fun = function(x) post$mean[1] + post$mean[2] * x)
  print(p)
  
  posterior_ls[[i]] = post
}
names(posterior_ls) = colnames(biomass)
posterior <- bind_rows(posterior_ls, .id = "sp")
saveRDS(posterior, "outputs/linear_singlepop/posterior.rds")

poptrends = posterior |>
  select(c(sp, variable, mean, sd)) |>
  dplyr::filter(variable %in% c("b_Intercept", "b_time")) |>
  tidyr::pivot_wider(names_from = variable, values_from = c(mean, sd)) |>
  rename("Intercept" = "mean_b_Intercept",
         "x1" = "mean_b_time",
         "Intercept_SD" = "sd_b_Intercept",
         "x1_SD" = "sd_b_time",
         "pop" = "sp") |>
  relocate(pop, .after = "x1_SD")
poptrends$pop = as.factor(poptrends$pop)
saveRDS(poptrends, "outputs/linear_singlepop_population_trends.rds")


################################################################################

# calculate the average trend across populations ----
posterior_slopes = dplyr::filter(posterior, variable == "b_time")
posterior_int = dplyr::filter(posterior, variable == "b_Intercept")

# assemble all draws in one matrix
for(i in 1:npops){
  # Read in the model 
  m = readRDS(paste0("outputs/linear_singlepop/model_", i, ".rds"))
  m_data = m$data
  
  if(i > 1){
  # extract posterior draws in an array format
  all_draws <- rbind(all_draws, as_draws_matrix(m))
  } else {all_draws <- as_draws_matrix(m)}
}
saveRDS(all_draws, "outputs/linear_singlepop/all_draws.rds")
all_summary = summarise_draws(all_draws)

# note: this is a 90% credible interval
avg_trend = data.frame(
  "year" = time+1981,
  "biomass" = all_summary$mean[2]*time + all_summary$mean[1],
  "cilo" = all_summary$q5[2]*time + all_summary$q5[1],
  "cihi" = all_summary$q95[2]*time + all_summary$q95[1]
)
saveRDS(avg_trend, "outputs/linear_singlepop/avg_trend.rds")

# predict
pop_trends = list()
for(i in 1:npops){
  pop_trends[[i]] = data.frame(
    "year" = time+1981,
    "biomass" = posterior_slopes$mean[i]*time + posterior_int$mean[i],
    "cilo" = posterior_slopes$q5[i]*time + posterior_int$q5[i],
    "cihi" = posterior_slopes$q95[i]*time + posterior_int$q95[i]
  )
}
names(pop_trends) = colnames(biomass)
pop_trends = bind_rows(pop_trends, .id = "sp")
saveRDS(pop_trends, "outputs/linear_singlepop_pred_l.rds")

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
