# Script to make a generalized additive model per population, then average them
# This is the current model we will compare to hierarchical models

# libraries ----

library(here)
library(tidyr)
library(dplyr)
library(bayesGAM)
library(tidybayes)
library(posterior)
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

# gam per population ----

for(i in 1:npops){
  df = data.frame(
    "biomass" = biomass[,i],
    "time" = time
  )
  m = bayesGAM(biomass ~ time, data = df)
  saveRDS(m, paste0("outputs/gam_singlepop/model_", i, ".rds"))
}


## plot per population ----

posterior_ls = list()
for(i in 1:npops){
  # Read in the model 
  m = readRDS(paste0("outputs/gam_singlepop/model_", i, ".rds")) 
  
  # extract posterior draws in an array format
  draws_fit = m |> getStanResults() |> as_draws_matrix()
  posterior_ls[[i]] = posterior::summarize_draws(draws_fit)
}
names(posterior_ls) = colnames(biomass)
posterior <- bind_rows(posterior_ls, .id = "sp")
saveRDS(posterior, "outputs/gam_singlepop/posterior.rds")

summary_ls = list()
for(i in 1:npops){
  # Read in the model 
  m = readRDS(paste0("outputs/gam_singlepop/model_", i, ".rds")) 
  
  summary_ls[[i]] = as.data.frame(summary(m))
  summary_ls[[i]]$variable = rownames(summary_ls[[i]])
}
names(summary_ls) = colnames(biomass)
summary_df <- bind_rows(summary_ls, .id = "pop")
saveRDS(summary_df, "outputs/gam_singlepop/summary.rds")

poptrends = summary_df |>
  dplyr::select(c(pop, variable, mean, sd)) |>
  dplyr::filter(variable %in% c("(Intercept)", "time")) |>
  tidyr::pivot_wider(names_from = variable, values_from = c(mean, sd)) |>
  rename("Intercept" = "mean_(Intercept)",
         "x1" = "mean_time",
         "Intercept_SD" = "sd_(Intercept)",
         "x1_SD" = "sd_time") |>
  relocate(pop, .after = "x1_SD")
poptrends$pop = as.factor(poptrends$pop)
saveRDS(poptrends, "outputs/gam_singlepop_population_trends.rds")


################################################################################

# assemble all draws in one matrix
for(i in 1:npops){
  # Read in the model 
  m = readRDS(paste0("outputs/gam_singlepop/model_", i, ".rds"))

  if(i > 1){
    # extract posterior draws in an array format
    draws_fit = m |> getStanResults() |> as_draws_matrix()
    all_draws <- rbind(all_draws, draws_fit)
  } else {all_draws <- m |> getStanResults() |> as_draws_matrix()}
}
saveRDS(all_draws, "outputs/gam_singlepop_alldraws.rds")
all_summary = summarise_draws(all_draws)

# predict each population's trend from the model
pop_trends = list()
for(i in 1:npops){
  
  # Read in the model 
  m = readRDS(paste0("outputs/gam_singlepop/model_", i, ".rds"))
  
  pp = bayesGAM::posterior_predict(m, draws = 100)

  pop_trends[[i]] = data.frame(
    "year" = time+1981,
    "biomass" = pp@pp$yrep |> apply(2, mean),
    "cilo" = pp@pp$yrep |> apply(2, quantile, prob = .05),
    "cihi" = pp@pp$yrep |> apply(2, quantile, prob = .95)
  )
}
names(pop_trends) = colnames(biomass)
pop_trends = bind_rows(pop_trends, .id = "sp")
saveRDS(pop_trends, "outputs/gam_singlepop_pred_l.rds")

## plot the average trend  ----

avg_trend = pop_trends |>
  group_by(year) |>
  summarise(biomass = mean(biomass),
            cilo = mean(cilo),
            cihi = mean(cihi))
saveRDS(avg_trend, "outputs/gam_singlepop/avg_trend.rds")

# convert biomass to long
biomass_l = as.data.frame(biomass) |> 
  mutate("year" = time+1981) |>
  tidyr::pivot_longer(cols = -year, names_to = "sp", values_to = "biomass")

ggplot() +
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = sp), size = .3, alpha = .4, col = "grey") +
  geom_ribbon(data = avg_trend,
              aes(x = year, ymin = cilo, ymax = cihi), alpha = .3, fill = "red") + # credible interval?
  geom_ribbon(data = pop_trends,
             aes(x = year, ymin = cilo, ymax = cihi, group = sp), alpha = .01, fill = "blue") +
  geom_line(data = pop_trends,
            aes(x = year, y = biomass, group = sp), linewidth = .1, col = "blue") +
  geom_line(data = avg_trend,
          aes(x = year, y = biomass), linewidth = 1, col = "red") 
