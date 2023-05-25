# Script to make a multivariate generalized additive model 
# but without a year effect

# libraries ----

library(here)
library(dplyr)
library(bayesGAM)
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

biomass <- Year_Geom_Means_all %>% apply(2, scale, center = TRUE)
matplot(biomass, type = "l")

################################################################################

# gam on all populations ----

m = bayesGAM(biomass ~ time_m[,1] + 1)
saveRDS(m, paste0("outputs/gam_multipop.rds"))

# extract posterior draws in an array format
draws_fit = m |> getStanResults() |> as_draws_matrix()
posterior_df = posterior::summarize_draws(draws_fit)
saveRDS(posterior_df, "outputs/gam_multipop_posterior.rds")

# extract the model summary
summary_df = as.data.frame(summary(m))
summary_df$variable = rownames(summary_df)
# assign the population name to the numbers in the parameter names
pops = data.frame(
  "pop" = colnames(biomass),
  "number" = 1:ncol(biomass)
)
pop_numbers = data.frame(
  "number" = readr::parse_number(summary_df$variable)
)
pop_numbers = left_join(pop_numbers, pops)
summary_df$pop = pop_numbers$pop
saveRDS(summary_df, "outputs/gam_multipop_summary.rds")

# format into the population trends
poptrends = summary_df |>
  tidyr::separate(variable, sep = "\\[", into = c("variable", "pop_number")) |>
  select(c(pop, variable, mean, sd)) |>
  dplyr::filter(variable %in% c("beta_(Intercept)", "beta_time_m")) |>
  tidyr::pivot_wider(names_from = variable, values_from = c(mean, sd)) |>
  rename("Intercept" = "mean_beta_(Intercept)",
         "x1" = "mean_beta_time_m",
         "Intercept_SD" = "sd_beta_(Intercept)",
         "x1_SD" = "sd_beta_time_m") |>
  relocate(pop, .after = "x1_SD")
poptrends$pop = as.factor(poptrends$pop)
saveRDS(poptrends, "outputs/gam_multipop_population_trends.rds")


# predict from the model
pp = posterior_predict(m, draws = 100)

pop_trends_ls = list()
for(i in 1:npops){
  pop_trends_ls[[i]] = data.frame(
    "year" = time+1981,
    "biomass" = pp@pp[[i]] |> apply(2, mean),
    "cilo" = pp@pp[[i]] |> apply(2, quantile, prob = .05),
    "cihi" = pp@pp[[i]] |> apply(2, quantile, prob = .95)
  )
}
names(pop_trends_ls) = colnames(biomass)
pop_trends = bind_rows(pop_trends_ls, .id = "sp")
saveRDS(pop_trends, "outputs/gam_multipop_pred_l.rds")

# this is a mean of a mean. should fix this ----
avg_trend = pop_trends |>
  group_by(year) |>
  summarise(biomass = mean(biomass),
            cilo = mean(cilo),
            cihi = mean(cihi))
saveRDS(avg_trend, "outputs/gam_multipop_avg_trend.rds")

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
