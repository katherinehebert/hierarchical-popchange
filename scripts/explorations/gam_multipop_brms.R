# Script to run a multivariate GAM on population trends using the package brms

library(here)
library(dplyr)
library(tidyr)
library(brms)
library(tidybayes)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())

set.seed(12)

# data ----

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) |> as.numeric()
time <- scale(time, center = TRUE, scale = FALSE) # must run m1-5 with this change!!!
time_m <- as.matrix(time)
year_effect <- as.character(time)

# maybe center time?

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

biomass <- Year_Geom_Means_all |> 
  apply(2, scale, center = TRUE) |>
  as.data.frame()
matplot(biomass, type = "l")

# format into long
biomass$time = time
dat = pivot_longer(biomass, cols = -c(time), names_to = "series", values_to = "y")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time)
dat$year_effect = as.factor(dat$time)
data_train = dat

### check number of knots to use for the time spline ----
## using mgcv for a quick look 

find_k <- mgcv::gam(y ~ s(time, k = 10), data = data_train)
mgcv::k.check(find_k)
plot(find_k) # the default 10 seems to be ok. not great, but adding more doesn't improve this much

### write and run the model -----

# this ignores the species names for now. will make multivariate later

m1 <- brm(bf(y ~ s(time)),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, thin = 10, # made this lighter just for now. was 4000 and 1000 and thin 10
          refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m1, "models/brms_m1.rds")
summary(m1)
# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# exquisite breakdown of what the summary means

msms = conditional_smooths(m1)
plot(msms)

## posterior predictive checks -----

# density plot overlay of the original response values (the thick black line) 
# with 10 draws from the posterior distribution of the model
pp_check(m1)

# empirical cumulative distribution function of the observations and random draws 
# from the model posterior
pp_check(m1, type = "ecdf_overlay")


### adding species effect ----

#### write and run the model -----

m2 <- brm(bf(y ~ s(time, bs = "tp") + s(series, bs = "re")),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, # made this lighter just for now. was 4000 and 1000
          thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m2, "models/brms_m2.rds")

summary(m2) # the series smooth doesn't add much info to the model
# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# exquisite breakdown of what the summary means

msms = conditional_smooths(m2)
plot(msms)

## posterior predictive checks -----

# density plot overlay of the original response values (the thick black line) 
# with 10 draws from the posterior distribution of the model
pp_check(m2)

# empirical cumulative distribution function of the observations and random draws 
# from the model posterior
pp_check(m2, type = "ecdf_overlay")


### adding time random effect ----

m3 <- brm(bf(y ~ s(time, bs = "tp") + 
               s(year_effect, bs = "re")),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, # made this lighter just for now. was 4000 and 1000
          thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m3, "models/brms_m3.rds")
summary(m3) # year effect isn't adding much, it seems

# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# exquisite breakdown of what the summary means

msms = conditional_smooths(m3)
plot(msms)


## posterior predictive checks -----

# density plot overlay of the original response values (the thick black line) 
# with 10 draws from the posterior distribution of the model
pp_check(m3)

# empirical cumulative distribution function of the observations and random draws 
# from the model posterior
pp_check(m3, type = "ecdf_overlay")


## adding year random effect and a time smoother per species level -----
# this one is extremely slow...

m4 <- brm(bf(y ~ s(time, bs = "tp") + 
               s(time, by = series) +
               s(year_effect, bs = "re")),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, # made this lighter just for now. was 4000 and 1000
          thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m4, "models/brms_m4.rds") # not converged. needs more time or better priors
summary(m4) 

# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# exquisite breakdown of what the summary means

msms = conditional_smooths(m4)
plot(msms)


## posterior predictive checks -----

# density plot overlay of the original response values (the thick black line) 
# with 10 draws from the posterior distribution of the model
pp_check(m4)

# empirical cumulative distribution function of the observations and random draws 
# from the model posterior
pp_check(m4, type = "ecdf_overlay")



## year and species random effects -----

m5 <- brm(bf(y ~ s(time, bs = "tp") + 
               s(series, bs = "re") +
               s(year_effect, bs = "re")),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, # made this lighter just for now. was 4000 and 1000
          thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m5, "models/brms_m5.rds") 
summary(m5) 

# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# exquisite breakdown of what the summary means

msms = conditional_smooths(m5)
plot(msms)


## posterior predictive checks -----

# density plot overlay of the original response values (the thick black line) 
# with 10 draws from the posterior distribution of the model
pp_check(m5)

# empirical cumulative distribution function of the observations and random draws 
# from the model posterior
pp_check(m5, type = "ecdf_overlay")



## adding a time smoother per species level w/o year effect -----
# this one is extremely slow...

m6 <- brm(bf(y ~ s(time, by = series)),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, # made this lighter just for now. was 4000 and 1000
          thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m6, "models/brms_m6.rds") # not converged. needs more time or better priors
summary(m6) 

# https://fromthebottomoftheheap.net/2018/04/21/fitting-gams-with-brms/
# exquisite breakdown of what the summary means

msms = conditional_smooths(m6)
plot(msms)


## posterior predictive checks -----

# density plot overlay of the original response values (the thick black line) 
# with 10 draws from the posterior distribution of the model
pp_check(m6)

# empirical cumulative distribution function of the observations and random draws 
# from the model posterior
pp_check(m6, type = "ecdf_overlay")


# get var-covar matrix
vcm = vcov(m6, correlation = TRUE)
# Plot as heatmap
corrplot::corrplot(vcm, 
                   type = "lower",
                   method = "color", 
                   tl.cex = .4, tl.col = "black") # seems weak...



## time with a smoother for time per species
# to run
m7 <- brm(bf(y ~ s(time) + s(time, by = series)),
          data = data_train, 
          family = gaussian(), 
          cores = 4, seed = 17,
          iter = 3000, warmup = 1000, # made this lighter just for now. was 4000 and 1000
          thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))
saveRDS(m7, "models/brms_m7.rds") # not converged. needs more time or better priors
summary(m7) 


