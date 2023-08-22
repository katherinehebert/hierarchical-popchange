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

## plot the temporal smooth ----

source("~/Documents/GitHub/hierarchical-lpi/plot_mvgam_smooth_custom.R", echo=TRUE)
temporary = plot_mvgam_smooth_custom(mod1, smooth = "time") # this is essentially the LPI curve
temporary
