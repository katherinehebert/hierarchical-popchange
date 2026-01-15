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

# save the scaling information for rescaling later
biomass_scaling = list(
  "biomass" = biomass,
  "YMeans" = YMeans,
  "baselines" = baselines,
  "weightings" = weightings,
  "time" = time,
  "npops" = npops,
  "tsl" = tsl
)
saveRDS(biomass_scaling, "outputs/biomass_scaling.rds")


# format into long
biomass$time = time
dat = pivot_longer(biomass, cols = -c(time), names_to = "series", values_to = "y")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time)
data_train = dat
saveRDS(data_train, "outputs/data_train.rds")

################################################################################
# check number of knots ----

knots = ceiling(tsl/4)

kcheck_gam = mgcv::gam(y ~ s(time, bs = "tp", k = 9) + s(series, bs = 're', k = npops),
                       data = data_train)
summary(kcheck_gam)
mgcv::k.check(kcheck_gam)
mgcv::gam.check(kcheck_gam)


################################################################################

# hierarchical gam on all populations ----

# prepare the priors ----

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
write.csv(test_priors, "outputs/test_priors.csv")

# look at the priors
plot(mvgam_prior, type = 'smooths')

png("figures/trend_priors.png", width = 1000, height = 1300, type = "cairo")
par(mfrow = c(6,5))
for(i in 1:ncol(Year_Geom_Means_all)){
  plot(mvgam_prior, type = 'trend', series = i)
} 
dev.off()

png("figures/re_priors.png", width = 1000, height = 700, type = "cairo")
plot(mvgam_prior, type = 're')
dev.off()

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
              samples = 10000,
              parellel = TRUE
)
#saveRDS(mod1, paste0("outputs/gam_hierarchical_gp.rds")) 
saveRDS(mod1, paste0("outputs/gam_hierarchical_gp_2025.rds")) 
