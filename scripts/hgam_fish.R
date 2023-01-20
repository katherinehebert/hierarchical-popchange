# Script to build hierarchical GAM on fish population trends

library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(ggplot2)

theme_set(ggpubr::theme_pubr())

#### import & set up data ####

# from https://github.com/eric-pedersen/groundfish-data-analysis
# note: these are mean biomasses over the whole sampled area with an SE
# not actually raw abundances. But for this purpose, should be ok
load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")

time = rownames(Year_Geom_Means_all) %>% as.numeric()

# get number of populations & length of time series
npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

# scale all columns
Year_Geom_Means_all = apply(Year_Geom_Means_all, 2, scale) %>% as.data.frame()
Year_Geom_Means_all$time = as.integer(time) 
df = Year_Geom_Means_all %>% 
  pivot_longer(cols = -time, names_to = "species", "values_to" = "N")
df$species <- factor(df$species)

# build a GI model: global smoother plus group-level smoothers with differing wiggliness

modGI <- gam(N ~ 
               # global smooth
               s(time, bs = c("tp"), m = 2, k = 30) + 
               # smooth per species
               s(time, by = species, bs = c("tp"), m = 1, k = 30) +
               # random intercept
               s(species, bs = c("re")), 
             data = df, 
             method="REML", 
             family="gaussian")

gam.check(modGI) # not great
k.check(modGI)

# quartz(width = 12, height = 12)
# draw(modGI)

summary(modGI)

# drawing posteriors for the predicted response values:

beta <- coef(modGI)   # vector of model parameters

Vb <- vcov(modGI)     # default is the bayesian covariance matrix
heatmap(Vb, Rowv = NA, Colv = NA)

# Draw the estimated smooth 
# plot(modGI, seWithMean = TRUE, unconditional = TRUE)

# Now let's explore working with smooths in a little more detail. First evaluate the smooth over a range of values for `year`
time_sm <- smooth_estimates(modGI, "s(time)")
# add a confidence interval
time_sm <- confint(modGI, "s(time)")
time_sm

ggplot(data = filter(time_sm, smooth == "s(time)"), 
       aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .4) +
  geom_line(aes(y = est)) +
  labs(y = "N", x = "") +
  ggpubr::theme_pubr()

ggplot(data = time_sm, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = smooth), alpha = .1) +
  geom_line(aes(y = est, col = smooth)) +
  geom_ribbon(data = filter(time_sm, smooth == "s(time)"),
              aes(x = time, ymin = lower, ymax = upper), alpha = .4) +
  geom_line(data = filter(time_sm, smooth == "s(time)"),
            aes(x = time, y = est)) +
  labs(y = "N", x = "") +
  ggpubr::theme_pubr() +
  theme(legend.position = "none") 


# sample from the posterior distribution
modGI_post <- smooth_samples(modGI, 's(time)', n = 100, seed = 42) 

post_summary <- modGI_post %>% 
  group_by(.x1) %>% 
  summarise(mu = mean(value, na.rm = TRUE),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975),
            variance = var(value)) %>% ungroup()

time_fit <- evaluate_smooth(modGI, 's(time)')

draw(time_fit) + 
  geom_line(data = modGI_post, 
            aes(x = .x1, y = value, group = draw),
            alpha = 0.3, colour = 'red') 
# why does the global smooth start out at a higher intercept?
