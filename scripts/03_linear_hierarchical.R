# Script to make a hierarchical linear model on all populations
# but with a yearly random effect

# libraries ----

library(here)
library(dplyr)
library(tidyr)
library(HMSC)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggridges)
theme_set(theme_pubr())

set.seed(12)

# data ----

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

biomass <- Year_Geom_Means_all %>% apply(2, scale, center = FALSE)
matplot(biomass, type = "l")
saveRDS(biomass, "outputs/biomass.rds")

################################################################################

# build the model with HMSC

# format data
XData <- as.matrix(time)
YData <- biomass
randomEff <- as.factor(1:tsl) 

# prepare data
formDat <- as.HMSCdata(X = XData, 
                       Y = YData, 
                       Random = randomEff, 
                       interceptX = TRUE)

# run model
m <- hmsc(formDat,
          family = "gaussian",
          niter = 15000,
          nburn = 5000,
          thin = 5,
          verbose = FALSE)
# save model
saveRDS(m, "outputs/linear_hierarchical.RDS")

# get Rsquared
HMSC::Rsquared(m) # this is more than double the multipop model!

# extract the parameters and format them into an `mcmc` object.
coefs = coef(m)

# check for parameter convergence
coda::traceplot(as.mcmc(m, parameters = "meansParamX"))
# coda::traceplot(as.mcmc(m, parameters = "paramX"))

################################################################################

# Predict and plot model over time 

# use a function that is modified to generate quantiles aorund the predictions
# this is in a PR on the HMSC package, not officially in the package
# https://github.com/guiblanchet/HMSC/pull/27

source("~/Documents/GitHub/covariationLPI/scripts/FUN_predict.hmsc.R", echo=TRUE)

# prepare and format data
formPredDat = as.HMSCdata(X = XData, Random = randomEff)

# predict population size over time from the model ----
pred <- predictQ.hmsc(m,
                newdata = formPredDat,
                type = "response")
for(i in 1:7){
  colnames(pred[[i]]) <- colnames(YData)
}
matplot(pred$mu, type = "l")

pred_l <- pred$mu |> as.data.frame() |>
  mutate("year" = time+1981) |>
  pivot_longer(cols = -year, values_to = "pred_biomass", names_to = "pop")

pred_q025_l <- pred$q025 |> as.data.frame() |>
  mutate("year" = time+1981) |>
  pivot_longer(cols = -year, values_to = "pred_biomass_q025", names_to = "pop")

pred_q975_l <- pred$q975 |> as.data.frame() |>
  mutate("year" = time+1981) |>
  pivot_longer(cols = -year, values_to = "pred_biomass_q975", names_to = "pop")

pred_l <- left_join(pred_l, pred_q025_l) |> left_join(pred_q975_l)
# save model
saveRDS(pred_l, "outputs/linear_hierarchical_pred_l.RDS")

# convert biomass to long, to plot in the background
biomass_l = as.data.frame(biomass) |> 
  mutate("year" = time+1981) |>
  tidyr::pivot_longer(cols = -year, names_to = "pop", values_to = "biomass")
saveRDS(biomass_l, "outputs/biomass_l.rds")

# plot predicted trend over the real trend for each population
ggplot() + 
  geom_ribbon(data = pred_l, aes(x = year, 
                                 ymin = pred_biomass_q025, 
                                 ymax = pred_biomass_q975, fill = pop), 
              alpha = .4) +
  geom_line(data = pred_l, aes(x = year, y = pred_biomass)) +
  geom_line(data = biomass_l, aes(x = year, y = biomass, col = pop), lty = 2) +
  theme(legend.position = "none") +
  facet_wrap(~pop)


# plot population slopes ----

# summarise parameters per population
paramX_mean = m$results$estimation$paramX |> apply(1:2, mean) # this is the same as coef$meansParamX
paramX_sd = m$results$estimation$paramX |> apply(1:2, sd)
colnames(paramX_sd) = paste0(colnames(paramX_sd), "_SD")

population_trends = as.data.frame(cbind(paramX_mean, paramX_sd))
population_trends$pop = rownames(population_trends)
population_trends$pop <- factor(population_trends$pop, 
                                levels = population_trends$pop[order(population_trends$x1)])
saveRDS(population_trends, "outputs/linear_hierarchical_population_trends.rds")

data_breaks <- data.frame(start = c(-1, -0.1, 0.1),  # Create data with breaks
                          end = c(-0.1, 0.1, 1),
                          colors = factor(c(1,3,2)))

(plot_poptrends = 
    ggplot(data = population_trends) +
    geom_rect(data = data_breaks,
              aes(xmin = start,
                  xmax = end,
                  ymin = - Inf,
                  ymax = Inf,
                  fill = colors), alpha = .2) +
    geom_segment(aes(x = x1 - x1_SD,
                     xend = x1 + x1_SD,
                     y = pop, yend = pop)) +
    geom_point(aes(y = pop, x = x1), size = 3) +
    labs(y = "", x = "alpha i") +
    geom_vline(xintercept = 0, lwd = .3) +
    coord_cartesian(xlim = c(-1, 1)))

# try it as a density plot
df = as.data.frame(m$results$estimation$paramX[,2,]) |>
  mutate("pop" = rownames(m$results$estimation$paramX)) |>
  pivot_longer(cols = -pop, names_to = "iteration", values_to = "slope")

df$pop <- factor(df$pop, levels = levels(population_trends$pop))

ggplot(data = df) +
  geom_rect(data = data_breaks,
            aes(xmin = start,
                xmax = end,
                ymin = - Inf,
                ymax = Inf,
                fill = colors), alpha = .7) +
  geom_density_ridges(aes(x = slope, y = pop), col = "grey10", fill = "white")
  
################################################################################
# plot the average temporal trend ----

# average slope with CI
paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.025, .5, .975))

# LPI-Community ----
coefs$meansParamX # community's variance (temporal trends) 
coefs$varX # community's variance (temporal trends) what is this???????

# overview of slopes in the community
m$results$estimation$paramX[,2,] |> hist(col = "grey70", border = "white")
abline(v = paramX_overall_mean[2], col = "blue", lwd = 2)
abline(v = paramX_overall_quantile[2,2], col = "red")
abline(v = paramX_overall_quantile[1,2], col = "red", lty = 2)
abline(v = paramX_overall_quantile[3,2], col = "red", lty = 2)

m$results$estimation$paramX[,2,] |> density() |> plot(frame = F, main = "Distribution of temporal slopes")
abline(v = paramX_overall_mean[2], col = "blue", lwd = 2)
abline(v = paramX_overall_quantile[2,2], col = "red")
abline(v = paramX_overall_quantile[1,2], col = "red", lty = 2)
abline(v = paramX_overall_quantile[3,2], col = "red", lty = 2)

# plotting the predicted trend across all populations
overall_pred = predictQ.hmsc(m, type = "response")
plot(overall_pred$overall_mu, type = "l", ylim = c(-.3, 3), lwd = 2)
lines(overall_pred$overall_q025)
lines(overall_pred$overall_q975)

# same plot, but with the population trends behind the trend line
# in ggplot
df_overall = data.frame(
  "year" = time+1981,
  "mu" = overall_pred$overall_mu,
  "q025" = overall_pred$overall_q025,
  "q975" = overall_pred$overall_q975
)
# save 
saveRDS(df_overall, "outputs/linear_hierarchical_df_overall.RDS")

ggplot() +
  geom_ribbon(data = pred_l,
              aes(x = year, group = pop, ymin = pred_biomass_q025, ymax = pred_biomass_q975), alpha = .3, fill = "salmon") +
  geom_line(data = pred_l,
            aes(x = year, y = pred_biomass, group = pop), linewidth = .2) +
  geom_ribbon(data = df_overall, aes(x = year, ymin = q025, ymax = q975), alpha = .6, fill = "dodgerblue4") +
  geom_line(data = df_overall, aes(x = year, y = mu), linewidth = 1, col = "white")

################################################################################

# plot species associations ----

# extract all estimated species-to-species associations matrix
assoMat <- corRandomEff(m)

# Average
siteMean <- apply(assoMat[, , , 1], 1:2, mean)
siteLower <- apply(assoMat[, , , 1], 1:2, quantile, prob = c(0.025))
siteUpper <- apply(assoMat[, , , 1], 1:2, quantile, prob = c(0.975))

# black out the weakest associations (that overlap with 0)
#siteMean[which(siteLower <= 0 & siteUpper >= 0)] <- NA

plot.new()
# Plot as heatmap
corrplot::corrplot(siteMean, type = "lower",
                   method = "color", 
                   tl.cex = .4, tl.col = "black", 
                   lowCI.mat = siteLower, 
                   uppCI.mat = siteUpper, na.label.col = "grey")


# plot as an ordination to see associations between years, or between species (or both)
source("scripts/biplot.hmsc.R")
plot.new()
biplot.hmsc(m, Random = 1, display = "sites", type = c("text", "points"), randomeff_names = as.character(time+1981))
biplot.hmsc(m, Random = 1, display = "species", type = "text", species_names = colnames(YData))
