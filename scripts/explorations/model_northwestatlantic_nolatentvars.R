# Script to model Northwest Atlantic population trends with HMSC without latent variables

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

# custom functions
source("scripts/predictQ.hmsc.R") # this is a customized version of HMSC's predict
# this one produces quantiles of the predictions in addition to the mean prediction

# data ----

# from https://github.com/eric-pedersen/groundfish-data-analysis
# note: these are mean biomasses over the whole sampled area with an SE
# not raw abundances. but they are useful as population biomasses
load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")

time = rownames(Year_Geom_Means_all) %>% as.numeric()
time_uncentered = time

# get number of populations & length of time series
npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

################################################################################

## build the model with HMSC ----

# format data

XData <- as.matrix(time)

YData <- as.matrix(Year_Geom_Means_all) 

YData <- apply(YData, 2, scale, center = FALSE)

# instead of centering on the mean, center on the baseline biomass
for(i in 1:ncol(YData)){
  YData[,i] = YData[,i] - YData[1,i]
}

# prepare data
formDat <- as.HMSCdata(X = XData, 
                       Y = YData, 
                       interceptX = TRUE,
                       scaleX = TRUE)

# set flat priors
priors <- as.HMSCprior(data = formDat, family = "gaussian")

# run model
m <- hmsc(formDat,
          family = "gaussian",
          niter = 15000,
          nburn = 1000,
          thin = 1,
          verbose = FALSE,
          priors = priors)
# save model
saveRDS(m, "models/hmsc_nolatentvars.RDS")


################################################################################

## Model convergence & explanatory power ----

# extract the parameters and format them into an `mcmc` object.
mcmcMeansParamX <- as.mcmc(m, parameters = "meansParamX")

# check for parameter convergence 
bayesplot::mcmc_trace(mcmcMeansParamX)
ggsave("figures/traceplots_meansparamX_nolatentvars.png", width = 10, height = 4)

# Explanatory power ----
(R2comm <- Rsquared(m, averageSp = TRUE, type = "ols"))
(R2pops <- Rsquared(m, averageSp = FALSE, type = "ols") |> as.data.frame())
write.csv(R2pops, "outputs/R2_populationmodels_nolatentvars.csv")
R2pops$species = rownames(R2pops)

################################################################################

## Predict the model ----

# prepare and format data
formPredDat = as.HMSCdata(X = XData)

# predict population size over time from the model

# original HMSC predict function
pred_all <- predict(m,
                    newdata = formPredDat,
                    type = "response")
# convert to long format, for plotting later
pred_poptrends = pred_all |> as.data.frame() |> 
  mutate(year = time) |>
  tidyr::pivot_longer(cols = -year)

# with quantiles
pred <- predictQ.hmsc(m,
                      newdata = formPredDat,
                      type = "response")

# make a df for plotting later
pred_avgtrend = data.frame(
  year = time,
  mu = pred$overall_mu,
  q025 = pred$overall_q025,
  q975 = pred$overall_q975
)
write.csv(pred_avgtrend, "outputs/pred_avgtrend_nolatentvars.csv", row.names = FALSE)

# plot the model predictions

# convert data to long format for plotting
YData_l = YData |> as.data.frame() |> 
  mutate(year = time) |> pivot_longer(cols = c(-year))

(A = ggplot(data = YData_l) +
    geom_line(aes(x = year, y = value, group = name), lwd = .1, col = "#005b96") +
    geom_ribbon(data = pred_avgtrend, 
                aes(ymin = q025, ymax = q975, x = year),
                alpha = .3, fill = "#6497b1") +
    geom_line(data = pred_avgtrend, aes(x = year, y = mu), lwd = 1, col = "#03396c") +
    labs(y = "Biomass", x = "Year", title = "Observations") +
    coord_cartesian(ylim = c(-3,4)))


(B = ggplot() +
    geom_ribbon(data = pred_avgtrend, 
                aes(ymin = q025, ymax = q975, x = year),
                alpha = .3, fill = "#6497b1") +
    geom_line(data = pred_poptrends,
              aes(x = year, y = value, group = name), col = "#005b96", lwd = .1) +
    geom_line(data = pred_avgtrend, aes(x = year, y = mu), lwd = 1, col = "#03396c") +
    labs(x = "Year", y = "Biomass (scaled)", title = "Predictions") +
    ggpubr::theme_pubr() +
    coord_cartesian(ylim = c(-3,4))) # %>% plotly::ggplotly()

A + B + plot_annotation(tag_levels = "a")
ggsave("figures/hmsc_prediction_nolatentvars.png", width = 10.22, height = 5.22)


################################################################################

## Get parameter estimates ----

# extract the parameters and format them into an `mcmc` object.
coefs = coef(m)

# get variance of the parameter estimates for each population
poptrend_var = apply(m$results$estimation$paramX, 1:2, var) |> as.data.frame()

# prepare to plot population temporal slopes
temp = as.data.frame(coefs$paramX)
temp$pop = rownames(temp)
temp$pop = gsub("_", " ", temp$pop) |> stringr::str_to_sentence()
temp$pop <- factor(temp$pop, levels = temp$pop[order(temp$x1)])
temp <- cbind(temp, poptrend_var$x1)
colnames(temp)[4] = "x1_var"
write.csv(temp, "outputs/poptrends_nolatentvars.csv", row.names = FALSE)

(plot_poptrends = 
    ggplot(data = temp) +
    geom_segment(aes(y = pop, yend = pop, 
                     x = x1 - sqrt(x1_var), xend = x1 + sqrt(x1_var))) +
    geom_point(aes(y = pop, x = x1, fill = x1), size = 4, pch = 21) +
    labs(y = "", x = "Temporal slope estimates (\u03b1)", fill = "\u03b1") +
    geom_vline(xintercept = 0, lwd = .3, lty = 2) +
    geom_vline(xintercept = coefs$meansParamX[2]) +
    coord_cartesian(xlim = c(-1, 1))) +
  scale_fill_distiller(palette = "Spectral", limits = c(-1, 1), direction = 1) +
  theme_pubr() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major.x = element_line())
# ggsave("figures/hmsc_populationslopes.png", width = 8.6, height = 7)


# plot density of parameter estimates for the whole community ----
z1 <- apply(m$results$estimation$paramX, 3L, c)
z2 <- expand.grid(dimnames(m$results$estimation$paramX)[1:2])
coef_alliterations = data.frame(z2, z1) |> 
  dplyr::filter(Var2 == "x1") |>
  tidyr::pivot_longer(cols = -c(Var1, Var2))

(plot_trenddensity = ggdensity(coef_alliterations, x = "value", 
                               add = "mean", rug = FALSE, add.params = list(linetype = 1)) +
    labs(x = "Temporal slope estimates (\u03b1)", y = "Density") +
    geom_vline(xintercept = 0, lwd = .3, lty = 2) +
    theme(legend.position = "none", panel.grid.major.x = element_line())) +
  coord_cartesian(ylim = c(0, 1.4), xlim = c(-1, 1)) +
  annotate("text", x = coefs$meansParamX[2]-0.3, y = 1.35, 
           label = paste0("\u03bc = ", round(coefs$meansParamX[2], digits = 3),
                          "\n \u03c3\u00b2 = ", round(coefs$varX[2,2], digits = 3)))

# arrange the plot panels and save
(plot_trenddensity + 
    coord_cartesian(xlim = c(-1,1), ylim = c(0, 1.8)) + 
    theme(axis.title.y = element_text(vjust = -65)) +
    annotate("text", x = coefs$meansParamX[2]-0.2, y = 1.5, 
             label = paste0("\u03bc = ", round(coefs$meansParamX[2], digits = 3),
                            "\n \u03c3\u00b2 = ", round(coefs$varX[2,2], digits = 3)))) / 
  (plot_poptrends + coord_cartesian(xlim = c(-1,1)) +
     theme(axis.text.y = element_text(face = "italic"),
           panel.grid.major.x = element_line()) + 
     scale_fill_distiller(palette = "Spectral", limits = c(-1, 1), direction = 1)) + 
  plot_layout(heights = (c(1,3))) +
  plot_annotation(tag_levels = "a")
ggsave("figures/distribution_poptrends_nolatentvars.png", width = 8.23, height = 8.5)


################################################################################

## Summarise the model with population and community indices ----

(i_pop = coefs$meansParamX[2])
(i_comm = coefs$varX[2,2])

################################################################################

## Plot model predictions versus the LPI for comparison ----

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

# rlpi_temp = rlpi_poptrends[,grep("X", colnames(rlpi_poptrends))]
# 
# # calculate an index trend for each population
# rlpi_poplpi = rlpi_temp*0 # get an empty matrix ready
# for(i in 1:nrow(rlpi_temp)){
#   rlpi_poplpi[i,1] = 1
#   for(t in 2:ncol(rlpi_temp)){
#     rlpi_poplpi[i,t] = rlpi_poplpi[i,t-1]*10^rlpi_temp[i,t]
#   }
# }

# # create long version for plotting
# rlpi_poplpi$species = rownames(rlpi_poplpi)
# rlpi_poplpi <- pivot_longer(rlpi_poplpi, cols = c(-species)) |>
#   rename("time" = "name")
# rlpi_poplpi$time = readr::parse_number(rlpi_poplpi$time)

(A = ggplot() +
    geom_ribbon(data = rlpi_results, 
                aes(ymin = CI_low, ymax = CI_high, x = time), alpha = .3, fill = "#6497b1") + 
    geom_line(data = rlpi_results, aes(x = time, y = LPI_final)) +
    geom_hline(yintercept = 1, lty = 2, lwd = .5) +
    labs(title = "Living Planet Index", x = "Year", y = "LPI") +
    coord_cartesian(ylim = c(-0.2, 1.8)))

(B = ggplot(
  data = pred_avgtrend) +
    geom_ribbon(aes(ymin = q025, ymax = q975, x = year), 
                alpha = .3, fill = "#6497b1") + 
    geom_line(aes(x = year, y = mu)) +
    geom_hline(yintercept = 0, lty = 2, lwd = .5) +
    labs(title = "Hierarchical model", x = "Year", y = "Average biomass") +
    coord_cartesian(ylim = c(-3, 2)))

# (C = ggplot() +
#     geom_line(data = rlpi_poptrends, aes(x = time, y = value, group = species),
#               lwd = .1) +
#     geom_line(data = rlpi_avgtrend, aes(x = time, y = value))
# )

B + A + plot_annotation(tag_levels = "a")
ggsave("figures/compare_lpi_hmsc_nolatentvars.png", width = 11.7, height = 5.16)


################################################################################

## Curiosities ----

# # import species traits to explore how trends differ by group
# Heike_Traits <- readr::read_csv("~/Documents/GitHub/groundfish-data-analysis/data/Heike_Traits.csv")
# temp$DFO_clean_2014_name = rownames(temp)
# temp_traits = left_join(temp, Heike_Traits)
# 
# # reorder corrplot labels by vertical position in the water column
# verticalposition_siteMean = siteMean[order(temp_traits$vertical.position, decreasing = TRUE),
#                                      order(temp_traits$vertical.position, decreasing = TRUE)]
# corrplot::corrplot(verticalposition_siteMean, type = "lower",
#                    method = "color", 
#                    #tl.cex = 2.5, cl.cex = 3, 
#                    tl.col = "black", 
#                    lowCI.mat = siteLower, uppCI.mat = siteUpper,
#                    font = 3)
# 
# # reorder corrplot labels by trophic level
# trophic_siteMean = siteMean[order(temp_traits$trophic.level, decreasing = TRUE),
#                             order(temp_traits$trophic.level, decreasing = TRUE)]
# corrplot::corrplot(trophic_siteMean, type = "lower",
#                    method = "color", 
#                    #tl.cex = 2.5, cl.cex = 3, 
#                    tl.col = "black", 
#                    lowCI.mat = siteLower, uppCI.mat = siteUpper,
#                    font = 3)
# 
# # plot temporal slopes by grouping
# boxplot(temp_traits$x1 ~ temp_traits$trophic.group)
# boxplot(temp_traits$x1 ~ temp_traits$vertical.position)
# 
# # plot the relationship between temporal slope and body size
# plot(log(temp_traits$WeightData) ~ temp_traits$x1)
# tiny_lm = lm(log(temp_traits$WeightData) ~ temp_traits$x1)
# summary(tiny_lm)
# 
# plot(log(temp_traits$length) ~ temp_traits$x1)
# tiny_lm = lm(log(temp_traits$length) ~ temp_traits$x1)
# summary(tiny_lm)
