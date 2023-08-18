# make a model with the dataset at each time step?

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

model_poptrend <- function(year){
  
  Year_Geom_Means_all[1:which(rownames(Year_Geom_Means_all) == year),]
  
  time = rownames(Year_Geom_Means_all) %>% as.numeric()
  # get number of populations & length of time series
  npops <- ncol(Year_Geom_Means_all)
  tsl <- nrow(Year_Geom_Means_all)
  
  ## build the model with HMSC ----
  
  # format data
  
  XData <- as.matrix(time)
  
  YData <- as.matrix(Year_Geom_Means_all) 
  
  # instead of centering on the mean, center on the baseline biomass
  for(i in 1:ncol(YData)){
    YData[,i] = YData[,i] - YData[1,i]
  }
  
  # scale
  YData <- apply(YData, 2, scale, center = FALSE)
  
  # # log-transform abundance
  # YData <- apply(YData, 2, function(x) log(x+0.001))
  # hist(YData)
  
  # set random effect per year
  randomEff <- as.factor(1:tsl) 
  
  # prepare data
  formDat <- as.HMSCdata(X = XData, 
                         Y = YData, 
                         Random = randomEff, 
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
  saveRDS(m, paste0("models/hmsc_", year, ".RDS"))
  
  return(m)
}

# apply model to each year
m_year <- lapply(time[-c(1:5)], model_poptrend)

# get coefficient estimates
m_coefs <- lapply(m_year, coef)

# pull out the mus
mus = list()
for(i in 1:length(m_coefs)){
  mus[[i]] = m_coefs[[i]]$meansParamX[2]
}

mus |> unlist() |> plot(type = "l", ylim = c(-1,1))

