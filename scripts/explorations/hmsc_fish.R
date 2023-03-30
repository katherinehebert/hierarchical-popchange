# Script to build hierarchical linear model on fish population trends

library(dplyr)
library(HMSC)
library(ggplot2)
library(tidyr)
theme_set(ggpubr::theme_pubr())

#### import & set up data ####

library(readr)
LPD <- read_csv("~/Documents/GitHub/synchrony/data_raw/LPD2022_public.csv")

# keep only the DFO dataset
df <- LPD[grep("DFO \\(2016\\)", LPD$Citation, value = FALSE),]
# 2016 is also maybe a good candidate

# convert LPD to long format
df <- pivot_longer(df, cols = 33:103,
                   names_to = "year_obs", 
                   values_to = "obs_value") %>%
  subset(select = c(ID, Genus, Species, year_obs, obs_value)) 
df$obs_value <- as.numeric(df$obs_value)
df <- na.omit(df) #NAs are empty years. delete them

# reformat some columns
df$year_obs <- parse_number(df$year_obs)
df$scientific_name <- paste0(df$Genus, "_", df$Species)

# get number of populations & length of time series
time = unique(df$year_obs) %>% as.numeric()
npops <- length(unique(df$ID))
tsl <- length(unique(df$year_obs))

ggplot(df) +
  geom_line(aes(x = year_obs, y = log(obs_value), col = scientific_name, group = ID)) +
  labs(x = "Year", y = "Population size (log)") +
  theme(legend.position = "none")

# convert to wide format
df_l <- df
df_w <- pivot_wider(subset(df, select = -c(Genus, Species, scientific_name)),
                    names_from = ID, 
                    values_from = obs_value)

#### build HMSC model ####

# format data
XData <- as.matrix(df_w$year_obs)
YData <- as.matrix(subset(df_w, select = -c(year_obs)))
randomEff <- as.factor(df_w$scientific_name) # assign a location to pairings of pops

# prepare data
formDat <- as.HMSCdata(X = XData, Y = YData, Random = randomEff)

# run model
m <- hmsc(formDat,
                family = "gaussian",
                niter = 15000,
                nburn = 5000,
                thin = 5,
                verbose = FALSE)
# save model
saveRDS(m, "models/hmsc.RDS")

#### explore the model ####

# extract the parameters and format them into an `mcmc` object.
mcmcMeansParamX <- as.mcmc(m, parameters = "meansParamX")
# check for parameter convergence (visually)
traceplot(mcmcMeansParamX)

# model explanatory power
R2 <- Rsquared(m, averageSp = FALSE)
R2comm <- Rsquared(m, averageSp = TRUE)

# average model per species across iterations
avg <- as.data.frame(apply(m$results$estimation$paramX, 1:2, mean))

# get 95% confidence intervals
ci_lo <- as.data.frame(apply(m$results$estimation$paramX, 1:2, quantile, probs = 0.025))
ci_hi <- as.data.frame(apply(m$results$estimation$paramX, 1:2, quantile, probs = 0.975))

# plot all populations' models
ggplot(data = as.data.frame(XData), aes(x = V1)) +
  geom_abline(data = avg, 
              aes(slope = x1, intercept = Intercept),
              lwd = 0.5) +
  # geom_abline(data = ci_lo,
  #             aes(slope = x1, intercept = Intercept),
  #             lty = 2, col = "blue") +
  # geom_abline(data = ci_hi,
  #             aes(slope = x1, intercept = Intercept),
  #             lty = 2, col = "red") +
  ylim(c(-1,2)) +
  theme_classic()

#### covariance matrix ####

# extract all estimated species-to-species associations matrix
assoMat <- corRandomEff(m)
# Average
siteMean <- apply(assoMat[, , , 1], 1:2, mean)

# plot as a heatmap
heatmap(siteMean, Rowv = NA, Colv = NA)

# build matrix of colours for chordDiagram (red : neg, blue: pos)
siteDrawCol <- matrix(NA, nrow = nrow(siteMean),
                      ncol = ncol(siteMean))
siteDrawCol[which(siteMean < 0, arr.ind=TRUE)]<-"red"
siteDrawCol[which(siteMean > 0, arr.ind=TRUE)]<-"blue"

# chord diagram
circlize::chordDiagram(siteMean, symmetric = TRUE,
                       annotationTrack = c("grid"),
                       grid.col = "grey", 
                       col = siteDrawCol)


#### predict and plot model over time ####

# prepare and format data
formPredDat = as.HMSCdata(X = XData, Random = randomEff)

# predict values
pred <- predict(m,
                newdata = formPredDat,
                type = "response")

#### get mean trend across populations ####

# take the mean N per year
pred_mean = apply(pred, 1, mean)
plot(pred_mean ~ XData, 
     col = "white", xlab = "years", ylab = "N",
     ylim = c(min(pred), max(pred)))
for(i in 1:ncol(pred)) {lines(pred[,i] ~ XData, 
                              col = Colours[i], lwd = .5)}
lines(pred_mean ~ XData, lwd = 2)


#### calculate rate of change ####

# calculate log-ratio population change
calc_dt <- function(N){
  dt = c(1) # initial value
  for(i in 2:length(N)){
    dt[i] = log10(N[i]/N[i-1])
  }
  return(dt)
}

dt = apply(pred, 2, calc_dt)

dt_bar = apply(dt, 1, mean, na.rm = TRUE)
dt_ci = apply(dt, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
dt_var = apply(dt, 1, var, na.rm = TRUE)

calc_lpi <- function(rate){
  lpi = 1
  for(t in 2:tsl){
    lpi[t] = lpi[t-1]*10^rate[t]
  }
  return(lpi)
}

lpi = calc_lpi(dt_bar)
lpi_pops = apply(dt, 2, calc_lpi)

# plot lpi and all population lpis
matplot(lpi_pops, x = time, 
        type = "l", lwd = 0.2, lty = 1, 
        ylim = c(0, 5))
lines(x = time, y = lpi, lwd = 3)

# plot lpi with CI from the variance
plot(x = time, y = lpi, type = "l", ylim = c(0,2))
lines(x = time, y = lpi - 1.96*sqrt(dt_var), lty = 2)
lines(x = time, y = lpi + 1.96*sqrt(dt_var), lty = 2)

# plot the variance
plot(x = time, y = dt_var, type = "l",
     ylab = "Variance of growth rates")
