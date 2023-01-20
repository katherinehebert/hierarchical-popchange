# Script to build hierarchical linear model on fish population trends

library(dplyr)
library(HMSC)
library(ggplot2)

#### import & set up data ####

# from https://github.com/eric-pedersen/groundfish-data-analysis
# note: these are mean biomasses over the whole sampled area with an SE
# not actually raw abundances. But for this purpose, should be ok
load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")

time = rownames(Year_Geom_Means_all) %>% as.numeric()

# get number of populations & length of time series
npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

#### build HMSC model ####

# format data
XData <- as.matrix(time)
YData <- as.matrix(Year_Geom_Means_all)
randomEff <- as.factor(1:tsl) # assign a location to pairings of pops

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
