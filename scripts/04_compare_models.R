# Script to compare the models

# loading the data ----

biomass_l <- readRDS("outputs/biomass_l.rds")


## comparing the variation explained ----
# (i.e. which one has less residual variance)

## loading the models

linear_multi <- readRDS("outputs/linear_multipop.RDS")
linear_hier <- readRDS("outputs/linear_hierarchical.RDS")

# biomass ~ time
Rsquared(linear_multi)
# + a latent variable for year
Rsquared(linear_hier) # major improvement!!


## comparing the population trends ----

## loading the slope summaries
#poptrends_single <- readRDS("outputs/linear_singlepop/posterior.rds") ---- to fix
poptrends_multi = readRDS("outputs/linear_multipop_population_trends.rds")
poptrends_hier = readRDS("outputs/linear_hierarchical_population_trends.rds")

# slope_sd = poptrends_single$beta_year |> apply(1, sd, na.rm = TRUE) ---- to fix
# slope_mu = poptrends_single$beta_year |> mean(na.rm = TRUE) ---- to fix

# join together into a long data frame
poptrends_multi$model = "multivariate"
poptrends_hier$model = "hierarchical"
population_trends = rbind(poptrends_multi, poptrends_hier) |> as.data.frame()

# Create data with breaks to plot "trend zones" (decline, stable, growing)
data_breaks <- data.frame(start = c(-1, -0.1, 0.1),  
                          end = c(-0.1, 0.1, 1),
                          colors = factor(c(1,3,2)),
                          labels = c("decline", "stable", "growing"))

(plot_poptrends = 
    ggplot(data = population_trends) +
    # geom_rect(data = data_breaks,
    #           aes(xmin = start,
    #               xmax = end,
    #               ymin = - Inf,
    #               ymax = Inf,
    #               fill = colors), alpha = .2) +
    geom_segment(aes(x = x1 - x1_SD,
                     xend = x1 + x1_SD,
                     y = pop, yend = pop, col = model)) +
    geom_point(aes(y = pop, x = x1, col = model), size = 3) +
    labs(y = "", x = "temporal slope") +
    geom_vline(xintercept = 0, lwd = .3) +
    coord_cartesian(xlim = c(-1, 1))) +
  scale_color_manual(values = c("#cc5c76", "#1d457f"))

## comparing the model's overall trend prediction ----

## loading the slope summaries
pred_multi = readRDS("outputs/linear_multipop_pred_l.rds")
pred_hier = readRDS("outputs/linear_hierarchical_pred_l.rds")
df_overall_single <- readRDS("outputs/linear_singlepop/avg_trend.rds")
df_overall_multi = readRDS("outputs/linear_multipop_df_overall.RDS")
df_overall_hier = readRDS("outputs/linear_hierarchical_df_overall.RDS")

ggplot() +
  # raw data
  geom_line(data = biomass_l, aes(x = year, y = biomass, group = pop), lwd = .3, col = "grey85") +
  # # average of single-population models
  # geom_ribbon(data = df_overall_single, aes(x = year, ymin = cilo, ymax = cihi), alpha = .3, fill = "goldenrod1") +
  # geom_line(data = df_overall_single, aes(x = year, y = biomass), linewidth = 1, col = "goldenrod1") +
  # multipopulation model without year effect
  geom_ribbon(data = df_overall_multi, aes(x = year, ymin = q025, ymax = q975), 
              alpha = .6, fill = "#cc5c76") +
  geom_line(data = df_overall_multi, aes(x = year, y = mu), 
            linewidth = 1, col = "white") +
  # hierarchical model
  geom_ribbon(data = df_overall_hier, aes(x = year, ymin = q025, ymax = q975), 
              alpha = .6, fill = "#1d457f") +
  geom_line(data = df_overall_hier, aes(x = year, y = mu), 
            linewidth = 1, col = "white") +
  labs(x = "Year", y = "Biomass")


## compare distribution of slope estimates ====

linear_multi$results$estimation$paramX[,2,] |> 
  density() |> 
  plot(col = "#cc5c76", frame = F, main = "Distribution of temporal slopes", ylim = c(0, 1.5))
m = linear_multi
# average slope with CI
paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.025, .5, .975))
abline(v = paramX_overall_mean[2], col = "#cc5c76", lwd = 2)
abline(v = paramX_overall_quantile[1,2], col = "#cc5c76", lty = 2)
abline(v = paramX_overall_quantile[3,2], col = "#cc5c76", lty = 2)

linear_hier$results$estimation$paramX[,2,] |> density() |> lines(col =  "#1d457f")
m = linear_hier
# average slope with CI
paramX_overall_mean = m$results$estimation$paramX |> apply(2, mean) # this is the same as coef$paramX
paramX_overall_sd = m$results$estimation$paramX |> apply(2, sd)
paramX_overall_quantile = m$results$estimation$paramX |> apply(2, quantile, prob = c(.025, .5, .975))
abline(v = paramX_overall_mean[2], col = "#1d457f", lwd = 2)
abline(v = paramX_overall_quantile[1,2], col = "#1d457f", lty = 2)
abline(v = paramX_overall_quantile[3,2], col = "#1d457f", lty = 2)
