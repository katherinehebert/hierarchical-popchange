# Script to compare the hierarchical model output with the Living Planet Index

# libraries ----

library(here)
library(dplyr)
library(tidyr)
library(mvgam)
library(tidybayes)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggdist)
library(gganimate)
theme_set(theme_pubr() +
            theme(panel.grid.major.x = element_line()))

set.seed(12)

# read the previous data objects and the model object --------------------------

biomass_scaling = readRDS("outputs/biomass_scaling.rds")
mod1 = readRDS(paste0("outputs/gam_hierarchical_gp.rds"))
biomass = biomass_scaling$biomass
derivs = readRDS("outputs/derivs.rds")
derivs_pops_df = readRDS("outputs/derivs_pops_df.rds")
temporal_trend_yearly = readRDS("outputs/temporal_trend_yearly.rds")
avg_deriv_trend = readRDS("outputs/gam_hierarchical_df_overall_deriv.rds")

## Plot the LPI for the assemblage ---------------------------------------------

rlpi_poptrends <- readr::read_csv("default_infile_pops_lambda.csv", 
                                  col_select = -c(1,3))
rlpi_results <- readRDS("outputs/rlpi_results.rds")
rlpi_results$time <- as.integer(rownames(rlpi_results))
rlpi_results <- rlpi_results[-nrow(rlpi_results),] # remove last value

(rlpi_index = ggplot() +
    geom_ribbon(data = rlpi_results, 
                aes(ymin = CI_low, ymax = CI_high, x = time), alpha = .3, fill = "#6497b1") + 
    geom_line(data = rlpi_results, aes(x = time, y = LPI_final), col = "#03396c") +
    geom_hline(yintercept = 1, lwd = .3) +
    labs(y = "Living Planet Index", x = "") +
    coord_cartesian(ylim = c(0, 2))+
    theme(panel.grid.major.x = element_line()))


################################################################################
## Calculate the Living Planet Index with the growth rates from the hierarchical model

# input:
# dt: population growth rates
# cilo: lower confidence interval of dt (from model or bootstrapping)
# cihi: upper confidence interval of dt (from model or bootstrapping)

# function to calculate LPI value 
calclpi <- function(dt, dt_cilo, dt_cihi){
  
  LPI = data.frame(lpi = NA, cilo = NA, cihi = NA)
  LPI[1,] = 1 # initial value is 1 (no standard error)
  
  for(i in 2:length(dt)){
    LPI[i,"lpi"] <- LPI[i-1,"lpi"]*10^dt[i]
    # confidence intervals
    LPI[i,"cilo"] <- LPI[i-1,"lpi"]*10^(dt_cilo[i])
    LPI[i,"cihi"] <- LPI[i-1,"lpi"]*10^(dt_cihi[i])
  }
  return(LPI)
}

# calculate the LPI for the whole assemblage, based on the HGAM results
dgam_lpi =  calclpi(avg_deriv_trend$avg_trend, avg_deriv_trend$cilo, avg_deriv_trend$cihi)
dgam_lpi$time = rlpi_results$time

# Plot the two index trends ####################################################

palette_col = c("DGAM" = "maroon", "LPI" = "#03396c")
(index_plot = ggplot() +
    geom_ribbon(data = rlpi_results, 
                 aes(ymin = CI_low, ymax = CI_high, x = time, fill = "LPI"), 
                 alpha = .2) +
    geom_line(data = rlpi_results,
              aes(y = LPI_final, x = time, col = "LPI")) +
    geom_ribbon(data = dgam_lpi, 
                aes(ymin = cilo, ymax = cihi, x = time, fill = "DGAM"), 
                alpha = .2) +
    geom_line(data = dgam_lpi,
              aes(y = lpi, x = time, col = "DGAM")) +
    geom_hline(yintercept = 1, lwd = .3) +
    labs(x = "",
         y = "Living Planet Index",
         col = "Methodology", fill = "Methodology") +
    scale_color_manual(values = palette_col) +
    scale_fill_manual(values = palette_col) +
    coord_cartesian(ylim = c(0, 2)) +
    theme(panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line()))

ggsave("figures/compare_lpi_mvgam.png", width = 5.8, height = 5.2)


## Alternatively, calculate the LPI based on the time smoother -----------------
## which represents the common temporal trend across species

png("figures/time_smooth_deriv.png")
plot_mvgam_smooth(mod1, derivatives = FALSE)
dev.off()
