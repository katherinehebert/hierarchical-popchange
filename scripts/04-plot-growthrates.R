# Script to make figures of population growth rates:
## histogram of all population growth rates over all years
## animated histogram of all population growth rates per year
## dot plot of each species' growth rate over all years (with quantiles)


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
derivs_pops_df = readRDS("outputs/derivs_pops_df.rds")
temporal_trend_yearly = readRDS("outputs/temporal_trend_yearly.rds")
temporal_trend = readRDS("outputs/temporal_trend.rds")
avg_deriv_trend = readRDS("outputs/gam_hierarchical_df_overall_deriv.rds")

# plot a histogram of the average derivatives across all species ---------------

derivs_without1981 = dplyr::filter(derivs_pops_df, year != 1981)

# Plot the distribution of all growth rates, animated per year -----------------

# make a df with yearly means
derivs_without1981$frame = as.character(derivs_without1981$year)
derivs_without1981_summary = derivs_without1981 |>
  group_by(year) |>
  summarise(mu = mean(value),
            lower = mean(value) - sd(value),
            upper = mean(value) + sd(value))
derivs_without1981_summary$frame = as.character(derivs_without1981_summary$year)

(plot_trenddensity_anim = 
    ggplot(data = derivs_without1981) +
    geom_histogram(aes(x = value, fill = after_stat(x)), 
                   col = "black", linewidth = .2, bins = 19) + 
    geom_vline(data = derivs_without1981_summary, aes(xintercept = mu)) +
    geom_vline(data = derivs_without1981_summary, aes(xintercept = lower), lty = 2) +
    geom_vline(data = derivs_without1981_summary, aes(xintercept = upper), lty = 2) +
    theme(panel.grid.major.x = element_line(),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          title = element_text(size = 20),
          legend.position = "none") +
    scale_y_sqrt() +
    labs(x = "Taux de croissance annuel", 
         y = "Fréquence", 
         fill = "",
         title = "Année: {closest_state}") +
    scale_fill_distiller(palette = "RdYlGn", 
                         direction = 1, 
                         limits = c(-.4,.4)) +
    coord_cartesian(xlim = c(-.4, .4)) +
    transition_states(frame, state_length = 50)) 
anim_save("figures/trend_density.gif", plot_trenddensity_anim)

# Plot the average derivative per species with 90% quantiles -------------------

(plot_poptrends = ggplot(data = temporal_trend[order(temporal_trend$mu_deriv),]) +
   geom_segment(aes(x = lower, xend = upper,
                    y = species, yend = species), lwd = .3) +
   geom_point(aes(x = mu_deriv, y = species, fill = mu_deriv), 
              size = 4, pch = 21) +
   geom_vline(xintercept = 0, lwd = .3, lty = 2) +
   labs(x = "Median rate of change", 
        y = "",
        fill = "") +
   scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-.1,.1)) +
   coord_cartesian(xlim = c(-.4, .4)) +
   theme_pubr() +
   theme(axis.text.y = element_text(face = "italic"),
         panel.grid.major.x = element_line()) 
)

################################################################################
# 3-panel plot: histogram, mean trend, coherence -------------------------------
################################################################################

col_lims = max(abs(min(avg_deriv_trend$cilo)), max(avg_deriv_trend$cihi))

## Panel A - Histogram of all population growth rates of the full time series ----
(plot_trenddensity = 
   ggplot(data = derivs_without1981) +
   geom_histogram(aes(x = value, fill = after_stat(x)), 
                  col = "black", linewidth = .2, bins = 19) + 
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE)) +
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE) - sd(derivs_without1981$value, na.rm = TRUE), lty = 2) +
   geom_vline(xintercept = mean(derivs_without1981$value, na.rm = TRUE) + sd(derivs_without1981$value, na.rm = TRUE), lty = 2) +
   theme(panel.grid.major.x = element_line()) +
   scale_y_sqrt() +
   labs(x = expression(alpha), 
        y = "Frequency", 
        fill = "α",
        title = "Community Change Index") +
   scale_fill_distiller(palette = "RdYlGn", 
                        direction = 1, 
                        limits = c(-col_lims, col_lims)) +
   coord_cartesian(xlim = c(-.4, .4))) 
# variance
var(derivs_without1981$value, na.rm = TRUE)

## Panel B - Trend of the distribution of rates of change through time ---------
(avg_derivative_pointplot = 
    ggplot(data = avg_deriv_trend, aes(x = year)) +
    ggpattern::geom_area_pattern(aes(y = cihi),
                                 pattern = "gradient", 
                                 fill = "#ffffbf",
                                 pattern_fill  = "#ffffbf",
                                 pattern_fill2 = "#77c15c") + 
    ggpattern::geom_area_pattern(aes(y = cilo),
                                 pattern = "gradient", 
                                 fill = "#ffffbf",
                                 pattern_fill  = "#f98c59",
                                 pattern_fill2 = "#ffffbf") + 
    geom_line(aes(y = avg_trend, col = avg_trend), lwd = .1) +
    geom_line(aes(y = avg_trend), lwd = .8, col = "black") +
    geom_line(aes(y = cilo), lwd = .1, col = "white") +
    geom_line(aes(y = cihi), lwd = .1, col = "white") +
    geom_hline(yintercept = 0, lwd = .3) +
    scale_color_gradientn(colours = c("#d73027", "#f98c59", "#fee08b", "#ffffbf", "#d9ef8b", "#77c15c", "#1a9850"), 
                          limits = c(-col_lims, col_lims)) +
    labs(x = "",
         title = "Community Abundance Index",
         y = expression(mu[alpha]),
         col = expression(mu[alpha]),
         ) +
    coord_cartesian(ylim = c(-.5, .5)) +
    theme(panel.grid.major = element_line(),
          legend.position = "right") +
    scale_y_continuous(labels = scales::percent))

## Panel C - Trend of the distribution's variance through time -----------------

(plot_coherence = ggplot(data = temporal_trend_yearly) +
    geom_line(aes(y = sd^2, x = as.numeric(year), col = sd), lwd = 1.5) +
    geom_point(aes(y = sd^2, x = as.numeric(year), col = sd), 
               size = .9) +
    labs(title = "Community Stability Index", 
         y = expression(sigma^2),
         x = "",
         col = expression(sigma^2)) +
   colorspace::scale_color_continuous_sequential("YlGnBu",begin = .2) +
    #scale_color_distiller(palette = "YlGnBu", direction = 1, limits = c(0, .17),) +
    #coord_cartesian(ylim = c(0, .18)) +
    theme_pubr() +
    theme(panel.grid.major = element_line())
)

mean(temporal_trend_yearly$sd[2:10])

# Arrange the plot panels and save ---------------------------------------------
((plot_trenddensity + theme(legend.position = "top",
                            legend.key.width = unit(2, "cm"))) /
    ((avg_derivative_pointplot + theme(legend.position = "none"))  + 
    (plot_coherence + theme(legend.position = "none")))) +
  plot_annotation(tag_levels = "a")
ggsave("figures/assemblagevariability.png", width = 7.5, height = 7)


# Plot growth rate versus variance

df = 
  data.frame(
    "year" = as.numeric(temporal_trend_yearly$year),
    "mu" = avg_deriv_trend$avg_trend,
    "var" = temporal_trend_yearly$sd
  )

(plot_path = ggplot(data = df) +
    geom_vline(xintercept = 0, lty = 2, linewidth = .3) +
    geom_hline(yintercept = max(df$var, na.rm = T)/2, lty = 2, linewidth = .3) +
    geom_path(aes(x = mu, y = var), linewidth = .2, col = "black") +
    geom_point(aes(x = mu, y = var, fill = as.numeric(year)), size = 4, pch = 21) +
    geom_point(aes(x = mean(derivs_without1981$value, na.rm = TRUE),
                   y = sd(derivs_without1981$value, na.rm = TRUE)), pch = 8, size = 2) +
    ## text labels
    geom_text(aes(x = mean(derivs_without1981$value, na.rm = TRUE)-0.005,
                  y = sd(derivs_without1981$value, na.rm = TRUE),
                  label = "CCI"), 
              col = "black", hjust = 1) +
    colorspace::scale_fill_continuous_divergingx("PRGn", mid = 1997, rev = F) +
    coord_cartesian(xlim = c(-.08,.08)) +
    scale_y_reverse() +
    labs(x = expression(paste("Community Abundance Index ", (mu[alpha]))),
         y = expression(paste("Community Stability Index ", (sigma^2))),
         fill = "Year") +
    theme(legend.position = "right"))
ggsave("figures/assemblage_path.png", width = 5.16, height = 4.06)

(ggplot(data = df) +
    geom_vline(xintercept = 0, lty = 2, linewidth = .3) +
    geom_hline(yintercept = max(df$var, na.rm = T)/2, lty = 2, linewidth = .3) +
    geom_path(aes(x = mu, y = var), linewidth = .2, col = "black") +
    geom_point(aes(x = mu, y = var, fill = as.numeric(year)), size = 4, pch = 21) +
    geom_point(aes(x = mean(derivs_without1981$value, na.rm = TRUE),
                   y = sd(derivs_without1981$value, na.rm = TRUE)), pch = 8, size = 2) +
    ## text labels
    geom_text(aes(x = mu, y = var+0.004, label = as.character(year)),
              col = "black", hjust = 1) +
    colorspace::scale_fill_continuous_divergingx("PRGn", mid = 1997, rev = F, name = "Year") +
    coord_cartesian(xlim = c(-.08,.08)) +
    scale_y_reverse() )
ggsave("figures/assemblage_path_withyears.png", width = 10, height = 10)



# Arrange the plot panels and save ---------------------------------------------


(
  ((plot_trenddensity + theme(legend.position = "none")) +
     (avg_derivative_pointplot + theme(legend.position = "none"))) /
    ((plot_coherence + theme(legend.position = "none")) +
       (plot_path + theme(legend.position = "inside", 
                          legend.position.inside = c(.93,.55),
                          legend.key.height = unit(.5, "cm"))
       ))
) +
  plot_annotation(tag_levels = "a") #+
  #plot_layout(heights = c(1, 2))
ggsave("figures/community-indicators.png", width = 8.75, height = 8.05)
