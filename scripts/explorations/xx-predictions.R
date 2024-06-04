# Script to plot the predictions from the model 

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

# read the previous data objects and the model object ----

data_train = readRDS(paste0("outputs/data_train.rds"))
biomass_scaling = readRDS("outputs/biomass_scaling.rds")
predictions = readRDS("outputs/gam_hierarchical_pred_l.rds")

YMeans = biomass_scaling$YMeans
baselines = biomass_scaling$baselines

# plot the predictions next to the data ----

data_train = rename(data_train, "sp" = "series")
data_train$sp = gsub("_", " ", data_train$sp) |> stringr::str_to_sentence()
data_train$sp = gsub(" ", "\n", data_train$sp)
predictions$sp = gsub("_", " ", predictions$sp) |> stringr::str_to_sentence()
predictions$sp = gsub(" ", "\n", predictions$sp)
ggplot(data = predictions,
       aes(x = time, group = sp)) +
  geom_ribbon(aes(ymin = cilo, ymax = cihi, fill = sp), alpha = .2) +
  geom_line(aes(y = biomass)) +
  geom_line(data = data_train,
            aes(x = time+1981, y = y), lty = 2) +
  facet_wrap(~sp, scales = "free_y") +
  theme(legend.position = "none", strip.text = element_text(face = "italic")) +
  labs(x = "", y = "Biomass (kg per tow)")
ggsave("figures/mvgam_predictions_perspecies.png", width = 12.8, height = 7.63)


# # rescale the predictions to the response scale ----
# predictions_respscale = predictions |> group_by(sp) |> group_split()
# for(i in 1:length(predictions_respscale)){
#   predictions_respscale[[i]]$biomass = predictions_respscale[[i]]$biomass*YMeans[i]
#   predictions_respscale[[i]]$cilo = predictions_respscale[[i]]$cilo*YMeans[i]
#   predictions_respscale[[i]]$cihi = predictions_respscale[[i]]$cihi*YMeans[i]
# }
# predictions_respscale = bind_rows(predictions_respscale)
# 
# # summarise the predicted change in biomass over the whole time series, per species
# predictions_respscale_summary = predictions_respscale |>
#   group_by(sp) |>
#   summarise("mu" = mean(biomass),
#             "sd" = sd(biomass))
# predictions_respscale_summary$sp = factor(predictions_respscale_summary$sp,
#                                           levels = predictions_respscale_summary$sp[order(predictions_respscale_summary$mu)])
# 
# # # convert predictions back to proportion of the baseline biomass
# # predictions_respscale_perc = predictions |> group_by(sp) |> group_split()
# # for(i in 1:length(predictions_respscale_perc)){
# #   predictions_respscale_perc[[i]]$biomass = ((predictions_respscale_perc[[i]]$biomass*YMeans[i]) + baselines[i])/baselines[i]
# #   predictions_respscale_perc[[i]]$cilo = ((predictions_respscale_perc[[i]]$cilo*YMeans[i]) + baselines[i])/baselines[i]
# #   predictions_respscale_perc[[i]]$cihi = ((predictions_respscale_perc[[i]]$cihi*YMeans[i]) + baselines[i])/baselines[i]
# #   
# # }
# predictions_respscale_perc = bind_rows(predictions_respscale_perc)
# 
# # this may be the spot to mess with ----
# (plot_biomassdiff = ggplot(data = predictions_respscale_perc) +
#    geom_line(aes(y = biomass, x = time, group = sp)) +
#    geom_hline(yintercept = 0, lwd = .3, lty = 2) +
#    labs(y = "Population size relative to the baseline", 
#         x = "",
#         fill = "") +
#    #scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-20,1)) +
#    #coord_cartesian(xlim = c(-.5, .5)) +
#    theme_pubr() +
#    theme(axis.text.y = element_text(face = "italic"),
#          panel.grid.major.x = element_line()) 
# )
