# Script to convert R scripts to Rmd

# set working directory to folder where I want the reports
setwd("scripts/")

# knitr::spin("hmsc_fish.R", format = "Rmd", knit = FALSE)

knitr::spin("hgam_fish.R", format = "Rmd", knit = FALSE)
