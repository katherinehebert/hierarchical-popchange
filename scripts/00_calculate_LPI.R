# Script to calculate LPI on the data used to build models

library(rlpi)
library(ggplot2)
library(dplyr)
library(tidyr)

#### import & set up data ####

# from https://github.com/eric-pedersen/groundfish-data-analysis
# note: these are mean biomasses over the whole sampled area with an SE
# not actually raw abundances. But for this purpose, should be ok
load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")

time = rownames(Year_Geom_Means_all) %>% as.numeric()

# get number of populations & length of time series
npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

#### manipulate into rlpi format ####

df = Year_Geom_Means_all
df$year = time
df = df %>%
  pivot_longer(cols = 1:30, 
               names_to = "Binomial", 
               values_to = "popvalue") 
# assign population IDs
ids = data.frame(
  "Binomial" = unique(df$Binomial),
  "ID" = 1:length(unique(df$Binomial))
)
df = full_join(df, ids)
saveRDS(df, "data/fish_rlpiformat.rds")

#### run rlpi ####

temp = t(Year_Geom_Means_all) %>% as.data.frame()
colnames(temp) <- paste0("X", colnames(temp))
temp$Binomial <- rownames(temp)
temp = left_join(temp, subset(df, select = c("Binomial", "ID")))

# create infile
create_infile(pop_data_source = temp, 
              index_vector = TRUE, 
              start_col_name = "X1981", 
              end_col_name = "X2013")
df <- subset(df, select = c("Binomial", "ID", "year", "popvalue"))
write.table(df, "default_infile_pops.txt", row.names = FALSE)

# make LPI
lpi <- LPIMain("default_infile_infile.txt", use_weightings = 1, 
               REF_YEAR = 1981, PLOT_MAX = 2013, BOOT_STRAP_SIZE = 1000)

ggplot_lpi(lpi)  + theme(axis.text.x = element_text(angle = 0, hjust = .5))
ggsave("figures/rlpi_results.png", width = 6.45, height = 3.75)

#### save output ####
saveRDS(lpi, "outputs/rlpi_results.rds")
