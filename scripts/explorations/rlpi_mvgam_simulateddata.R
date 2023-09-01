# Script to calculate LPI on the data used to build models

library(rlpi)
library(ggplot2)
library(dplyr)
library(tidyr)

# Simulate dataset
set.seed(2)
dat <- sim_mvgam(T = 30,
                 n_series = 10,
                 seasonality = "hierarchical",
                 n_lv = 1,
                 trend_model = 'AR1', # each time step depends on the previous one
                 family = 'poisson',
                 mu_obs = 10,
                 freq = 10,
                 trend_rel = 1,
                 train_prop = 0.7)

dat_train = dat$data_train
dat_test = dat$data_test

ggplot(dat_train) +
  geom_line(aes(x = time, y = y, group = series), lwd = .5)

# get number of populations & length of time series
npops <- 10
tsl <- 30

#### manipulate into rlpi format ####

df = subset(dat_train, select = -c(year, season))
df = rename(df,
            "Binomial" = "series",
            "popvalue" = "y",
            "year" = "time")

# assign population IDs
ids = data.frame(
  "Binomial" = unique(df$Binomial),
  "ID" = 1:length(unique(df$Binomial))
)
df = full_join(df, ids)
saveRDS(df, "data/mvgam_simulateddata_rlpiformat.rds")

#### run rlpi ####
df_wide = pivot_wider(df, id_cols = c(year), names_from = Binomial, values_from = popvalue)
df_wide = as.matrix(df_wide) %>% t()
colnames(df_wide) = df_wide[1,]
df_wide = df_wide[-1,]
df_wide = as.data.frame(df_wide)
df_wide$Binomial = rownames(df_wide)
df_wide = right_join(distinct(subset(df, select = c("ID", "Binomial"))), df_wide)

# create infile
create_infile(pop_data_source = df_wide, 
              index_vector = TRUE, 
              start_col_name = "1", 
              end_col_name = "21")
df <- subset(df, select = c("Binomial", "ID", "year", "popvalue"))
write.table(df, "default_infile_pops.txt", row.names = FALSE)

# make LPI
lpi <- LPIMain("default_infile_infile.txt", use_weightings = 1, 
               REF_YEAR = 1, PLOT_MAX = 21, BOOT_STRAP_SIZE = 1000)

ggplot_lpi(lpi) + 
  theme(axis.text.x = element_text(angle = 0, hjust = .5))
ggsave("figures/rlpi_simulateddata.png", width = 6.45, height = 3.75)

#### save output ####
saveRDS(lpi, "outputs/rlpi_simulateddata_results.rds")
