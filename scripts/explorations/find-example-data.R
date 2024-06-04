# Script to find poisson-distributed example data from LPD

library(dplyr)
library(tidyr)
library(ggplot2)

lpd = read.csv("data/LPD2022_public.csv")

# cod dataset
#dfo = filter(lpd, Replicate == 0 & Citation == "DFO (2016). DFO Maritimes Research Vessel Trawl Surveys Fish Observations. Saint Andrews Biological Station, OBIS, Digital.")

# filter to remove replicates
lpd = filter(lpd, Replicate == 0 & Units == "Number of individuals")

NAO = lpd %>% filter(M_ocean == "North Atlantic Ocean")

df = NAO %>% pivot_longer(cols = 33:103,
                     names_to = "time", 
                     values_to = "y") 
df$y <- as.numeric(df$y)
df <- df[-which(is.na(df$y)),]
df$time <- gsub("X", "", df$time) %>% as.integer()

df2 = filter(df, Citation == "Theriault, M.-H. (2015). DFO Gulf Region Community Aquatic Monitoring Program (CAMP) Version 4. Bedford Institute of Oceanography. Dartmouth, NS, Canada, OBIS Canada Digital Collections.")
# subset to 2004-2010 because the data is most consistent
df2 = filter(df2, time > 2003 & time < 2011)

# subset to essential columns
df2 = subset(df2, select = c(ID, Binomial, time, y))

ggplot(data = df2) +
  geom_line(aes(x = time, y = log10(y+1), col = Binomial, group = ID))

# save the dataset
saveRDS(df2, "data/exampledata-CAMP2015.rds")
