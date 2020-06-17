# Animal Health Loss-Inefficiency Proportion Modeling

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)

### Data mgmt ----------------------------------------------------------------
dat_ineff <- read_csv('data/liv_ineff_loss5.csv')
dat_liv_sick <- read_csv('data/livestock_sick.csv')

# match on household ID
m_y <- unique(dat_ineff$IntDate.x)
m_ylist <- vector(mode = 'list', length = length(m_y))
names(m_ylist) <- m_y
for (key in m_y) {
  
  dl <- list(dat_ineff[dat_ineff$IntDate.x == key, ],
             dat_liv_sick[dat_liv_sick$IntDate == key, ]) %>% 
    lapply(function(d) distinct(d, HousehldID, .keep_all = TRUE))
  
  # reduce dfs in list to joined df
  m_ylist[[key]] <- dl %>%
    reduce(inner_join, by = 'HousehldID')
  
}

dat_ineff_sym <- bind_rows(m_ylist)

deathloss <- dat_ineff_sym %>% 
  select(contains('Dead')) %>% 
  rowSums()

dat_ineff_sym$aggloss <- deathloss

# write_csv(dat_ineff_sym, 'data/ineff_loss_sym4.csv')
