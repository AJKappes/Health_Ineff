# Clean PBASS Data for SFA analysis

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)

### data import -------------------------------------------------------------
dfiles <- paste('data/', list.files('data/'), sep = '')
dlist <- vector(mode = 'list', length = length(dfiles))
names(dlist) <- stringr::str_sub(dfiles, 6, -5)
for (d in 1:length(dlist)) {
  
  dlist[[d]] <- read_csv(dfiles[d])
  
}

### household matching over time ---------------------------------------------
m_y <- unique(dlist$liv_inv[['IntDate']])
m_ylist <- vector(mode = 'list', length = length(m_y))
names(m_ylist) <- m_y
for (key in m_y) {
  
  data <- list(dlist$crop_land[dlist$crop_land['IntDate'] == key, ],
               dlist$hh_demogs[dlist$hh_demogs['IntDate'] == key, ], 
               dlist$liv_inc[dlist$liv_inc['IntDate'] == key, ],
               dlist$liv_inv[dlist$liv_inv['IntDate'] == key, ],
               dlist$liv_exp2[dlist$liv_exp2['IntDate'] == key, ],
               dlist$milkcrop_prod[dlist$milkcrop_prod['IntDate'] == key, ],
               dlist$prod_capital[dlist$prod_capital['IntDate'] == key, ]) %>% 
    lapply(function(d) distinct(d, HousehldID, .keep_all = TRUE))
  
  # reduce dfs in list to one df with variable join
  m_ylist[[key]] <- data %>% reduce(inner_join, by = 'HousehldID')
  
}

# combine and write data ----------------------------------------------------
data <- bind_rows(m_ylist)
write_csv(data, path = 'data/joindatav3.csv')

