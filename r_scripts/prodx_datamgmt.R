# Livestock and Crop Production Data Mgmt

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)

data <- read_csv('data/joindatav2.csv')

# helper function
output_fun <- function(df) {
  
  total <- rowSums(df)
  
  if (is.na(sum(total))) {
    
    for (j in 1:ncol(df)) {
      
      naidx <- which(is.na(df[, j]))
      df[naidx, j] <- 0
      
    }
    
    print('NAs coerced to 0 value')
    out <- list(total = rowSums(df),
                df = df)
    return(out)
    
  } else {
    
    print('No NA conversion needed')
    out <- list(total = total,
                df = df)
    return(out)
    
  }
  
}


### Construct production data -------------------------------------------------

### Outputs ###

# Livestock
liv_outputs <- data %>%
  select(matches('Consumed|Eaten|Sold|Milk|Eggs'))

# Crops
crop_outputs <- data %>%
  select(contains('Produced')) %>%
  select(-matches('Milk|Eggs'))

# sorghum and pulses are characters -- fix this
chr_sorgh <- c('500GM', '08k')
chr_pulse <- c('0.5.', '6k', '1KG', '0/')
crop_outputs$HomeSorghumProduced[crop_outputs$HomeSorghumProduced %in% chr_sorgh] <- c(0.5, 8)
crop_outputs$HomePulsesProduced[crop_outputs$HomePulsesProduced %in% chr_pulse] <- c(0.5, 6, 1, 0)
crop_outputs <- sapply(crop_outputs, as.numeric)

# Total ag output measures
liv_outputs_total <- tibble(liv_out_total = output_fun(liv_outputs)[['total']])
crop_outputs_total <- tibble(crop_out_total = output_fun(crop_outputs)[['total']])

### Inputs ###

Land <- data['TotalAcres']
Labor <- data[, 'TotalHHMembers'] - rowSums(data[, grep('Under3', names(data))])

Capital_gen <- data %>%
  select(matches('Usable|Income|Mud|Stone'))
Capital_crop <- data %>%
  select(contains('Implements'))

TotalLivIncome <- data %>%
  select(contains('SaleValue')) %>%
  rowSums() %>%
  as.data.frame() %>%
  rename(., TotalLivInc = .)
Liv_exp <- data %>% 
  select(contains('Exp')) %>% 
  sapply(function(j) replace_na(j, 0))
Liv_treat_exp <- data %>% 
  select(matches('TreatCost|otherCost')) %>% 
  sapply(function(j) replace_na(j, 0))
TotalLivExp <- tibble(TotalLivExp = rowSums(Liv_exp))
TotalLivTreatExp <- tibble(TotalLivTreatExp = rowSums(Liv_treat_exp))  

Edu <- data %>%
  select(contains('Educ')) %>%
  sapply(function(j) ifelse(j %in% c(99, 77, 5), 0, j))
max_Edu <- apply(Edu, 1, max)
bin_Edu <- tibble(no_edu = ifelse(max_Edu == 1, 1, 0),
                  prim_edu = ifelse(max_Edu == 2, 1, 0),
                  sec_edu = ifelse(max_Edu == 3, 1, 0))

liv_inputs <- output_fun(bind_cols(Land, Labor, Capital_gen,
                                   TotalLivIncome, TotalLivExp, TotalLivTreatExp,
                                   bin_Edu))[['df']]

crop_inputs <- output_fun(bind_cols(Land, Labor, Capital_gen,
                                    Capital_crop, TotalLivIncome, bin_Edu))[['df']]

# village binary
# use village 2 as the base
vil <- unique(data[['VillageID.x']])
vil_var <- paste('vil', vil, sep = '')
bin_vil <- matrix(0, nrow = nrow(data), ncol = length(vil_var))
for (j in 1:(length(vil_var) - 1)) {
  
  bin_vil[, j] <- ifelse(data['VillageID.x'] == vil[j], 1, 0)
  
}
colnames(bin_vil) <- vil_var
bin_vil <- as_tibble(bin_vil[, -dim(bin_vil)[2]])

# combine all production data
liv_prodx <- bind_cols(liv_outputs_total, liv_inputs, bin_vil)
crop_prodx <- bind_cols(crop_outputs_total, crop_inputs, bin_vil)
write_csv(liv_prodx, 'data/liv_prodx.csv')
write_csv(crop_prodx, 'data/crop_prodx.csv')


