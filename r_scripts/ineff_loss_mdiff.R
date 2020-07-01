# Inefficiency - Livestock Health Loss Analysis

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)


### Data and helpers -------------------------------------------------------------

# aggregate data
dat <- read_csv('data/ineff_loss_sym4.csv')
dat$edu <- ifelse(dat$no_edu == 0, 1, 0)
dat$respskin <- ifelse(dat$RespDisorders == 1 | dat$SkinDisorders == 1,
                       1, 0)
dat$acctreat <- ifelse(dat$TotalLivTreatExp > 0, 1, 0)

# inefficiency difference function
get_ineff_diff <- function(inefflist) {
  
  out <- sapply(1:length(inefflist), function(v)
    inefflist[[v]][2] - inefflist[[v]][1])
  
  names(out) <- names(inefflist)
  
  return(out)
  
}

# names
comgraz_names <- c('noncommunal',
                   'communal',
                   'ineff_diff')

health_names <- c('without_sympt',
                  'with_sympt',
                  'ineff_diff')


### Aggreagate level stats -------------------------------------------------------

# general illness inspection
gen_ill_agg <- sapply(0:1, function(i)
  mean(dat[dat$GeneralIllness == i, ]$p_ineff))

# edu gen ill
#   having > edu, not sick, sick
#   not having edu, not sick, sick
#   inefficiencies all lower when having edu
gen_ill_agg_primedu <- lapply(0:1, function(i)
  sapply(0:1, function(j)
    mean(dat[dat['prim_edu'] == i &
               dat['sec_edu'] == 0 &
               dat['GeneralIllness'] == j, ]$p_ineff)))

gen_ill_agg_secedu <- lapply(0:1, function(i)
  sapply(0:1, function(j)
    mean(dat[dat['sec_edu'] == i &
               dat['no_edu'] == 0 &
               dat['GeneralIllness'] == j, ]$p_ineff)))

gen_ill_agg_noedu <- lapply(0:1, function(i)
  sapply(0:1, function(j)
    mean(dat[dat['edu'] == i &
               dat['GeneralIllness'] == j, ]$p_ineff)))

# access to livestock treatment
respskin_agg_acctreat <- lapply(0:1, function(i)
  sapply(0:1, function(j)
    mean(dat[dat['acctreat'] == i &
               dat['respskin'] == j, ]$p_ineff)))

# disorder inspection
disorders <- names(dat)[grep('Disorders', names(dat))]
dis_dlist_agg <- lapply(disorders, function(d)
  sapply(0:1, function(i)
    mean(dat[dat[d] == i, ]$p_ineff)))
names(dis_dlist_agg) <- disorders

# shared grazing - 1 shared, 0 not shared
dat$comgraz_bin <- ifelse(dat$UnsharedLand < 1, 1, 0)
comgraz_agg <- sapply(0:1, function(i)
  mean(dat[dat$comgraz_bin == i, ]$p_ineff))

genill_comgraz_agg <- sapply(0:1, function(i)
  mean(dat[dat['GeneralIllness'] == i &
         dat['comgraz_bin'] == 1, ]$p_ineff))

# resp given communal grazing
resp_comgraz_agg <- sapply(0:1, function(i)
  mean(dat[dat['RespDisorders'] == i &
             dat['comgraz_bin'] == 1, ]$p_ineff))

# skin given communal grazing
skin_comgraz_agg <- sapply(0:1, function(i)
  mean(dat[dat['SkinDisorders'] == i &
             dat['comgraz_bin'] == 1, ]$p_ineff))

# resp or skin given communal grazing
respskin_comgraz_agg <- sapply(0:1, function(i)
  mean(dat[dat['RespDisorders'] == i |
             dat['SkinDisorders'] == i &
             dat['comgraz_bin'] == 1, ]$p_ineff))


### Village level stats ----------------------------------------------------------

# subset by village
vill_id <- names(dat)[grep('vil', names(dat))]
vill_dlist <- lapply(vill_id, function(v) dat[dat[v] == 1, ])
names(vill_dlist) <- vill_id

# general illness
gen_ill_vill <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['GeneralIllness'] == i, ]$p_ineff)))

# communal grazing
vill_comgraz <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['comgraz_bin'] == i, ]$p_ineff)))

# treatment occurance
vill_treat <- lapply(vill_dlist, function(v)
  c(mean(v$TotalLivTreatExp), mean(v$AnimalTreated),
    mean(v$TotalLivTreatExp)*mean(v$AnimalTreated)))

vill_treat_bin <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['AnimalTreated'] == i, ]$p_ineff)))

treat_prop <- sapply(1:length(vill_treat), function(v)
  vill_treat[[v]][2])

treat_avgexp <- sapply(1:length(vill_treat), function(v)
  vill_treat[[v]][3])

# respiratory disorder
vill_resp <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['RespDisorders'] == i, ]$p_ineff)))

vill_resp_prop <- sapply(vill_dlist, function(v)
  mean(v$RespDisorders))

# skin disorder
vill_skin <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['SkinDisorders'] == i, ]$p_ineff)))

vill_skin_prop <- sapply(vill_dlist, function(v)
  mean(v$SkinDisorders))

# general illness given communal grazing
vill_genill_comgraz <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['GeneralIllness'] == i &
             v['comgraz_bin'] == 1, ]$p_ineff)))

# respiratory disease given communal grazing
vill_resp_comgraz <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['RespDisorders'] == i &
             v['comgraz_bin'] == 1, ]$p_ineff)))

# skin disease given communal grazing
vill_skin_comgraz <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['SkinDisorders'] == i &
             v['comgraz_bin'] == 1, ]$p_ineff)))

# resp or skin given communal grazing
vill_respskin_comgraz <- lapply(vill_dlist, function(v)
  sapply(0:1, function(i)
    mean(v[v['RespDisorders'] == i |
             v['SkinDisorders'] == i &
             v['comgraz_bin'] == 1, ]$p_ineff)))

# ineff diffs
agg_dis_ineff_diff <- get_ineff_diff(dis_dlist_agg)
gen_ill_vill_diff <- get_ineff_diff(gen_ill_vill)
comgraz_vill_diff <- get_ineff_diff(vill_comgraz)
resp_vill_diff <- get_ineff_diff(vill_resp)
skin_vill_diff <- get_ineff_diff(vill_skin)
genill_comgraz_diff <- get_ineff_diff(vill_genill_comgraz)
resp_comgraz_diff <- get_ineff_diff(vill_resp_comgraz)
skin_comgraz_diff <- get_ineff_diff(vill_skin_comgraz)
respskin_comgraz_diff <- get_ineff_diff(vill_respskin_comgraz)

### Tables ------------------------------------------------------------

# aggregate
agg_cond_table <- rbind(genill_comgraz_agg, resp_comgraz_agg,
                  skin_comgraz_agg, respskin_comgraz_agg)
agg_cond_diff <- agg_cond_table[, 2] - agg_cond_table[, 1]
agg_cond_table <- cbind(agg_cond_table, agg_cond_diff)
colnames(agg_cond_table) <- comgraz_names

agg_dis_diff_table <- bind_rows(dis_dlist_agg,
                                agg_dis_ineff_diff) %>%
  t()
colnames(agg_dis_diff_table) <- health_names

# village
vill_genill_table <- bind_rows(gen_ill_vill,
                               gen_ill_vill_diff) %>% 
  t()
colnames(vill_genill_table) <- c('no_illness',
                                 'illness',
                                 'ineff_diff')

vill_comgraz_table <- bind_rows(vill_comgraz,
                                comgraz_vill_diff) %>% 
  t()
colnames(vill_comgraz_table) <- comgraz_names

vill_resp_table <- bind_rows(vill_resp,
                             resp_vill_diff) %>% 
  t()
colnames(vill_resp_table) <- health_names

vill_skin_table <- bind_rows(vill_skin,
                             skin_vill_diff) %>% 
  t()
colnames(vill_skin_table) <- health_names

vill_gencom_table <- bind_rows(vill_genill_comgraz,
                               genill_comgraz_diff) %>% 
  t()
colnames(vill_gencom_table) <- c('no_illness',
                                 'illness',
                                 'ineff_diff')

vill_respcom_table <- bind_rows(vill_resp_comgraz,
                                resp_comgraz_diff) %>% 
  t()
colnames(vill_respcom_table) <- health_names

vill_skincom_table <- bind_rows(vill_skin_comgraz,
                                skin_comgraz_diff) %>% 
  t()
colnames(vill_skincom_table) <- health_names

vill_respskincom_table <- bind_rows(vill_respskin_comgraz,
                                    respskin_comgraz_diff) %>% 
  t()
colnames(vill_respskincom_table) <- health_names


### Plots --------------------------------------------------------------

# plot(treat_prop, gen_ill_vill_diff)
# plot(treat_prop, genill_comgraz_diff)
# plot(treat_prop, respskin_comgraz_diff)
# plot(treat_prop, resp_comgraz_diff)
# plot(treat_prop, skin_comgraz_diff)


### Testing ------------------------------------------------------------
