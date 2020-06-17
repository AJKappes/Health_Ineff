# Stochastic Frontier Analysis Modeling

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)
library(ggplot2)
library(cowplot)
library(latex2exp)

liv_prodx <- read_csv('data/liv_prodx_wloss4.csv')
crop_prodx <- read_csv('data/crop_prodx2.csv')

### Data mgmt ---------------------------------------------------------------

# address extreme value error in data input
sub_idx <- function(data, extreme, q) {
  
  var_qs <- map_dbl(data[, extreme], function(j) quantile(j, q))
  idx <- map2(extreme, var_qs, function(j, i) which(data[, j] > i))
  idx <- unique(unlist(idx))
  return(data[-idx, ])
  
}

liv_loss_vars <- names(liv_prodx)[grep('Dead', names(liv_prodx))]

liv_ext_vars <- c('OffFarmNetIncome', 'TotalAcres', 'UsablePhones',
                  'liv_out_total', 'TotalLivExp', 'TotalLivTreatExp',
                  'TotalLivInc')

crop_ext_vars <- c('crop_out_total', 'TotalAcres', 'UsablePhones',
                   'OffFarmNetIncome', 'DraftImplements', 'HandImplements',
                   'TotalLivInc')

liv_reduc_vars <- names(liv_prodx)[grep('liv_out|Acres|HHMem|Exp|Inc|vil',
                                        names(liv_prodx))]

crop_reduc_vars <- names(crop_prodx)[grep('crop_out|Acres|HHMem|Usable|
                                            Mud|Buildings|Implements|Inc',
                                            names(crop_prodx))]

# remove zero production
liv_noidx <- which(liv_prodx['liv_out_total'] == 0)

# livestock and crop production
liv_ext_prodx <- sub_idx(liv_prodx[-liv_noidx, ], liv_ext_vars, 0.9)
liv_loss <- liv_ext_prodx[liv_loss_vars]
liv_comgraz <- liv_ext_prodx['UnsharedLand']
liv_edu <- liv_ext_prodx %>% select(contains('edu'))
liv_dates_hhID <- liv_ext_prodx[, c('IntDate.x', 'HousehldID')]
liv_prodx <- liv_ext_prodx[liv_reduc_vars]
crop_prodx <- sub_idx(crop_prodx, crop_ext_vars, 0.9)[crop_reduc_vars]


### Frontier functions -------------------------------------------------

# get initial model values
init_fun <- function(data) {
  
  dep <- names(data)[1]
  vars <- names(data)[-1]
  
  fx <- as.formula(paste(dep, paste(vars, collapse = '+'), sep = '~'))
  mod <- lm(fx, data = data)
  
  beta_init <- mod$coefficients
  err <- mod$residuals
  sig_init <- sqrt(1/nrow(data)*sum(err^2))
  
  out <- list(beta = beta_init, error = err, sig = sig_init)
  return(out)
  
}


# half normal SFA loglikelihood
sfa_ll <- function(theta, data, indicate){
  
  # data and dimensions
  y <- data[, 1] %>% as.matrix()
  X <- cbind(1, data[, -1]) %>% as.matrix()
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # theta collection of params(beta, sigma, lambda)
  beta <- theta[1:p]
  # sig_u <- theta[p+1]
  # sig_v <- theta[p+2]
  # sig <- (sig_u^2 + sig_v^2)^(1/2)
  # lam <- sig_u/sig_v
  sig <- theta[p+1]
  lam <- theta[p+2]
  
  # likelihood variable, e
  # e = v - u
  # u: firm inefficiency
  # v: firm specific idiosyncratic shock v ~ N(0, sig_v^2)
  
  if (indicate == 'fixed' | indicate == 'hetero') {
    
    # fixed effect setup
    e <- y - X%*%beta
    
    # skewed negative errors => switch sign(e) '-'
    ll <- -n*log(sig) + n/2*log(2/pi) -
      1/2*sum((-e/sig)^2) +
      sum(log(pnorm(e*lam/sig)))
    
  } else if (indicate == 'random') {
    
    # random effects setup
    n_sim <- 100
    ll_vec <- c()
    h_idx <- grep('Inc|vil', names(data))
    h_mus <- c(X[, h_idx]%*%theta[h_idx])
    sd_mus <- sd(h_mus)
    for (i in 1:n_sim) {
  
      w <- rnorm(dim(X)[1], mean = h_mus, sd = sd_mus)
      e <- y - X%*%beta - w
      ll <- -n*log(sig) + n/2*log(2/pi) -
        1/2*sum((e/sig)^2) +
        sum(log(pnorm(-e*lam/sig)))
  
      ll_vec[i] <- ll
  
    }
    
    ll <- mean(ll_vec)
    
  } else {
    
    stop('\nMust indicate modeling structure:
         fixed, hetero, or random')
    
  }
  
  return(-ll)
  
}

sfa_ineff <- function(data, indicate) {
  
  y <- data[, 1] %>% as.matrix()
  X <- cbind(1, data[, -1]) %>% as.matrix()
  p <- dim(X)[2]
  
  beta_init <- init_fun(data)[['beta']]
  sig_init <- init_fun(data)[['sig']]
  lam_init <- 1
  theta <- c(beta_init, sig_init, lam_init)
  
  oprout <- optim(theta, sfa_ll, method = 'BFGS',
                  data = data, indicate = indicate)
  mles <- oprout$par
  conv <- oprout$convergence
  cat('Convergence exit status:', conv,
      '\n')
  
  sig_mle <- mles[p + 1]
  lam_mle <- mles[p + 2]
  e <- y - X%*%mles[1:p]
  
  if (indicate == 'fixed' | indicate == 'random') {
    
    z <- e*lam_mle/sig_mle
    
  }
  
  if (indicate == 'hetero') {
    
    h_idx <- grep('Inc|vil', names(data))
    h_mus <- X[, h_idx]%*%theta[h_idx]
    z <- e*lam_mle/sig_mle - h_mus/(sig_mle*lam_mle)
    
  }
  
  u_e <- (sig_mle*lam_mle/(1 + lam_mle^2))*((dnorm(z)/(1 - pnorm(z))) - z)
  p_u_e <- (max(u_e) - u_e)/max(u_e)
  
  out <- list(ineff = u_e, prop_ineff = p_u_e,
              mle = mles, error = e)
  return(out)
  
}


### Inefficiency results -------------------------------------------------------------

set.seed(7)
sfa_models <- list()
models <- c('random', 'hetero', 'fixed')
for (mod in models) {
  
  cat('Model iteration:', mod, '\n')
  sfa_models[[mod]] <- sfa_ineff(liv_prodx, mod)
  
}

# parameters
ols_beta <- init_fun(liv_prodx)$beta %>% round(6)
ols_sig <- init_fun(liv_prodx)$sig %>% round(6)

params_data <- data.frame(ols = c(ols_beta, ols_sig, '-'),
                          rand_mle = round(sfa_models[['random']]$mle, 6),
                          hetero_mle = round(sfa_models[['hetero']]$mle, 6))
  

row.names(params_data) <- c('const', 'total_acres', 'labor', 'off_farm_income',
                            'livestock_income', 'livestock_expense',
                            'livestock_treat_expense', 'vil68', 'vil55',
                            'vil10', 'vil28', 'vil53', 'vil13', 'vil35',
                            'vil49', 'vil67', 'sigma', 'lambda')

# ineff summary
u_data <- matrix(0, nrow = 4, ncol = length(models))
p_data <- matrix(0, nrow = 4, ncol = length(models))
for (j in 1:length(models)) {
  
  i_x <- sfa_models[[models[j]]]$ineff
  p_x <- sfa_models[[models[j]]]$prop_ineff
  u_data[, j] <- c(mean(i_x), sd(i_x), min(i_x), max(i_x))
  p_data[, j] <- c(mean(p_x), sd(p_x), min(p_x), max(p_x))
  
  if (j == length(models)) {
    
   ineff_data <- data.frame(cbind(u_data, p_data)) %>% round(6) 
   colnames(ineff_data) <- c('rand_dev', 'hetero_dev', 'fixed_dev',
                             'rand_ineff', 'hetero_ineff', 'fixed_ineff')
   row.names(ineff_data) <- c('mean', 'sd', 'min', 'max')
    
  }
  
}

# plotting
# e_plot <- ggplot(ineff_data) +
#   geom_density(aes(e), fill = 'gray', colour = 'gray', alpha = .35) +
#   labs(title = 'Full error component',
#        x = TeX('$e_i = v_i - u_i$'),
#        y = 'Density') +
#   theme(plot.title = element_text(family = 'serif', hjust = 0.5),
#         axis.title = element_text(family = 'serif'))
# 
# u_e_plot <- ggplot(ineff_data) +
#   geom_density(aes(u_e), fill = 'black', colour = 'gray', alpha = .5, size = 0.1) +
#   labs(title = TeX('Inefficiency component'),
#        x = TeX('$u_i|e_i$'),
#        y = '') +
#   theme(plot.title = element_text(family = 'serif', hjust = 0.5),
#         axis.title = element_text(family = 'serif'))
#   
# liv_plots <- plot_grid(e_plot, u_e_plot)

# write ineff and loss
liv_data_out <- cbind(sfa_models[['random']]$ineff, sfa_models[['random']]$prop_ineff,
                      liv_dates_hhID, liv_prodx, liv_loss, liv_comgraz, liv_edu)
names(liv_data_out)[1:2] <- c('ineff', 'p_ineff')
# write_csv(liv_data_out, 'data/liv_ineff_loss5.csv')


### Testing ---------------------------------------------------------------
