# Inefficiency - Livestock health model
# Bayesian bin-log normal approximation

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)


# Data mgmt ---------------------------------------------------------------

dt <- read_csv('data/ineff_loss_sym4.csv')
dt$PoorCondition[which(dt$PoorCondition == 99)] <- 0

# response label:
## indicates whether household inefficiency is
## greater than the mean (0) or less than the mean (1)
y <- ifelse(dt$p_ineff >= mean(dt$p_ineff), 0, 1)

edu <- ifelse(dt$prim_edu == 1 |
                dt$sec_edu == 1, 1, 0)

respskin <- ifelse(dt$RespDisorders == 1 |
                     dt$SkinDisorders == 1, 1, 0)

# features:
## household education level, 1 primary or sec edu, 0 otherwise
## livestock general illness, 1 illness, 0 otherwise
## respiratory/skin disorder, 1 if, 0 otherwise
X <- matrix(cbind(1,
                  edu,
                  dt$GeneralIllness,
                  #dt$PoorCondition,
                  respskin),
            nrow = nrow(dt))



# Model funcs -------------------------------------------------------------

# pseudodata for normal z
# parameters: eta as linear combination X'B for logistic link
#             sig2 for y -> z ~ N(eta, sig2)
# beta linear weight estimates

z_gen <- function(X, B, n = length(y)) {
  
  B <- matrix(B)
  eta <- X%*%B
  z <- eta + (1 + exp(eta)^2)/exp(eta)*(y/n - exp(eta)/(1 + exp(eta)))
  return(z)
  
}


sig2_gen <- function(X, B, n = length(y)) {
  
  B <- matrix(B)
  eta <- X%*%B
  s <- 1/n*(1 + exp(eta))^2/exp(eta)
  return(s)
  
}

beta_approx <- function(X, W, y) {
  
  y <- matrix(y)
  b <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y
  return(b)
  
}

beta_var <- function(X, B, n = length(y)) {
  
  B <- matrix(B)
  eta <- X%*%B
  v <- -(-n*exp(eta)/(1 + exp(eta))^2)
  vmat <- solve(t(X)%*%diag(c(v))%*%X)
  return(vmat)
  
}


# Optimization ------------------------------------------------------------

# iterative normal approximation

# set convergence error and constrain iterations
err <- 1
tol <- 1E-10
itr <- 1
max_itr <- 100

# inital beta vals and storage
beta_m <- rep(0, dim(X)[2])
beta_names <- paste('b', 0:(dim(X)[2] - 1),
                    sep = '')
beta_vals <- data.frame()


# routine
while(err > tol & itr <= max_itr) {
  
  z <- z_gen(X, beta_m)
  sig2 <- sig2_gen(X, beta_m)
  W <- diag(1/as.vector(sig2))
  beta_mp1 <- beta_approx(X, W, y)
  beta_vals[itr, beta_names] <- c(beta_mp1)
  
  err <- max(abs(beta_mp1 - beta_m))
  
  cat('Iteration', itr, 'max convergence error:', err,
      '\n')
  
  if (err < tol) {
    
    beta_mp1var <- diag(beta_var(X, beta_mp1))
    
  }
  
  beta_m <- beta_mp1
  itr <- itr + 1
  
}



# Posteriors --------------------------------------------------------------

post_list <- vector(mode = 'list', length = length(beta_mp1))
beta_post <- lapply(1:length(post_list),
                    function(i)
                      post_list[[i]] = rnorm(n = 1000,
                                             mean = beta_mp1[i],
                                             sd = sqrt(beta_mp1var[i])))
# Parameter CIs
tcrit <- qt(1-.05/2, dim(X)[1] - dim(X)[2])
lb <- sapply(beta_post,
             function(i) mean(i) - tcrit*sqrt(var(i)))
ub <- sapply(beta_post,
             function(i) mean(i) + tcrit*sqrt(var(i)))
param_table <- data.frame(odds_coef = beta_mp1,
                          lb = lb,
                          ub = ub) %>% 
  exp() %>% 
  round(4)
param_names <- c('int', 'edu', 'gen_illness', 'resp_skin_disorder')
row.names(param_table) <- param_names

 

# TODO --------------------------------------------------------------------

# graph beta convergence (extend iterations)
# generate beta posteriors
# beta inference

















