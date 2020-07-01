# Inefficiency - Livestock health model
# Bayesian bin-log normal approximation

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)


# Data mgmt ---------------------------------------------------------------

dt <- read_csv('data/ineff_loss_sym4.csv')
dt$PoorCondition[which(dt$PoorCondition == 99)] <- 0

y <- ifelse(dt$p_ineff >= mean(dt$p_ineff), 0, 1)

respskin <- ifelse(dt$RespDisorders == 1 |
                     dt$SkinDisorders == 1, 1, 0)
X <- matrix(cbind(1,
                  dt$prim_edu,
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
  
  beta_m <- beta_mp1
  itr <- itr + 1
  
}



# TODO --------------------------------------------------------------------

# graph beta convergence (extend iterations)
# generate beta posteriors
# beta inference

















