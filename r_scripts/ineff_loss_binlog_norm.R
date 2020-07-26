# Inefficiency - Livestock health model
# Bayesian bin-log normal approximation

remove(list = objects())
setwd('~/research/africa/Health_Ineff/')
library(tidyverse)
library(dotwhisker)
library(latex2exp)


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

comgraz <- ifelse(dt$UnsharedLand == 0, 1, 0)

vil_names <- names(dt)[grep('vil', names(dt))]
# excludes village 67 as reference village
vils <- dt[, vil_names[1:length(vil_names) - 1]]
# village interaction under livestock general illness
vil_illness <- sapply(vils,
                      function(i) i*dt$GeneralIllness)

# features:
## household education level, 1 primary or sec edu, 0 otherwise
## livestock general illness, 1 illness, 0 otherwise
## respiratory/skin disorder, 1 if, 0 otherwise
## animal treated, 1 if, 0 otherwise
X <- cbind(1,
           # dt$no_edu,
           dt$GeneralIllness,
           dt$no_edu*dt$GeneralIllness,
           respskin,
           dt$no_edu*respskin,
           comgraz*dt$GeneralIllness,
           comgraz*respskin,
           dt$AnimalTreated,
           vils) %>% 
  as.matrix(nrow = nrow(dt))

param_names <- c('intercept',
                 # 'no_edu',
                 'gen_illness', 'gen_illness_no_edu',
                 'resp_skin_disorder', 'resp_skin_disorder_no_edu',
                 'comgraz_gen_illness', 'comgraz_resp_skin_disorder',
                 'treated',
                 names(vils))

tcrit <- qt(1-.05/2, dim(X)[1] - dim(X)[2])


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

# beta convergence plot

# ggplot(beta_vals, aes(1:nrow(beta_vals))) +
#   geom_line(aes(y = b0, colour = 'intercept')) +
#   geom_line(aes(y = b1, colour = 'gen_illness')) +
#   geom_line(aes(y = b2, colour = 'gen_illness_no_edu')) +
#   geom_line(aes(y = b3, colour = 'resp_skin_disorder')) +
#   geom_line(aes(y = b4, colour = 'resp_skin_disorder_no_edu')) +
#   geom_line(aes(y = b5, colour = 'comgraz_gen_illness')) +
#   geom_line(aes(y = b6, colour = 'comgraz_resp_skin_disorder')) +
#   geom_line(aes(y = b7, colour = 'treated')) +
#   geom_line(aes(y = b8, colour = 'vil68')) +
#   geom_line(aes(y = b9, colour = 'vil55')) +
#   geom_line(aes(y = b10, colour = 'vil10')) +
#   geom_line(aes(y = b11, colour = 'vil28')) +
#   geom_line(aes(y = b12, colour = 'vil53')) +
#   geom_line(aes(y = b13, colour = 'vil13')) +
#   geom_line(aes(y = b14, colour = 'vil35')) +
#   geom_line(aes(y = b15, colour = 'vil49')) +
#   # scale_color_manual(name = '',
#   #                    values = c('b_0' = 'black', 'b_1' = 'orange',
#   #                               'b_2' = 'blue', 'b_3' = 'red')) +
#   labs(title = 'Figure 1. Estimator Convergence (diffuse prior)',
#        x = 'Iteration',
#        y = TeX('$\\beta$')) +
#   theme(plot.title = element_text(family = 'serif', hjust = 0.5),
#         axis.title = element_text(family = 'serif'),
#         legend.text = element_text(family = 'serif')) +
#   theme_bw()


# Posteriors --------------------------------------------------------------

post_list <- vector(mode = 'list', length = length(beta_mp1))

# Log_odds parameters
## normal posterior - MAP is mean or median
beta_post <- lapply(1:length(post_list),
                    function(i)
                      post_list[[i]] = rnorm(n = 1000,
                                             mean = beta_mp1[i],
                                             sd = sqrt(beta_mp1var[i])))
beta_map <- sapply(beta_post,
                   function(i) mean(i))
beta_se <- sapply(beta_post,
                  function(i) sqrt(var(i)))
lb <- beta_map - tcrit*beta_se
ub <- beta_map + tcrit*beta_se

param_table <- data.frame(log_odds_coef = beta_map,
                          lb = lb,
                          ub = ub) %>% 
  round(4)
row.names(param_table) <- param_names

# Odds parameters
beta_odds_post <- lapply(beta_post,
                         function(i) exp(i))
beta_odds_map <- sapply(beta_odds_post,
                        function(i) mean(i))
beta_odds_se <- sapply(beta_odds_post,
                       function(i) sqrt(var(i)))

lb_odds <- beta_odds_map - tcrit*beta_odds_se
ub_odds <- beta_odds_map + tcrit*beta_odds_se

param_odds_table <- data.frame(odds_coef = beta_odds_map,
                               lb = lb_odds,
                               ub = ub_odds) %>% 
  round(4)
row.names(param_odds_table) <- param_names


# Parameter plot ----------------------------------------------------------

# Log odds
plt_df <- data.frame(term = param_names,
                     estimate = beta_map,
                     std.error = beta_se)

# plt_df %>% 
#   dwplot() +
#   theme_bw() +
#   theme(legend.position = 'none') +
#   ggtitle('Estimated Log-Odds (Coefficients) Impacts on Being \nBelow Mean Household Livestock Production Inefficiency') +
#   theme(plot.title = element_text(family = 'serif', hjust = 0.5),
#         axis.title = element_text(family = 'serif'),
#         axis.text = element_text(family = 'serif')) +
#   labs(x = 'Log-odds Coefficient Estimate') +
#   geom_vline(xintercept = 0, colour = 'grey60', linetype = 2)
# 
# # Odds
plt_odds_df <- data.frame(term = param_names,
                          estimate = beta_odds_map,
                          std.error = beta_odds_se)

# plt_odds_df %>% 
#   dwplot() +
#   theme_bw() +
#   theme(legend.position = 'none') +
#   ggtitle('Odds Impacts on Being Below Mean Household Livestock Production Inefficiency') +
#   theme(plot.title = element_text(family = 'serif', hjust = 0.5),
#         axis.title = element_text(family = 'serif'),
#         axis.text = element_text(family = 'serif')) +
#   labs(x = 'Log-odds Coefficient Estimate') +
#   geom_vline(xintercept = 1, colour = 'grey60', linetype = 2)


# Counterfactual analysis -------------------------------------------------

get_prob <- function(vals) {
  
  B <- matrix(beta_map)
  x <- t(matrix(vals))
  
  p <- 1/(1 + exp(-x%*%B))
  return(p)
  
}

get_prob_vec <- function(j_vals) {
  
  vec <- sapply(1:length(vils_vallist),
                function(i)
                  get_prob(vils_vallist[[i]][[j_vals]]))
  return(vec)
  
}

int_vec <- c(1, rep(0, length(param_names) - 1))
vils_vallist <- vector(mode = 'list', length = length(names(vils)))
names(vils_vallist) <- names(vils)
vils_vallist <- lapply(vils_vallist,
                       function(i)
                         i = list(bench_vals = int_vec,
                                  gen_illvals = int_vec,
                                  gen_ill_noeduvals = int_vec,
                                  rs_disvals = int_vec,
                                  rs_dis_noeduvals = int_vec,
                                  comgraz_illval = int_vec,
                                  comgraz_rsval = int_vec,
                                  treated_genillval = int_vec,
                                  treated_rsval = int_vec))

for (i in 1:length(vils_vallist)) {
  
  vil_idx <- grep('vil', param_names)[i]
  
  vils_vallist[[i]]$bench_vals[grep('vil', param_names)[i]] <- 1
  vils_vallist[[i]]$gen_illvals[c(grep('^gen_illness$', param_names),
                                  vil_idx)] <- 1
  vils_vallist[[i]]$gen_ill_noeduvals[c(grep('^gen_illness', param_names),
                                        vil_idx)] <- 1
  vils_vallist[[i]]$rs_disvals[c(grep('^resp_skin_disorder$', param_names),
                                 vil_idx)] <- 1
  vils_vallist[[i]]$rs_dis_noeduvals[c(grep('^resp_skin', param_names),
                                       vil_idx)] <- 1
  vils_vallist[[i]]$comgraz_illval[c(grep('^comgraz_gen', param_names),
                                     vil_idx)] <- 1
  vils_vallist[[i]]$comgraz_rsval[c(grep('^comgraz_resp', param_names),
                                    vil_idx)] <- 1
  vils_vallist[[i]]$treated_genillval[c(grep('treat|^gen_illness$', param_names),
                                        vil_idx)] <- 1
  vils_vallist[[i]]$treated_rsval[c(grep('treat|^resp_skin_disorder$', param_names),
                                    vil_idx)] <- 1
  
}

# events P(being below inefficiency mean)
prob_events <- c('bench_prob', 'gen_ill_prob', 'gen_ill_noedu_prob',
                 'respskin_prob', 'respskin_noedu_prob',
                 'comgraz_genill_prob', 'comgraz_respskin_prob',
                 'treat_genill_prob', 'treat_respskin_prob')

prob_list <- lapply(1:length(prob_events),
                    function(i)
                      get_prob_vec(i))
names(prob_list) <- prob_events

prob_events_df <- bind_rows(prob_list) %>% 
  data.frame() %>% 
  round(4)
row.names(prob_events_df) <- names(vils)

# probability event differences
prob_diff_df <- prob_events_df %>%
  transmute(gen_ill_diff = gen_ill_prob - bench_prob,
            gen_ill_noedu_diff = gen_ill_noedu_prob - bench_prob,
            respskin_diff = respskin_prob - bench_prob,
            respskin_noedu_diff = respskin_noedu_prob - bench_prob,
            comgraz_genill_diff = comgraz_genill_prob - bench_prob,
            comgraz_respskin_diff = comgraz_respskin_prob - bench_prob,
            treat_genill_diff = treat_genill_prob - bench_prob,
            treat_respskin_diff =  treat_respskin_prob - bench_prob)
row.names(prob_diff_df) <- names(vils)

# TODO --------------------------------------------------------------------

