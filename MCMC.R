# 2021-01-05
# Model validation
#

library(tidyverse)
library(lubridate)
library(ggpubr)
library(cowplot)
library(parallel)
# settings --------
#et= "increasing"
#et= "decreasing"
et = "uniform"

n_sim = 10 
inf_ref = 14
T_max = 14
TimeFrame = 202008:202112
time = 202112
n_used = 20000
p_i_used = 3/1000

iter <- 10000; burnin <- 3000

# functions --------
#

# iter <- 20000; burnin <- 3000

#
# function to clean false test negative data


#

linearfunctionP<-function(x,a,b, direction){
  if(direction == "decreasing"){
    p = (seq(a,b,length.out=x))/sum(seq(a,b,length.out=x))
  }
  
  if(direction == "increasing"){
    p = (b+1-seq(a,b,length.out=x))/sum(seq(a,b,length.out=x))
  }
  
  if(direction == "uniform"){
    p = rep(1/T_max, T_max)
  }
  
  return(p)
  
}

#
# function to clean false test negative data
#

getFTN = function(p_tn_vec, infectious_period){

ftn = tibble(p_tn_vec) %>%
  select(days_since_exposure, fnr_med) %>%
  expand(days_since_exposure = -50:35) %>%
  left_join(p_tn_vec) 

ftn = ftn[!duplicated(ftn), ] %>%
  mutate(fnr_med = ifelse(is.na(fnr_med),
                          1,
                          fnr_med)) %>%
  mutate(cum_prob_infectious = ifelse(is.na(cum_prob_infectious),
                                      0,
                                      cum_prob_infectious)) %>%
  mutate(days_since_exposure_pre_departure = days_since_exposure,
         days_since_exposure_test = days_since_exposure,
         days_since_exposure_test_2 = days_since_exposure,
         days_since_exposure_test_3 = days_since_exposure,
         days_since_exposure_test_4 = days_since_exposure,
         days_since_exposure_q3 = days_since_exposure,
         days_since_exposure_q7 = days_since_exposure,
         days_since_exposure_q14 = days_since_exposure,
         days_since_exposure_q21 = days_since_exposure)

ftn}

#
# function to simulate data; also used to simulate the travel control effectiveness in the ms

sampleInfection = function(n, p_i, n_inf, T_max, ftn, incubation_period, inf_ref, direction,time){

if(n_inf > 0 ){
  
  
  exp_prob = linearfunctionP(x=T_max, a = 1, b = T_max, direction = direction)
  t = sample(1:T_max,n_inf,replace=TRUE,prob = exp_prob)
  #t = sample(1:T_max, n_inf, replace = TRUE)
  
  # generate incubation period (days since exposure)
  peak_day = sample(incubation_period$days_since_exposure[1:20],
                    n_inf,
                    replace = TRUE,
                    prob = incubation_period$prob_symptom_onset[1:20])
  
  # generate symptomatic
  p_asymp = rnorm(1,
                  mean = 0.1954023,
                  sd = 0.05046811)
  while (p_asymp >1 || p_asymp <0) {
    p_asymp = rnorm(1,
                    mean = 0.1954023,
                    sd = 0.05046811)}
  asympt = rbinom(n = n_inf, size = 1, p = p_asymp)
  
  n_asympt = sum(asympt)
  n_sympt = n_inf - n_asympt
  
  # matrix of events times
  if (time %in% 202008:202104 | time == 202112){
    dat = tibble(
      t = t,
      peak_day = peak_day,
      days_since_exposure_pre_departure = t - 2,
      days_since_exposure_test = t,
      days_since_exposure_test_2 = t + 7,
      days_since_exposure_test_3 = t + 14,
      days_since_exposure_test_4 = t + 21,
      days_since_exposure_q3 = t + 3,
      days_since_exposure_q7 = t + 7,
      days_since_exposure_q14 = t + 14,
      days_since_exposure_q21 = t + 21,
      asympt = asympt
    )}
  
  if (time %in% 202105:202108){
    dat = tibble(
      t = t,
      peak_day = peak_day,
      days_since_exposure_pre_departure = t - 2,
      days_since_exposure_test = t,
      days_since_exposure_test_2 = t + 5,
      days_since_exposure_test_3 = t + 10,
      days_since_exposure_test_4 = t + 20,
      days_since_exposure_q3 = t + 3,
      days_since_exposure_q7 = t + 7,
      days_since_exposure_q14 = t + 14,
      days_since_exposure_q21 = t + 21,
      asympt = asympt
    )}
  
  if (time %in% 202109:202111){
    dat = tibble(
      t = t,
      peak_day = peak_day,
      days_since_exposure_pre_departure = t - 2,
      days_since_exposure_test = t,
      days_since_exposure_test_2 = t + 6,
      days_since_exposure_test_3 = t + 12,
      days_since_exposure_test_4 = t + 21,
      days_since_exposure_q3 = t + 3,
      days_since_exposure_q7 = t + 7,
      days_since_exposure_q14 = t + 14,
      days_since_exposure_q21 = t + 21,
      asympt = asympt
    )}
  # all infected individuals who intended to travel
  
  # pass symptom screen and board
  # 1:Pass 48 hour negative----------
  dat = dat %>% 
    mutate(no_sympt_when_travel = ifelse(asympt == 1 | 
                                           (asympt == 0 & peak_day >= t), 
                                         1, 
                                         0)) %>% # all asym and show sym after travel can board
    mutate(escape_sympt_screen = ifelse(no_sympt_when_travel == 0, 
                                        sample(0:1, n(), replace = TRUE, prob = c(0.7, 0.3)),
                                        0))  # show sym before travel have 30% to board; 1 = yes 0 = no
  
  # test sensitivity and infectiousness 48 hours before departure -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_pre_departure', 
                      'fnr_med')], 
              by.x = 'days_since_exposure_pre_departure',
              by.y = 'days_since_exposure_pre_departure',
              all.x = TRUE) 
  colnames(dat)[colnames(dat) %in% c('fnr_med')] = 'fnr_med_test_pre_departure'  
  dat = dat %>%
    mutate(
      excape_48hr_test = rbinom(1:n_inf, 1, fnr_med_test_pre_departure), 
      board = ifelse(excape_48hr_test == 1,no_sympt_when_travel + escape_sympt_screen,0)
    )
  
  # test sensitivity and infectiousness on arrival -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_test', 
                      'fnr_med', 
                      'cum_prob_infectious')], 
              by.x = 'days_since_exposure_test',
              by.y = 'days_since_exposure_test',
              all.x = TRUE) 
  colnames(dat)[colnames(dat) %in% c('fnr_med', 
                                     'cum_prob_infectious')] = c('fnr_med_test',
                                                                 'cum_prob_infectious_test')
  
  # test sensitivity on 2nd test -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_test_2', 
                      'fnr_med')], 
              by.x = 'days_since_exposure_test_2',
              by.y = 'days_since_exposure_test_2',
              all.x = TRUE) 
  colnames(dat)[colnames(dat) %in% c('fnr_med')] = 'fnr_med_test_2'
  
  # test sensitivity on 3rd test -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_test_3', 
                      'fnr_med')], 
              by.x = 'days_since_exposure_test_3',
              by.y = 'days_since_exposure_test_3',
              all.x = TRUE) 
  colnames(dat)[colnames(dat) %in% c('fnr_med')] = 'fnr_med_test_3'  
  
  # test sensitivity on 4th test -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_test_4', 
                      'fnr_med')], 
              by.x = 'days_since_exposure_test_4',
              by.y = 'days_since_exposure_test_4',
              all.x = TRUE) 
  colnames(dat)[colnames(dat) %in% c('fnr_med')] = 'fnr_med_test_4'  
  
  # infectiousness on day 3 -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_q3', 
                      'cum_prob_infectious')], 
              by.x = 'days_since_exposure_q3',
              by.y = 'days_since_exposure_q3',
              all.x = TRUE)
  colnames(dat)[colnames(dat) %in% c('cum_prob_infectious')] = 'cum_prob_infectious_q3'
  
  # infectiousness on day 7 -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_q7', 
                      'cum_prob_infectious')], 
              by.x = 'days_since_exposure_q7',
              by.y = 'days_since_exposure_q7',
              all.x = TRUE)
  colnames(dat)[colnames(dat) %in% c('cum_prob_infectious')] = 'cum_prob_infectious_q7'
  
  # infectiousness on day 14 -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_q14', 
                      'cum_prob_infectious')], 
              by.x = 'days_since_exposure_q14',
              by.y = 'days_since_exposure_q14',
              all.x = TRUE)
  colnames(dat)[colnames(dat) %in% c('cum_prob_infectious')] = 'cum_prob_infectious_q14'
  
  # infectiousness on day 21 -------%
  dat = merge(dat, 
              ftn[ ,c('days_since_exposure_q21', 
                      'cum_prob_infectious')], 
              by.x = 'days_since_exposure_q21',
              by.y = 'days_since_exposure_q21',
              all.x = TRUE)
  colnames(dat)[colnames(dat) %in% c('cum_prob_infectious')] = 'cum_prob_infectious_q21'
  
  
  #
  # generate test-negative events
  #
  
  dat = dat %>%
    
    # test on-arrival
    mutate(tn_on_arrival = rbinom(n = n_inf, 
                                  size = 1, 
                                  prob = fnr_med_test)) %>%
    # those not boarding can't be false-negative
    mutate(tn_on_arrival = ifelse(board == 0,
                                  0,
                                  tn_on_arrival)) %>%
    
    # test on-day 5
    mutate(tn_on_q_2 = rbinom(n = n_inf, 
                              size = 1, 
                              prob = fnr_med_test_2)) %>%
    # on-arrival test-positives can't be false-negative on day 5
    mutate(tn_on_q_2 = ifelse(tn_on_arrival == 0,
                              0,
                              tn_on_q_2)) %>%
    
    # test on-day 12
    mutate(tn_on_q_3 = rbinom(n = n_inf, 
                              size = 1, 
                              prob = fnr_med_test_3)) %>%
    # on-arrival test-positives can't be false-negative on day 11
    mutate(tn_on_q_3 = ifelse(tn_on_arrival == 0,
                              0,
                              tn_on_q_3)) %>%
    
    # test on-day 19
    mutate(tn_on_q_4 = rbinom(n = n_inf, 
                              size = 1, 
                              prob = fnr_med_test_4)) %>%
    # on-arrival or day 12 test-positives can't be false-negative on day 19
    mutate(tn_on_q_4 = ifelse(tn_on_arrival == 0 | 
                                tn_on_q_3 == 0,
                              0,
                              tn_on_q_4))
  
  
  
  #
  # generate infectiousness
  #
  
  dat = dat %>%
    
    # infectiousness on-arrival
    mutate(infectious_on_arrival = ifelse(days_since_exposure_test < inf_ref &
                                            tn_on_arrival == 1,
                                          1, 0)) %>%
    mutate(remained_infectious_on_arrival = ifelse(days_since_exposure_test < inf_ref &
                                                     tn_on_arrival == 1,
                                                   1 - cum_prob_infectious_test,
                                                   0)) %>%
    # infectiousness on day 3
    mutate(infectious_q3 = ifelse(days_since_exposure_q3 < inf_ref &
                                    tn_on_arrival == 1,
                                  1, 0)) %>%
    mutate(remained_infectious_q3 = ifelse(days_since_exposure_q3 < inf_ref &
                                             tn_on_arrival == 1,
                                           1 - cum_prob_infectious_q3,
                                           0)) %>%
    
    # infectiousness on day 3 and test on arrival
    mutate(infectious_test_q3 = ifelse(tn_on_arrival == 1,
                                       infectious_q3, 
                                       0)) %>%
    mutate(remained_infectious_test_q3 = ifelse(tn_on_arrival == 1,
                                                remained_infectious_q3, 
                                                0)) %>%
    
    # infectiousness on day 7
    mutate(infectious_q7 = ifelse(days_since_exposure_q7 < inf_ref &
                                    tn_on_arrival == 1,
                                  1, 0)) %>%
    mutate(remained_infectious_q7 = ifelse(days_since_exposure_q7 < inf_ref &
                                             tn_on_arrival == 1,
                                           1 - cum_prob_infectious_q7,
                                           0)) %>%
    
    # infectiousness on day 7 and test on day 5
    mutate(infectious_test_q7 = ifelse(tn_on_q_2 == 1,
                                       infectious_q7, 
                                       0)) %>%
    mutate(remained_infectious_test_q7 = ifelse(tn_on_q_2 == 1,
                                                remained_infectious_q7, 
                                                0)) %>%
    
    
    # infectiousness on day 14
    mutate(infectious_q14 = ifelse(days_since_exposure_q14 < inf_ref &
                                     tn_on_arrival == 1,
                                   1, 0)) %>%
    mutate(remained_infectious_q14 = ifelse(days_since_exposure_q14 < inf_ref &
                                              tn_on_arrival == 1,
                                            1 - cum_prob_infectious_q14,
                                            0)) %>%
    
    # infectiousness on day 14 and test on day 12
    mutate(infectious_test_q14 = ifelse(tn_on_q_3 == 1,
                                        infectious_q14, 
                                        0)) %>%
    mutate(remained_infectious_test_q14 = ifelse(tn_on_q_3 == 1,
                                                 remained_infectious_q14, 
                                                 0)) %>%
    
    
    # infectiousness on day 21
    mutate(infectious_q21 = ifelse(days_since_exposure_q21 < inf_ref &
                                     tn_on_arrival == 1,
                                   1, 0)) %>%
    mutate(remained_infectious_q21 = ifelse(days_since_exposure_q21 < inf_ref &
                                              tn_on_arrival == 1,
                                            1 - cum_prob_infectious_q21,
                                            0))  %>%
    
    # infectiousness on day 21 and test on day 19
    mutate(infectious_test_q21 = ifelse(tn_on_q_4 == 1,
                                        infectious_q21, 
                                        0)) %>%
    mutate(remained_infectious_test_q21 = ifelse(tn_on_q_4 == 1,
                                                 remained_infectious_q21, 
                                                 0))
  
  
  #
  # output
  #
  
  if (time %in% 202002:202012){
    rt = 
      tibble(
        n_travel = n,
        p_i = p_i,
        
        n_inf = n_inf,
        n_boarding = sum(dat$board),
        n_asympt = sum(dat$asympt),
        
        n_tn_on_arrival = sum(dat$tn_on_arrival),
        
        infectious_on_arrival = sum(dat$infectious_on_arrival),
        infectious_q3 = sum(dat$infectious_q3),
        infectious_test_q3 = sum(dat$infectious_test_q3),
        infectious_q7 = sum(dat$infectious_q7),
        infectious_test_q7 = sum(dat$infectious_test_q7),
        infectious_q14 = sum(dat$infectious_q14),
        infectious_test_q14 = sum(dat$infectious_test_q14),
        infectious_q21 = 0,
        infectious_test_q21 = 0,
        
        remained_infectious_on_arrival = sum(dat$remained_infectious_on_arrival),
        remained_infectious_q3 = sum(dat$remained_infectious_q3),
        remained_infectious_test_q3 = sum(dat$remained_infectious_test_q3),
        remained_infectious_q7 = sum(dat$remained_infectious_q7),
        remained_infectious_test_q7 = sum(dat$remained_infectious_test_q7),
        remained_infectious_q14 = sum(dat$remained_infectious_q14),
        remained_infectious_test_q14 = sum(dat$remained_infectious_test_q14),
        remained_infectious_q21 = 0,
        remained_infectious_test_q21 = 0
        
      )}
  if (time %in% 202101:202203){
    rt = tibble(
      n_travel = n,
      p_i = p_i,
      
      n_inf = n_inf,
      n_boarding = sum(dat$board),
      n_asympt = sum(dat$asympt),
      
      n_tn_on_arrival = sum(dat$tn_on_arrival),
      
      infectious_on_arrival = sum(dat$infectious_on_arrival),
      infectious_q3 = sum(dat$infectious_q3),
      infectious_test_q3 = sum(dat$infectious_test_q3),
      infectious_q7 = sum(dat$infectious_q7),
      infectious_test_q7 = sum(dat$infectious_test_q7),
      infectious_q14 = sum(dat$infectious_q14),
      infectious_test_q14 = sum(dat$infectious_test_q14),
      infectious_q21 = sum(dat$infectious_q21),
      infectious_test_q21 = sum(dat$infectious_test_q21),
      
      remained_infectious_on_arrival = sum(dat$remained_infectious_on_arrival),
      remained_infectious_q3 = sum(dat$remained_infectious_q3),
      remained_infectious_test_q3 = sum(dat$remained_infectious_test_q3),
      remained_infectious_q7 = sum(dat$remained_infectious_q7),
      remained_infectious_test_q7 = sum(dat$remained_infectious_test_q7),
      remained_infectious_q14 = sum(dat$remained_infectious_q14),
      remained_infectious_test_q14 = sum(dat$remained_infectious_test_q14),
      remained_infectious_q21 = sum(dat$remained_infectious_q21),
      remained_infectious_test_q21 = sum(dat$remained_infectious_test_q21)
      
    )
  }
  
} else {
  
  rt = tibble(
    n_travel = n,
    p_i = p_i,
    n_inf = n_inf,
    n_boarding = 0,
    n_asympt = 0,
    n_tn_on_arrival = 0,
    infectious_on_arrival = 0,
    infectious_q3 = 0,
    infectious_test_q3 = 0,
    infectious_q7 = 0,
    infectious_test_q7 = 0,
    infectious_q14 = 0,
    infectious_test_q14 = 0,
    infectious_q21 = 0,
    infectious_test_q21 = 0,
    remained_infectious_on_arrival = 0,
    remained_infectious_q3 = 0,
    remained_infectious_test_q3 = 0,
    remained_infectious_q7 = 0,
    remained_infectious_test_q7 = 0,
    remained_infectious_q14 = 0,
    remained_infectious_test_q14 = 0,
    remained_infectious_q21 = 0,
    remained_infectious_test_q21 = 0
  )
  
}


rt



}

#
# log prior for parameters
#
logprior <- function(current, data, p_tn_vec){
  
  #
  # function to calculate the log prior
  #
  
  # p1 = dnorm(current$n, 
  #            mean = data$Arrivals, 
  #            sd = 10^4, 
  #            log = TRUE)
  p1 = dunif(current$n, 
             min = 0, 
             max = data$Arrivals*2, 
             log = TRUE)
  p2 = dunif(current$p_i/10^4, 0, 1, log=TRUE)
  
  p_3 = dnorm(current$p_asymt,
              mean = 0.1954023,
              sd = 0.05046811,
              log = TRUE)
  
  p_4 = dnorm(current$p_fnr,
              mean = p_tn_vec$fnr_med,
              sd = p_tn_vec$fnr_sd,
              log = TRUE)
  
  current$logprior <- sum(p1) + sum(p2) + p_3 + sum(p_4)
  
  
  current
  
}


#
# log likelihood function
#
loglikelihood <- function(current, data, p_t, incubation_period, p_tn_vec){
  
  p_tp_vec = tibble(
    
    days_since_exposure = p_tn_vec$days_since_exposure,
    tpr_med = 1 - current$p_fnr
    
  )
  
  #
  # generate infections from travelers
  # 
  n_inf = current$n *current$p_i/10^4
  
  # course of disease (days since exposure)
  t = 1:T_max
  p_symptom_onset = incubation_period$prob_symptom_onset[1:T_max]
  p_tp = p_tn_vec$fnr_med
  
  # probability to board and test
  p_healthy = 1 - current$p_i/10^4
  p_asymt = current$p_asymt
  
  p_symt_board = (1 - p_asymt)*sum(unlist(lapply(t, function(i){
    
    # i - exposure day
    p_tem = ifelse(t >= i, 1, 0.3) # incubation period >= exposure; if yes 1, if no 0.3
    p_t[i]*sum(p_symptom_onset*p_tem)
    
  })))
  
  p_neg_48_board =sum(unlist(lapply(t, function(i){
    p_t[i]*(1-p_tp_vec$tpr_med[p_tp_vec$days_since_exposure == (i-2)])
    })))
  
  p_test = p_healthy + (p_asymt + p_symt_board)*p_neg_48_board *current$p_i/10^4
  
  # number of test positives
  p_tp_asymt = p_asymt*sum(unlist(lapply(t, function(i){
    
    p_t[i]*sum(p_symptom_onset*p_tp_vec$tpr_med[p_tp_vec$days_since_exposure %in% t])
    # p_t*sum(p_symptom_onset * as.numeric(t >= i))
    
  })))
  
  p_tp_symt = (1 - p_asymt)*sum(unlist(lapply(t, function(i){
    
    # i - exposure day
    p_tem = ifelse(t >= i, 1, 0.3) # incubation period >= exposure; if yes 1, if no 0.3
    p_t[i]*sum(p_symptom_onset*p_tem*p_tp_vec$tpr_med[p_tp_vec$days_since_exposure %in% t])
    
  })))
  
  
  p_tp = (p_tp_asymt + p_tp_symt)*p_neg_48_board
  
  
  
  #
  # function to calculate likelihood in log scale
  #
  
  N_pos = current$n * current$p_i/10^4 * p_tp
  N_case = current$n * current$p_i/10^4 *(p_asymt + p_symt_board)*p_neg_48_board
  
  ll <- dbinom(data$Arrivals, round(current$n), p_test, log=TRUE) +
    dpois(data$Positives, N_pos, log=TRUE)
  
  
  current$loglikelihood <- sum(ll)
  
  current
  
}

bad_config <- function(current){
  
  #
  # function for parameter constrain
  #
  
  ok_config <- TRUE
  if(!all(current$n > 0) |
     !all(current$p_i >= 0) | 
     !all(current$p_i < 10^4) |
     current$p_asymt < 0 |
     current$p_asymt > 1 |
     !all(current$p_fnr >= 0) | 
     !all(current$p_fnr <= 1)) ok_config <- FALSE
  bad <- FALSE; if(!ok_config) bad <- TRUE
  bad
  
}


#
# Metropolis-Hastings method for updating the parameter values
#

metropolis <- function(old, current, data, p_t, incubation_period, p_tn_vec){
  
  REJECT <- bad_config(current)
  
  if(!REJECT){
    current <- logprior(current, data, p_tn_vec)
    current <- loglikelihood(current, data, p_t, incubation_period, p_tn_vec)
    lu <- log(runif(1))
    if(lu > (current$logprior + sum(current$loglikelihood) - old$logprior - sum(old$loglikelihood))) REJECT <- TRUE
    
  }
  
  if(REJECT){current <- old}
  
  current
  
}


#
# main function for MCMC process
#

mcmc <- function(data,
                 p_t, incubation_period, p_tn_vec, 
                 MCMC_interation=1000, BURNIN_interation=100){
  
  par.length = nrow(data)*2
  
  current <- list(n = data$Arrivals+20,
                  p_i = rep(0.5, par.length/2),# death/severe
                  p_asymt = 0.21,
                  p_fnr = p_tn_vec$fnr_med,
                  logprior=0, loglikelihood=0,
                  infection=0, boosting=0)
  
  
  interation.length = MCMC_interation+BURNIN_interation
  
  dump <- list(n = matrix(NA, 
                          ncol = par.length/2,
                          nrow = interation.length),
               p_i = matrix(NA, 
                            ncol = par.length/2,
                            nrow = interation.length),
               p_asymt = rep(NA, interation.length),
               p_fnr = matrix(NA, 
                              ncol = nrow(p_tn_vec),
                              nrow = interation.length),
               accept = matrix(NA,
                               ncol = par.length + 1 + nrow(p_tn_vec), 
                               nrow = interation.length),
               likelihood = c(),
               infection = c()) 
  
  
  sigma <- list(# n = rep(10, par.length/2),
    n = floor(data$Arrivals/1000)+2,
     p = rep(1.5, par.length/2),
    p_asymt = 0.05,
    p_fnr = rep(0.02, nrow(p_tn_vec)))
  
  current <- logprior(current, data, p_tn_vec)
  current <- loglikelihood(current, data, p_t, incubation_period, p_tn_vec)
  
  for(interation in (-BURNIN_interation+1):MCMC_interation){
    
    for(i_para in 1:2){
      
      for(ii_para in 1:(par.length/2)){
        
        old <- current
        current[[i_para]][ii_para] <- rnorm(1, 
                                            current[[i_para]][ii_para], 
                                            sigma[[i_para]][ii_para])
        current <- metropolis(old, current, data, p_t, incubation_period, p_tn_vec)
        
      }
      
      
    }
    
    i_para = 3
    old <- current
    current[[i_para]] <- rnorm(1, 
                               current[[i_para]], 
                               sigma[[i_para]])
    
    current <- metropolis(old, current, data, p_t, incubation_period, p_tn_vec)
    
    for(i_para in 4){
      
      for(ii_para in 1:nrow(p_tn_vec)){
        
        old <- current
        current[[i_para]][ii_para] <- rnorm(1, 
                                            current[[i_para]][ii_para], 
                                            sigma[[i_para]][ii_para])
        current <- metropolis(old, current, data, p_t, incubation_period, p_tn_vec)
        
      }
      
      
    }
    
    
    rr <- interation + BURNIN_interation
    
    for(i_para in 1:length(sigma)){
      
      if(i_para != 3){
        
        dump[[i_para]][rr, ] <- current[[i_para]]
        
        if(i_para < 3){
          col.used = ((i_para - 1)*par.length/2 + 1):(i_para * par.length/ 2)
        } else {
          col.used = (par.length + 2):ncol(dump$accept)
        }
        
        if(rr > 1){
          dump$accept[rr, col.used] <- 1*(dump[[i_para]][rr, ] != dump[[i_para]][rr-1, ])
        }
        
      }
      
      
      if(i_para == 3) {
        
        dump[[i_para]][rr] <- current[[i_para]]
        
        if(rr > 1){
          dump$accept[rr, par.length + i_para - 2] <- 1*(dump[[i_para]][rr] != dump[[i_para]][rr-1])
        }
        
      }
      
    }
    
    
    dump$likelihood[rr] <- current$loglikelihood
  }
  
  dump
  
}

#
# summary results from mcmc
#

smrResults = function(data, estimates){
  
  # intend-to-travel
  # n = estimates$n[1001:3000, ]
n = lapply(estimates, function(x) x$n[(burnin+1):(iter+burnin)])
  n = do.call('cbind', n)
  
  smr = lapply(1:ncol(n), function(x){
    
    quantile(n[,x],
             c(.025, .5, .975))
    
  })
  smr = do.call('rbind', smr)
  
  colnames(smr) = c('n_025', 'n_50', 'n_975')
  
  
  # prevalence
  # prev = estimates$p_i[1001:3000, ]
  prev = lapply(estimates, function(x) x$p_i[(burnin+1):(iter+burnin)])
  prev = do.call('cbind', prev)
  
  smr2 = lapply(1:ncol(prev), function(x){
    
    quantile(prev[,x],
             c(.025, .5, .975))/10^4 * 1000
    
  })
  smr2 = do.call('rbind', smr2)
  
  colnames(smr2) = c('p_i_025', 'p_i_50', 'p_i_975')
  
  smr = cbind(data, smr, smr2)
  
  smr
  
}

#
# data --------
#

# extracted from https://github.com/cmmid/pcr-profile
p_tn_vec = read.csv("test_sensitivity_uk.csv",
                    header = TRUE) %>%
  mutate(days_since_exposure = days_since_peak_viral_load + 5) %>%
  select(days_since_exposure, fnr_med, fnr_lb,fnr_ub,fnr_sd)

# estimated from https://github.com/ehylau/COVID-19
incubation_period = read.csv('incubation_period.csv',
                             header = TRUE,
                             stringsAsFactors = FALSE) %>%
  mutate(prob_symptom_onset = ifelse(days_since_exposure <= 14,
                                     prob_symptom_onset/sum(prob_symptom_onset[1:14]),
                                     0)) %>%
  select(days_since_exposure, prob_symptom_onset)

# estimated from https://github.com/ehylau/COVID-19
infectious_period = read.csv("nature_med_infectiousness.csv",
                             header = TRUE,
                             stringsAsFactors = FALSE) %>%
  select(days_since_exposure,prob_infectious,cum_prob_infectious)


p_tn_vec = left_join(p_tn_vec, infectious_period)

ftn = getFTN(p_tn_vec, infectious_period)


data_model_validation = lapply(n_used, function(n){
  
  sim = lapply(p_i_used, function(p_i){
    
    sim = lapply(1:n_sim, function(x) {
      
      n_inf = rbinom(1, n, p_i)
      sim = sampleInfection(n, 
                            p_i,
                            T_max = 14,
                            n_inf = n_inf, 
                            ftn,
                            time = time,
                            incubation_period,
                            inf_ref = inf_ref,
                            direction = et) %>%
        mutate(
          Arrivals = n_travel - n_inf + n_boarding,
          Positives = n_boarding - n_tn_on_arrival
        ) %>%
        select(n_travel, 
               p_i, 
               n_inf, 
               n_boarding,
               Arrivals,
               Positives)
      
      sim
      
    })
    sim = do.call('rbind', sim)
    
    sim
    
  })
  
  sim = do.call('rbind', sim)
  
  sim
  
})

data = do.call("rbind", data_model_validation)





#bookmark3----

p_t<-linearfunctionP(T_max,1,T_max,et)
#p_t<-1/T_max
#logunif<-function(x,a,b){1/(x*log(b/a))}
#p_t<-logunif(1:T_max,1,T_max)

results<- lapply(1:nrow(data), function(kk) {
  
  print(Sys.time())
  #print(kk)
  results <- mcmc(data[kk, ], 
                  p_t, incubation_period, p_tn_vec, iter, burnin)
  print(Sys.time())
  results
  
})  

stopCluster(cl)



#
# summary results -------%
#



smr<-smrResults(data[1:nrow(data), ], results_new)
