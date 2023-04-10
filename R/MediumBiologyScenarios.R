#'
#' MediumBiologyScenario
#'
#'

source("AuxillaryFunctions.R")
library(dplyr)
library(ggplot2)
library(reshape2)
library(TMB)
library(tidyr)
## Pass the OM generated data to the TMB model
#sink(file = "compile_output.txt")
compile(file = file.path(DIR$tmb, "AgeStructuredModel.cpp"), flags = "-Wignored-attributes -O3")
#sink()
#dyn.unload(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))
dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))
#setwd(DIR$R)

fig_dir = file.path(DIR$fig, "MediumBiology")
if(!dir.exists(fig_dir))
  dir.create(fig_dir)

flat_bio = readRDS(file = file.path(DIR$data, "Medium_biology.RDS"))

this_bio = flat_bio

n_years = 60
years = (2020 - n_years + 1):2020
full_years = (min(years) - 1):max(years)
## observation temporal frequency
survey_year_obs = 1980:2020
survey_ages = this_bio$ages
fishery_year_obs = 1980:2020
fishery_ages = this_bio$ages

############
## Build a multinomial model to double check estimability of all parameters
## In this case we have 'good' data, annual data, no ageing error.
##  survey index cv = 0.05
##  year effective sample size = 1000
############
TMB_data = list()
TMB_data$ages = this_bio$ages
TMB_data$maxAgePlusGroup = 1
TMB_data$years = years
TMB_data$n_years = length(TMB_data$years)
TMB_data$n_ages = length(TMB_data$ages)
TMB_data$n_fisheries = 1
## No ageing error
TMB_data$ageing_error_matrix = matrix(0, nrow = TMB_data$n_ages, ncol = TMB_data$n_ages)
diag(TMB_data$ageing_error_matrix) = 1;

TMB_data$survey_year_indicator = as.integer(TMB_data$years %in% survey_year_obs)
TMB_data$survey_obs = rep(0, sum(TMB_data$survey_year_indicator))
TMB_data$survey_cv = rep(0.15, sum(TMB_data$survey_year_indicator))
TMB_data$survey_AF_obs = matrix(5, nrow = TMB_data$n_ages, ncol = sum(TMB_data$survey_year_indicator))

TMB_data$fishery_year_indicator = array(as.integer(TMB_data$years %in% fishery_year_obs), dim = c(TMB_data$n_years, TMB_data$n_fisheries))
TMB_data$fishery_AF_obs = array(5, dim = c(TMB_data$n_ages, length(fishery_year_obs), TMB_data$n_fisheries))

TMB_data$catches = array(1000, dim = c(TMB_data$n_years, TMB_data$n_fisheries))# this will be overriden in the simulate() call
TMB_data$F_method = 0
TMB_data$F_iterations = 4
TMB_data$F_max = 3

TMB_data$catch_indicator = array(1, dim = c(TMB_data$n_years, TMB_data$n_fisheries))
TMB_data$ycs_estimated = c(rep(1, n_years))
TMB_data$standardise_ycs = 0;

TMB_data$catchMeanLength = TMB_data$stockMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = TMB_data$n_years, nrow = TMB_data$n_ages)
TMB_data$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = TMB_data$n_years, nrow = TMB_data$n_ages)
TMB_data$natMor = this_bio$M
TMB_data$steepness = this_bio$h
TMB_data$stockRecruitmentModelCode = 2 ## BH
TMB_data$propZ_ssb = rep(0.5, TMB_data$n_years)
TMB_data$propZ_survey = rep(0.5, TMB_data$n_years)
TMB_data$sel_ato95_bounds = c(0.1,20)
TMB_data$sel_a50_bounds = c(0.1,60)
TMB_data$mean_weight_a = this_bio$a
TMB_data$mean_weight_b = this_bio$b
TMB_data$estimate_F_init = 0
TMB_data$estimate_init_age_devs = 0
TMB_data$n_init_age_devs = 1

## iniital fishery_probs
fishery_probs = c(rep(this_bio$M * 1.5, 20), rep(this_bio$M, 40) )
plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)

## The same parameters as OM, to check for consistency
OM_pars = list(
  ln_R0 = log(this_bio$R0),
  ln_ycs_est =  rnorm(sum(TMB_data$ycs_estimated),  -0.5 * this_bio$sigma_r * this_bio$sigma_r, this_bio$sigma_r),
  ln_sigma_r = log( this_bio$sigma_r),
  ln_extra_survey_cv = log(0.0001),
  ln_F_init = log(0.01),
  ln_init_age_devs = rep(log(0.1)),
  ln_sigma_init_age_devs = log(0.6),
  logit_f_a50 = logit_general(rep(this_bio$f_a50, TMB_data$n_fisheries), TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  logit_f_ato95 = logit_general(rep(this_bio$f_ato95, TMB_data$n_fisheries), TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
  logit_survey_a50 = logit_general(this_bio$s_a50, TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  logit_survey_ato95 = logit_general(this_bio$s_ato95, TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
  logit_surveyQ = qlogis(0.2),
  ln_F = array(log(fishery_probs), dim = c(TMB_data$n_fisheries,TMB_data$n_years)),
  ln_catch_sd = log(0.02),
  ln_Fmax = log(0.05),
  ln_F40 = log(0.05),
  ln_F35 = log(0.05),
  ln_F30 = log(0.05),
  ln_Fmsy = log(0.06),
  ln_F_0_1 = log(0.04)
  
)

# these parameters we are not estimating.
na_map = fix_pars(par_list = OM_pars, pars_to_exclude = c("ln_catch_sd", "ln_extra_survey_cv","ln_sigma_r", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs"))
OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)

OM_report = OM_obj$report()
OM_report$F30_nll
OM_report$F35_nll
OM_report$F40_nll
OM_report$F_0_1_nll

OM_report$B_msy
OM_report$SPR_msy
OM_report$YPR_msy
######
## Find deterministic Fmsy
## for simulation scenarios
## 20 years of over fishing 
## then fish at Fmsy for 40 years
######
fishery_sel = as.numeric(OM_report$fishery_selectivity)
waa = as.numeric(OM_report$stockMeanWeight[,60])
waa_catch = as.numeric(OM_report$catchMeanWeight[,60])

paa = as.numeric(TMB_data$propMat[,60])
ages = TMB_data$ages
M = TMB_data$natMor
Ninit = OM_report$N[,1]

F_30 = find_F_percent(target_spr = 30, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_35 = find_F_percent(target_spr = 35, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_40 = find_F_percent(target_spr = 40, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_max = find_F_max(fishery_sel, M, waa_catch, ages)
F_0.1 = find_F_0.1(fishery_sel, M, waa_catch, ages, derivative_method = 1)
F_0.1_alt = find_F_0.1(fishery_sel, M, waa_catch, ages, derivative_method = 2)
F_msy = find_F_msy(fishery_sel = fishery_sel, M = M, catch_waa = waa_catch, ssb_waa = waa, paa = paa, ages = ages, Ninit = OM_report$equilibrium_at_age, plus_group = T, prop_Z = 0.5, n_runs = 200)

###############
## Simulation
## Ramp first 20 years from 0.001 to  1.5 * F_msy
## 40 years at F_40
###############
## OM Fs
fishery_probs = c(seq(from = 0.01, to = F_msy$F_msy * 1.5, length = 20), rep(F_40$F_ref, 40))

png(filename = file.path(fig_dir, "F_pattern.png"), width = 6, height = 5, res = 250, units = "in")
plot(years, fishery_probs, type = "l", xlab = "Years", ylab = "F", lty = 2, lwd = 2, ylim = c(0,0.55))
dev.off()


OM_pars$ln_F = array(log(fishery_probs), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
# Set recruit devs at zero
OM_pars$ln_ycs_est =  rep(0,sum(TMB_data$ycs_estimated))# rnorm(sum(TMB_data$ycs_estimated),  0, this_bio$sigma_r)
TMB_data$F_method = 0
## re run the OM with true values so 
## we can check the log-likelihood calculations
OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_report = OM_obj$report()

# Plot SSB with deterministic recruitment
plot(full_years, OM_report$ssb / OM_report$B0 * 100, type = "l", lwd = 3, xlab = "Years", ylab = "Depletion", ylim = c(0,100));
# plot catches
plot(years, OM_report$pred_catches / 1000, type = "l", lwd = 3, xlab = "Years", ylab = "Depletion", ylim = c(0,10000));



## look at 100 simulations with stochastic recruitment
n_sims = 100
SSB_df = recruit_df = depletion_df = NULL
for(sim_iter in 1:n_sims) {
  ## simulate recruitment
  OM_pars$ln_ycs_est = rnorm(sum(TMB_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
  ## run model
  OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
  OM_report = OM_obj$report()
  ## save SSB
  this_ssb = data.frame(years = full_years, SSB = OM_report$ssb, sim_iter = sim_iter)
  SSB_df = rbind(SSB_df, this_ssb)
  ## this depletion
  this_depletion = data.frame(years = full_years, depletion = OM_report$depletion, sim_iter = sim_iter)
  depletion_df = rbind(depletion_df, this_depletion)
  ## save reccruitment
  this_recruit = data.frame(years = years, ycs = OM_report$ycs, sim_iter = sim_iter)
  recruit_df = rbind(recruit_df, this_recruit)
}

ggplot(SSB_df, aes(x = years, y = SSB, group = sim_iter)) +
  geom_line(linewidth = 1.1, alpha = 0.3) +
  theme_bw()+
  ylim(0, NA) +
  labs(x = "Year", y = "SSB")

ggplot(depletion_df, aes(x = years, y = depletion, group = sim_iter)) +
  geom_line(linewidth = 1.1, alpha = 0.3)+
  theme_bw()+
  ylim(0, NA) +
  labs(x = "Year", y = "Depletion")

ggplot(recruit_df, aes(x = years, y = ycs, group = sim_iter)) +
  geom_line(linewidth = 1.1, alpha = 0.3) +
  theme_bw()+
  ylim(0, NA) +
  labs(x = "Year", y = "Recruitment multiplier")

#############
## What about with non-equilibrium
## i.e. recruitment variation
#############
TMB_data$estimate_init_age_devs = 1
TMB_data$n_init_age_devs = 50
OM_pars$ln_sigma_init_age_devs = log(0.5)
OM_pars$ln_init_age_devs = rnorm(TMB_data$n_init_age_devs, 0, exp(OM_pars$ln_sigma_init_age_devs))
OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
OM_report = OM_obj$report()
OM_report$B0
OM_report$Binit

plot(ages, OM_report$N[,1], type = "l", lwd = 2, lty = 1)
lines(ages, OM_report$equilibrium_at_age, col = "red", lwd = 2, lty = 2)

## look at 100 simulations with stochastic recruitment
n_sims = 100
SSB_df = recruit_df = depletion_df = NULL
for(sim_iter in 1:n_sims) {
  ## simulate initial age deviations
  OM_pars$ln_init_age_devs = rnorm(TMB_data$n_init_age_devs, 0, exp(OM_pars$ln_sigma_init_age_devs))
  ## simulate recruitment
  OM_pars$ln_ycs_est = rnorm(sum(TMB_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
  ## run model
  OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
  OM_report = OM_obj$report()
  ## save SSB
  this_ssb = data.frame(years = full_years, SSB = OM_report$ssb, sim_iter = sim_iter)
  SSB_df = rbind(SSB_df, this_ssb)
  ## this depletion
  this_depletion = data.frame(years = full_years, depletion = OM_report$depletion, sim_iter = sim_iter)
  depletion_df = rbind(depletion_df, this_depletion)
  ## save reccruitment
  this_recruit = data.frame(years = years, ycs = OM_report$ycs, sim_iter = sim_iter)
  recruit_df = rbind(recruit_df, this_recruit)
}

ggplot(SSB_df, aes(x = years, y = SSB, group = sim_iter)) +
  geom_line(linewidth = 1.1, alpha = 0.3) +
  theme_bw()+
  ylim(0, NA) +
  labs(x = "Year", y = "SSB")

ggplot(depletion_df, aes(x = years, y = depletion, group = sim_iter)) +
  geom_line(linewidth = 1.1, alpha = 0.3)+
  theme_bw()+
  ylim(0, NA) +
  labs(x = "Year", y = "Depletion")

ggplot(recruit_df, aes(x = years, y = ycs, group = sim_iter)) +
  geom_line(linewidth = 1.1, alpha = 0.3) +
  theme_bw()+
  ylim(0, NA) +
  labs(x = "Year", y = "Recruitment multiplier")


## look at 100 simulations with stochastic recruitment
n_sims = 100
start_year = seq(from = min(years), to = min(years) + 20, by = 5)
#start_year = min(years) ## self test
is_self_test = T
SSB_df = recruit_df = depletion_df = NULL
survey_year_obs = years
fishery_year_obs = years
sim_survey_year_obs = years
sim_fishery_year_obs = years
TMB_data$survey_year_indicator = as.integer(TMB_data$years %in% survey_year_obs)
TMB_data$fishery_year_indicator = array(1, dim = c(length(TMB_data$years), TMB_data$n_fisheries))
## alter observation containers
TMB_data$survey_obs = rep(0, sum(TMB_data$survey_year_indicator))
TMB_data$survey_cv = rep(0.15, sum(TMB_data$survey_year_indicator))
TMB_data$survey_AF_obs = matrix(5, nrow = TMB_data$n_ages, ncol = sum(TMB_data$survey_year_indicator))
TMB_data$fishery_AF_obs = array(5, dim = c(TMB_data$n_ages, sum(TMB_data$fishery_year_indicator), TMB_data$n_fisheries))
convergence = matrix(T, ncol = length(start_year), nrow = n_sims)
estimate_ycs_after = 1961
MLE_pars_df = fishing_mortality_df = full_nll_df = SSB_df = recruit_df = depletion_df = reference_df = simple_metrics = OM_SSB_df = survey_abundance_df = survey_select_df = fishery_select_df = survey_AF_df = survey_mean_age_df =  fishery_AF_df = fishery_mean_age_df = NULL
## don't estimate the last 5 year class parameters
TMB_data$ycs_estimated[(length(TMB_data$ycs_estimated) - 9):length(TMB_data$ycs_estimated)] = 0
## simulate recruitment
OM_pars$ln_ycs_est = rnorm(sum(TMB_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
stochastic_recruitmnet = F ## does each simulation generate a new recruitment trend
TMB_data$n_init_age_devs = 1
OM_pars$ln_init_age_devs = 0
one_yr = two_yr = three_yr = four_yr = five_yr = six_yr = seven_yr = eight_yr = NULL
#'
#' Algorithm
#' iterate over each start year EM
#' Simulate the same OM for each EM n_sims times-set the seed
#'
#'
sim_start_time = Sys.time()
for(start_year_ndx in 1:length(start_year)) {
  cat("trialling start year ", start_year[start_year_ndx], "\n")
  EM_data = TMB_data
  EM_pars = OM_pars
  EM_data$years = start_year[start_year_ndx]:max(TMB_data$years)
  EM_year_ndx = which(TMB_data$years %in% EM_data$years)
  full_years = (min(EM_data$years) - 1):max(EM_data$years)
  
  EM_data$n_years = length(EM_data$years)
  sim_survey_year_obs = sim_survey_year_obs[sim_survey_year_obs %in% EM_data$years]
  sim_fishery_year_obs = sim_fishery_year_obs[sim_fishery_year_obs %in% EM_data$years]
  
  EM_data$survey_year_indicator = as.integer(EM_data$years %in% sim_survey_year_obs)
  EM_data$survey_obs = rep(10, sum(EM_data$survey_year_indicator))
  EM_data$survey_cv = rep(0.15, sum(EM_data$survey_year_indicator))
  EM_data$survey_AF_obs = matrix(10, nrow = EM_data$n_ages, ncol = sum(EM_data$survey_year_indicator))
  EM_data$fishery_year_indicator = array(as.integer(EM_data$years %in% sim_fishery_year_obs), dim = c(EM_data$n_years, EM_data$n_fisheries))
  EM_data$fishery_AF_obs = array(10, dim = c(EM_data$n_ages, length(sim_fishery_year_obs), EM_data$n_fisheries))
  EM_data$catches = array(1000, dim = c(EM_data$n_years, EM_data$n_fisheries))# this will be overriden in the simulate() call
  EM_data$catch_indicator = array(1, dim = c(EM_data$n_years, EM_data$n_fisheries))
  EM_data$ycs_estimated = TMB_data$ycs_estimated[EM_year_ndx]
  EM_data$ycs_estimated[which(EM_data$years < estimate_ycs_after)] = 0
  EM_data$catchMeanLength = EM_data$stockMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = EM_data$n_years, nrow = EM_data$n_ages)
  EM_data$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = EM_data$n_years, nrow = EM_data$n_ages)
  EM_data$propZ_ssb = rep(0.5, EM_data$n_years)
  EM_data$propZ_survey = rep(0.5, EM_data$n_years)
  EM_data$F_method = 1
  ## reset EM pars
  EM_pars$ln_ycs_est =  rnorm(sum(EM_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
  EM_pars$ln_F = OM_pars$ln_F 
  EM_pars$ln_init_age_devs = 0
  EM_data$n_init_age_devs = 1
  na_map = fix_pars(par_list = EM_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_extra_survey_cv", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
  #"logit_survey_ato95", "logit_f_ato95",
  for(sim_iter in 1:n_sims) {
    if(sim_iter %% 10 == 0)
      cat("simulation iterator ", sim_iter, "\n")
    
    ## keep the seed constant among start year runs so they are comparable
    ## also allows us to restart the simulation if one fails.
    set.seed(sim_iter)
    ## simulate initial age deviations
    ## simulate recruitment
    if(stochastic_recruitmnet) {
      OM_pars$ln_ycs_est = rnorm(sum(TMB_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
      OM_pars$ln_init_age_devs = rnorm(TMB_data$n_init_age_devs, 0, exp(OM_pars$ln_sigma_init_age_devs))
    }
    ## 
    #EM_pars$ln_ycs_est = OM_pars$ln_ycs_est[which(EM_data$ycs_estimated == 1)] ## fix at true values
    #EM_pars$logit_f_a50 =  OM_pars$logit_f_a50
    #EM_pars$logit_f_ato95 = OM_pars$logit_f_ato95
    #EM_pars$logit_survey_ato95 = OM_pars$logit_survey_ato95
    
    #EM_pars$logit_surveyQ =  OM_pars$logit_surveyQ
    ## if self-test don't do the initial age-devs or F-init
    if(is_self_test) {
      TMB_data$estimate_F_init = 0
      TMB_data$estimate_init_age_devs = 0
    }
    ## Build OM model
    OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    ## Simulate data
    OM_sim = OM_obj$simulate(complete = T)
    ## OM_true values
    OM_sim$F_method = 1
    OM_true_obj <- MakeADFun(OM_sim, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    OM_report = OM_true_obj$report()
    OM_nll = data.frame(total = OM_true_obj$fn(), catch = OM_report$catch_nll, recruitment = OM_report$recruit_nll, survey_ndx = OM_report$survey_index_nll, survey_AF = OM_report$survey_comp_nll, fishery_AF = OM_report$fishery_comp_nll)
    OM_nll$model = "OM"
    OM_nll$start_year = start_year[start_year_ndx]
    OM_nll$sim_iter = sim_iter
    ## save some OM information
    this_om_ssb = data.frame(years = full_years, SSB = c(OM_sim$ssb[1],OM_sim$ssb[EM_year_ndx]), sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    OM_SSB_df = rbind(OM_SSB_df, this_om_ssb)
    ## Fill EM data with simulated data
    
    EM_data$catches = matrix(OM_sim$catches[EM_year_ndx,], ncol = 1)
    #EM_data$catches = matrix(OM_sim$pred_catches[EM_year_ndx,], ncol = 1)
    
    EM_data$survey_obs = OM_sim$survey_obs[which(survey_year_obs %in% sim_survey_year_obs)]
    EM_data$survey_AF_obs = OM_sim$survey_AF_obs[,which(survey_year_obs %in% sim_survey_year_obs)]
    EM_data$fishery_AF_obs = array(OM_sim$fishery_AF_obs[, which(fishery_year_obs %in% sim_fishery_year_obs), 1], dim = c(length(EM_data$ages), length(sim_fishery_year_obs), 1))
    ## for debugging you can change simulated obs with the true values.
    use_expected_vals_for_observed = F
    if(use_expected_vals_for_observed) {
      EM_data$catches = matrix(OM_sim$pred_catches[EM_year_ndx,], ncol = 1)
      EM_data$survey_obs = OM_sim$survey_index_fitted[which(survey_year_obs %in% sim_survey_year_obs)]
      EM_data$survey_AF_obs = OM_sim$survey_AF_fitted[,which(survey_year_obs %in% sim_survey_year_obs)]
      EM_data$fishery_AF_obs = array(OM_sim$fishery_AF_fitted[, which(fishery_year_obs %in% sim_fishery_year_obs), 1], dim = c(length(EM_data$ages), length(sim_fishery_year_obs), 1))
    }
    
    EM_data$estimate_F_init = 0
    EM_data$estimate_init_age_devs = 0
    
    ## build EM
    EM_obj <- MakeADFun(EM_data, EM_pars, map = na_map, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM_obj)
    mle_spatial = nlminb(start = EM_obj$par, objective = EM_obj$fn, gradient  = EM_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM_obj$gr(mle_spatial$par))
                               h = optimHess(mle_spatial$par, fn = EM_obj$fn, gr = EM_obj$gr)
                               mle_spatial$par = mle_spatial$par - solve(h,g)
                               mle_spatial$objective = EM_obj$fn(mle_spatial$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    ## extract estimate pars and save them in data frame
    tmp_est_pars = data.frame(label = names(mle_spatial$par), value = mle_spatial$par, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    MLE_pars_df = rbind(MLE_pars_df, tmp_est_pars)
    EM_rep = EM_obj$report(mle_spatial$par)
    tmp_Bzero = data.frame(B0 = EM_rep$B0,Binit = EM_rep$Binit, R0 = EM_rep$R0, terminal_depletion = EM_rep$ssb[length(EM_rep$ssb)] / EM_rep$B0 * 100, sim_iter = sim_iter, start_year = start_year[start_year_ndx], Finit = EM_rep$F_init)
    simple_metrics = rbind(simple_metrics, tmp_Bzero)
    ## save some derived quantities
    ## save SSB
    this_ssb = data.frame(years = full_years, OM_SSB = c(OM_sim$ssb[1], OM_sim$ssb[EM_year_ndx + 1]), EM_SSB = EM_rep$ssb, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    SSB_df = rbind(SSB_df, this_ssb)
    ## this depletion
    this_depletion = data.frame(years = full_years, OM_depletion = c(OM_sim$depletion[1], OM_sim$depletion[EM_year_ndx + 1]), EM_depletion = EM_rep$depletion, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    depletion_df = rbind(depletion_df, this_depletion)
    ## save recruitment
    this_recruit = data.frame(years = EM_data$years, OM_ycs = OM_sim$ycs[EM_year_ndx], EM_ycs = EM_rep$ycs, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    recruit_df = rbind(recruit_df, this_recruit)
    ## selectivities
    this_srv_sel = data.frame(ages = ages, EM_selectivity = EM_rep$survey_selectivity, OM_selectivity = OM_sim$survey_selectivity, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    survey_select_df = rbind(survey_select_df, this_srv_sel)
    this_fish_sel = data.frame(ages = ages, EM_selectivity = EM_rep$fishery_selectivity, OM_selectivity = OM_sim$fishery_selectivity, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    fishery_select_df = rbind(fishery_select_df, this_fish_sel)
    ## fishing mortality
    tmp_F_df = data.frame(years = EM_data$years, EM_Fs = as.numeric(EM_rep$annual_Fs), OM_Fs = as.numeric(OM_sim$annual_Fs[1, EM_year_ndx]), sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    fishing_mortality_df = rbind(fishing_mortality_df, tmp_F_df)
    ## save reference points
    this_ref_df = data.frame(F30 = EM_rep$F30, F35 = EM_rep$F35, F40 = EM_rep$F40, F_0_1 = EM_rep$F_0_1, Fmsy = EM_rep$Fmsy, Fmax = EM_rep$Fmax, Bmsy = EM_rep$B_msy, B0 = EM_rep$B0, terminal_depletion = EM_rep$ssb[length(EM_rep$ssb)] / EM_rep$B0 * 100, 
                             sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    reference_df = rbind(reference_df, this_ref_df)
    
    ## survey abundance
    tmp_srv_abundance = data.frame(years = sim_survey_year_obs, SE = EM_rep$survey_sd, fitted = EM_rep$survey_index_fitted, observed = EM_data$survey_obs, OM_fitted = OM_sim$survey_index_fitted[which(survey_year_obs %in% sim_survey_year_obs)], sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    survey_abundance_df = rbind(survey_abundance_df, tmp_srv_abundance)
    ## survey Comp
    dimnames(EM_data$survey_AF_obs) = list(ages, sim_survey_year_obs)
    molten_fish_AF = reshape2::melt(EM_data$survey_AF_obs)
    molten_fish_AF_fitted = reshape2::melt(EM_rep$survey_AF_fitted)
    colnames(molten_fish_AF) = c("age", "year", "observed")
    molten_fish_AF$fitted = molten_fish_AF_fitted$value
    molten_fish_AF = molten_fish_AF %>% group_by(year) %>% mutate(Neff = sum(observed), observed_prop = observed / Neff)
    molten_fish_AF$start_year = start_year[start_year_ndx]
    molten_fish_AF$sim_iter = sim_iter
    ## get mean age values
    mean_age_AF= molten_fish_AF %>% ungroup() %>% group_by(year) %>% summarise(Ey = sum(age * fitted), Oy = sum(age * observed_prop), E_squared_y = sum(age^2 * fitted), Neff = mean(Neff)) %>% ungroup()
    mean_age_AF$Ry = mean_age_AF$Oy - mean_age_AF$Ey
    mean_age_AF$SEy = sqrt((mean_age_AF$E_squared_y - mean_age_AF$Ey^2) / mean_age_AF$Neff)
    mean_age_AF$'Std.res' <- (mean_age_AF$Oy - mean_age_AF$Ey)/mean_age_AF$SEy
    ## I think this is the final Francis weighting value TODO: to check
    Nmult <- 1 / var(mean_age_AF$'Std.res',na.rm=TRUE)
    # Find the adjusted confidence intervals
    mean_age_AF$ObsloAdj <- mean_age_AF$Oy - 2 * mean_age_AF$SEy / sqrt(Nmult)
    mean_age_AF$ObshiAdj <- mean_age_AF$Oy + 2 * mean_age_AF$SEy / sqrt(Nmult)
    mean_age_AF$start_year = start_year[start_year_ndx]
    mean_age_AF$sim_iter = sim_iter
    survey_AF_df = rbind(survey_AF_df, molten_fish_AF)
    survey_mean_age_df = rbind(survey_mean_age_df, mean_age_AF)
    
    ## Fishery Comp
    dimnames(EM_data$fishery_AF_obs) = list(ages, sim_fishery_year_obs)
    molten_fish_AF = reshape2::melt(EM_data$fishery_AF_obs)
    molten_fish_AF_fitted = reshape2::melt(EM_rep$fishery_AF_fitted)
    colnames(molten_fish_AF) = c("age", "year", "fishery_ndx", "observed")
    molten_fish_AF$fitted = molten_fish_AF_fitted$value
    molten_fish_AF = molten_fish_AF %>% group_by(year, fishery_ndx) %>% mutate(Neff = sum(observed), observed_prop = observed / Neff)
    molten_fish_AF$start_year = start_year[start_year_ndx]
    molten_fish_AF$sim_iter = sim_iter
    ## get mean age values
    mean_age_AF= molten_fish_AF %>% ungroup() %>% group_by(year, fishery_ndx) %>% summarise(Ey = sum(age * fitted), Oy = sum(age * observed_prop), E_squared_y = sum(age^2 * fitted), Neff = mean(Neff)) %>% ungroup()
    mean_age_AF$Ry = mean_age_AF$Oy - mean_age_AF$Ey
    mean_age_AF$SEy = sqrt((mean_age_AF$E_squared_y - mean_age_AF$Ey^2) / mean_age_AF$Neff)
    mean_age_AF$'Std.res' <- (mean_age_AF$Oy - mean_age_AF$Ey)/mean_age_AF$SEy
    ## I think this is the final Francis weighting value TODO: to check
    Nmult <- 1 / var(mean_age_AF$'Std.res',na.rm=TRUE)
    # Find the adjusted confidence intervals
    mean_age_AF$ObsloAdj <- mean_age_AF$Oy - 2 * mean_age_AF$SEy / sqrt(Nmult)
    mean_age_AF$ObshiAdj <- mean_age_AF$Oy + 2 * mean_age_AF$SEy / sqrt(Nmult)
    mean_age_AF$start_year = start_year[start_year_ndx]
    mean_age_AF$sim_iter = sim_iter
    fishery_AF_df = rbind(fishery_AF_df, molten_fish_AF)
    fishery_mean_age_df = rbind(fishery_mean_age_df, mean_age_AF)
    
    # negative log likelihood
    EM_nll = data.frame(total = mle_spatial$objective, catch = EM_rep$catch_nll, recruitment = EM_rep$recruit_nll, survey_ndx = EM_rep$survey_index_nll, survey_AF = EM_rep$survey_comp_nll, fishery_AF = EM_rep$fishery_comp_nll)
    EM_nll$model = "EM"
    EM_nll$start_year = start_year[start_year_ndx]
    EM_nll$sim_iter = sim_iter
    full_nll_df = rbind(full_nll_df, EM_nll, OM_nll)
    
    ## look at specific cohorts over time for biases. This may help us understand what
    ## dynamic is causing issues with SSB.
    one_yr_tmp = data.frame(OM_numbers = c(OM_report$N[1, 1], OM_report$N[1,EM_year_ndx + 1]), EM_numbers = EM_rep$N[1,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    one_yr = rbind(one_yr, one_yr_tmp)
    two_yr_tmp = data.frame(OM_numbers = c(OM_report$N[2, 1], OM_report$N[2,EM_year_ndx + 1]), EM_numbers = EM_rep$N[2,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    two_yr = rbind(two_yr, two_yr_tmp)
    three_yr_tmp = data.frame(OM_numbers = c(OM_report$N[3, 1], OM_report$N[3,EM_year_ndx + 1]), years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    three_yr = rbind(three_yr, three_yr_tmp)
    four_yr_tmp = data.frame(OM_numbers = c(OM_report$N[4, 1], OM_report$N[4,EM_year_ndx + 1]), EM_numbers = EM_rep$N[4,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    four_yr = rbind(four_yr, four_yr_tmp)
    five_yr_tmp = data.frame(OM_numbers = c(OM_report$N[5, 1], OM_report$N[5,EM_year_ndx + 1]), EM_numbers = EM_rep$N[5,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    five_yr = rbind(five_yr, five_yr_tmp)
    six_yr_tmp = data.frame(OM_numbers = c(OM_report$N[6, 1], OM_report$N[6,EM_year_ndx + 1]), EM_numbers = EM_rep$N[6,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    six_yr = rbind(six_yr, six_yr_tmp)
    seven_yr_tmp = data.frame(OM_numbers = c(OM_report$N[7, 1], OM_report$N[7,EM_year_ndx + 1]), EM_numbers = EM_rep$N[7,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    seven_yr = rbind(seven_yr, seven_yr_tmp)
    eight_yr_tmp = data.frame(OM_numbers = c(OM_report$N[8, 1], OM_report$N[8,EM_year_ndx + 1]), EM_numbers = EM_rep$N[8,], years = full_years, sim_iter = sim_iter, start_year = start_year[start_year_ndx])
    eight_yr = rbind(eight_yr, eight_yr_tmp)
  }
}
cat("Simulation run time = ", Sys.time() - sim_start_time, " minutes\n")
unique(names(EM_obj$par))

## look at fit to AFs
plot(ages, EM_data$survey_AF_obs[,1]/colSums(EM_data$survey_AF_obs)[1], type = "p")
lines(ages, EM_rep$survey_AF_fitted[,1], type = "l", col = "red", lwd = 2)
plot(ages, EM_data$survey_AF_obs[,2]/colSums(EM_data$survey_AF_obs)[2], type = "p")
lines(ages, EM_rep$survey_AF_fitted[,2], type = "l", col = "red", lwd = 2)
plot(ages, EM_data$survey_AF_obs[,3]/colSums(EM_data$survey_AF_obs)[23], type = "p")
lines(ages, EM_rep$survey_AF_fitted[,3], type = "l", col = "red", lwd = 2)
plot(ages, EM_data$survey_AF_obs[,4]/colSums(EM_data$survey_AF_obs)[4], type = "p")
lines(ages, EM_rep$survey_AF_fitted[,4], type = "l", col = "red", lwd = 2)

##########
# Look at Biases and Derived quantities
##########
ggplot(full_nll_df %>% pivot_longer(cols = c("total","catch","recruitment","survey_ndx","survey_AF","fishery_AF"), names_to ="component")) +
  geom_boxplot(aes(x = model, y = value, col = model, fill = model)) +
  facet_wrap(~component, scales = "free_y") +
  theme_bw()

## Plot SSB
ggplot() +
  geom_line(data = SSB_df, aes(x = years, y = EM_SSB, group = sim_iter, col = "EM"), linewidth = 1.1, alpha = 0.2) +
  geom_line(data = SSB_df, aes(x = years, y = OM_SSB, group = sim_iter, col = "OM"), linewidth = 1.1, alpha = 0.2) +
  ylim(0,NA) +
  theme_bw()
# Plot RE in SSB
SSB_df = SSB_df %>% group_by(years, sim_iter, start_year) %>% mutate(RE = (OM_SSB - EM_SSB) / OM_SSB * 100)
ggplot(data = SSB_df, aes(x = years, y = RE, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  #ylim(-50,50) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw();

# Plot RE for estimated 1 yr olds
one_yr = one_yr %>% group_by(years, sim_iter, start_year) %>% mutate(RE = (OM_numbers - EM_numbers) / OM_numbers * 100)
ggplot(data = one_yr, aes(x = years, y = RE, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  #ylim(-50,50) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw();
two_yr = two_yr %>% group_by(years, sim_iter, start_year) %>% mutate(RE = (OM_numbers - EM_numbers) / OM_numbers * 100)
ggplot(data = two_yr, aes(x = years, y = RE, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  #ylim(-50,50) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw();

# Plot RE in Fs
fishing_mortality_df = fishing_mortality_df %>% group_by(years, sim_iter, start_year) %>% mutate(RE = (OM_Fs - EM_Fs) / OM_Fs * 100)
ggplot(data = fishing_mortality_df, aes(x = years, y = RE, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw();

## Plot Recruitment
ggplot() +
  geom_line(data = recruit_df, aes(x = years, y = EM_ycs, group = sim_iter, col = "EM"), linewidth = 1.1, alpha = 0.2) +
  geom_line(data = recruit_df, aes(x = years, y = OM_ycs, group = sim_iter, col = "OM"), linewidth = 1.1, alpha = 0.2) +
  theme_bw()

recruit_df = recruit_df %>% group_by(years, sim_iter, start_year) %>% mutate(RE = (OM_ycs - EM_ycs) / OM_ycs * 100)
# Plot RE in Recruitment
ggplot(data = recruit_df, aes(x = years, y = RE, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()
## Fishery selectivity
ggplot(fishery_select_df, aes(x = ages, y = (OM_selectivity - EM_selectivity) / OM_selectivity * 100, group = ages)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()

## Survey selectivity
ggplot(survey_select_df, aes(x = ages, y = (OM_selectivity - EM_selectivity) / OM_selectivity * 100, group = ages)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()
## Plot survey index
survey_abundance_df$std_resids = (log(survey_abundance_df$observed/survey_abundance_df$fitted) + 0.5*survey_abundance_df$SE^2)/survey_abundance_df$SE
ggplot(data = survey_abundance_df, aes(x = years, y = std_resids, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  labs(title = "Survey index", x = "Years", y = "Pearson residuals") +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()

## Plot survey mean age resids
ggplot(data = survey_mean_age_df, aes(x = year, y = Std.res, group = year)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  labs(title = "Survey mean age index", x = "Years", y = "Pearson residuals") +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()
## Plot fishery mean age resids
ggplot(data = fishery_mean_age_df, aes(x = year, y = Std.res, group = year)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  labs(title = "Fishery mean age index", x = "Years", y = "Pearson residuals") +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()

## plot aggregated survey AF
summarise_obs = survey_AF_df %>% group_by(age, sim_iter) %>% summarise(aggregated_observed = sum(observed), aggregrated_expected = sum(Neff * fitted)) %>% ungroup() %>%
  group_by(age) %>% mutate(aggregated_resid = aggregated_observed - aggregrated_expected)

ggplot(data = summarise_obs, aes(x = age, y = aggregated_resid, group= age)) +
  geom_boxplot() +
  labs(title = "Survey aggregated residuals", x = "Years", y = "Aggregated residuals") +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()

## plot aggregated survey AF
summarise_obs = fishery_AF_df %>% group_by(age, sim_iter) %>% summarise(aggregated_observed = sum(observed), aggregrated_expected = sum(Neff * fitted)) %>% ungroup() %>%
  group_by(age) %>% mutate(aggregated_resid = aggregated_observed - aggregrated_expected)

ggplot(data = summarise_obs, aes(x = age, y = aggregated_resid, group= age)) +
  geom_boxplot() +
  labs(title = "Fishery aggregated residuals", x = "Years", y = "Aggregated residuals") +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw()

## look at simulated abundance vs predicted
ggplot(survey_abundance_df, aes(x = years)) +
  geom_boxplot(aes(y = observed, group = years)) +
  geom_line(linewidth = 1.1, aes(y = OM_fitted, col = "Expected")) +
  theme_bw()
ggplot(survey_abundance_df, aes(x = years)) +
  geom_boxplot(aes(y = fitted, group = years)) +
  geom_line(linewidth = 1.1, aes(y = OM_fitted, col = "Expected")) +
  theme_bw()
### Reference points
ggplot(reference_df, aes(x = start_year, y = F30, group = start_year)) +
  geom_boxplot()

ggplot(reference_df, aes(x = start_year, y = F35, group = start_year)) +
  geom_boxplot()

ggplot(reference_df, aes(x = start_year, y = F40, group = start_year)) +
  geom_boxplot()

ggplot(reference_df, aes(x = start_year, y = F_0_1, group = start_year)) +
  geom_boxplot()

ggplot(reference_df, aes(x = start_year, y = Fmsy, group = start_year)) +
  geom_boxplot()

ggplot(reference_df, aes(x = start_year, y = B0, group = start_year)) +
  geom_boxplot() +
  geom_hline(yintercept = OM_report$B0, linetype = "dashed", col ="red")# +
#ylim(c(1e7, 2e8))

ggplot(reference_df, aes(x = start_year, y = terminal_depletion, group = start_year)) +
  geom_boxplot()

ggplot(reference_df, aes(x = start_year, y = Fmax, group = start_year)) +
  geom_boxplot()

