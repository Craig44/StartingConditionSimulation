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
library(purrr)

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

TMB_data$rec_devs_sum_to_zero = 0
TMB_data$Q_r_for_sum_to_zero = Q_sum_to_zero_QR(length(TMB_data$years))
## iniital fishery_probs

fishery_probs = c(seq(from = 0.01, to = this_bio$M * 1.5, length = 20), rep(this_bio$M, 40) )
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
OM_pars$ln_sigma_init_age_devs = OM_pars$ln_sigma_r
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

OM_label = "Medium_OM1"
output_data = file.path(DIR$data, OM_label)
if(!dir.exists(output_data))
  dir.create(output_data)
OM_fig_dir = file.path(fig_dir, OM_label)
if(!dir.exists(OM_fig_dir))
  dir.create(OM_fig_dir)

# - OM1 - equilibrium N1
# - OM2 - N1 X dev_a
# - OM3 - N1 X exp(-(M + F_init))
# - OM4 - N1 X exp(-(M + F_init)) X dev_a

#'
#' Algorithm
#' iterate over each start year EM
#' Simulate the same OM for each EM n_sims times-set the seed
#'
#'
sim_start_time = Sys.time()
mle_lst_EM1 = mle_lst_EM2 = mle_lst_EM3 = mle_lst_EM4 = mle_lst_EM5 = mle_lst_EM6 = mle_lst_EM7 = list()
for(start_year_ndx in 1:length(start_year)) {
  cat("trialling start year ", start_year[start_year_ndx], "\n")
  EM5_data = EM6_data = EM_data = TMB_data
  EM_pars = OM_pars
  EM_data$years = start_year[start_year_ndx]:max(TMB_data$years)
  EM_year_ndx = which(TMB_data$years %in% EM_data$years)
  EM_excluded_year_ndx = which(!TMB_data$years %in% EM_data$years)
  
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
  ## taylor the EM5 
  EM5_data$survey_year_indicator = as.integer(EM5_data$years %in% sim_survey_year_obs)
  EM5_data$survey_obs = rep(10, sum(EM5_data$survey_year_indicator))
  EM5_data$survey_cv = rep(0.15, sum(EM5_data$survey_year_indicator))
  EM5_data$survey_AF_obs = matrix(10, nrow = EM5_data$n_ages, ncol = sum(EM5_data$survey_year_indicator))
  EM5_data$fishery_year_indicator = array(as.integer(EM5_data$years %in% sim_fishery_year_obs), dim = c(EM5_data$n_years, EM5_data$n_fisheries))
  EM5_data$fishery_AF_obs = array(10, dim = c(EM5_data$n_ages, length(sim_fishery_year_obs), EM5_data$n_fisheries))
  EM5_data$catches = array(1000, dim = c(EM5_data$n_years, EM5_data$n_fisheries))# this will be overriden in the simulate() call
  EM5_data$catch_indicator = array(1, dim = c(EM5_data$n_years, EM5_data$n_fisheries))
  EM5_data$ycs_estimated = TMB_data$ycs_estimated
  EM5_data$ycs_estimated[which(EM5_data$years < estimate_ycs_after)] = 0
  EM5_data$catchMeanLength = EM5_data$stockMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = EM5_data$n_years, nrow = EM5_data$n_ages)
  EM5_data$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = EM5_data$n_years, nrow = EM5_data$n_ages)
  EM5_data$propZ_ssb = rep(0.5, EM5_data$n_years)
  EM5_data$propZ_survey = rep(0.5, EM5_data$n_years)
  EM5_data$F_method = 1
  ## EM6 data inputs
  EM6_data = EM5_data
  EM7_data = EM5_data
  EM5_pars = EM6_pars = EM7_pars = EM_pars
  EM5_pars$ln_ycs_est =  rnorm(sum(EM5_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
  EM6_pars$ln_ycs_est =  rnorm(sum(EM5_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
  EM7_pars$ln_ycs_est =  rnorm(sum(EM5_data$ycs_estimated),  -0.5*this_bio$sigma_r^2, this_bio$sigma_r)
  
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
    EM5_data$catches = matrix(OM_sim$catches, ncol = 1)
    EM6_data$catches = matrix(OM_sim$catches, ncol = 1)
    EM7_data$catches = matrix(OM_sim$catches, ncol = 1)
    ## Doubel for EM5
    EM5_data$catches[EM_excluded_year_ndx,] = EM5_data$catches[EM_excluded_year_ndx,] * 1.25 ## +25%
    EM6_data$catches[EM_excluded_year_ndx,] = EM6_data$catches[EM_excluded_year_ndx,] * 0.75 ## -25%
    
    #EM_data$catches = matrix(OM_sim$pred_catches[EM_year_ndx,], ncol = 1)
    
    EM7_data$survey_obs = EM5_data$survey_obs = EM6_data$survey_obs = EM_data$survey_obs = OM_sim$survey_obs[which(survey_year_obs %in% sim_survey_year_obs)]
    EM7_data$survey_AF_obs =EM5_data$survey_AF_obs = EM6_data$survey_AF_obs = EM_data$survey_AF_obs = OM_sim$survey_AF_obs[,which(survey_year_obs %in% sim_survey_year_obs)]
    EM7_data$fishery_AF_obs = EM5_data$fishery_AF_obs = EM6_data$fishery_AF_obs = EM_data$fishery_AF_obs = array(OM_sim$fishery_AF_obs[, which(fishery_year_obs %in% sim_fishery_year_obs), 1], dim = c(length(EM_data$ages), length(sim_fishery_year_obs), 1))
    

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
    
    ## EM2
    ## estimate F-init only 
    EM2_data = EM_data
    EM2_data$estimate_F_init = 1
    EM2_data$estimate_init_age_devs = 0
    
    ## build EM 2
    na_map_EM2 = fix_pars(par_list = EM_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_extra_survey_cv", "ln_init_age_devs", "ln_sigma_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
    EM2_obj <- MakeADFun(EM2_data, EM_pars, map = na_map_EM2, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM2_obj)
    mle_spatial_EM2 = nlminb(start = EM2_obj$par, objective = EM2_obj$fn, gradient  = EM2_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM2_obj$gr(mle_spatial_EM2$par))
                               h = optimHess(mle_spatial_EM2$par, fn = EM2_obj$fn, gr = EM2_obj$gr)
                               mle_spatial_EM2$par = mle_spatial_EM2$par - solve(h,g)
                               mle_spatial_EM2$objective = EM2_obj$fn(mle_spatial_EM2$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation EM2 ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    EM2_rep = EM2_obj$report(mle_spatial_EM2$par)
    

    ## EM3
    ## estimate init age devs NOT F-init
    EM3_data = EM_data
    EM3_data$estimate_F_init = 0
    EM3_data$estimate_init_age_devs = 1
    EM3_data$n_init_age_devs = 25
    EM3_pars = EM_pars
    EM3_pars$ln_init_age_devs = rep(log(1), EM3_data$n_init_age_devs)
    na_map_EM3 = fix_pars(par_list = EM3_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_F_init","ln_extra_survey_cv", "ln_sigma_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
    ## build EM 2 
    EM3_obj <- MakeADFun(EM3_data, EM3_pars, map = na_map_EM3, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM3_obj)
    mle_spatial_EM3 = nlminb(start = EM3_obj$par, objective = EM3_obj$fn, gradient  = EM3_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM3_obj$gr(mle_spatial_EM3$par))
                               h = optimHess(mle_spatial_EM3$par, fn = EM3_obj$fn, gr = EM3_obj$gr)
                               mle_spatial_EM3$par = mle_spatial_EM3$par - solve(h,g)
                               mle_spatial_EM3$objective = EM3_obj$fn(mle_spatial_EM3$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation EM3 ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    EM3_rep = EM3_obj$report(mle_spatial_EM3$par)
    
    ## EM4
    ## estimate F-init and init age devs
    EM4_data = EM_data
    EM4_data$estimate_F_init = 1
    EM4_data$estimate_init_age_devs = 1
    EM4_data$n_init_age_devs = 25
    EM4_pars = EM_pars
    EM4_pars$ln_init_age_devs = rep(log(1), EM4_data$n_init_age_devs)
    na_map_EM4 = fix_pars(par_list = EM4_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_extra_survey_cv", "ln_sigma_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
    ## build EM 2 
    EM4_obj <- MakeADFun(EM4_data, EM4_pars, map = na_map_EM4, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM4_obj)
    mle_spatial_EM4 = nlminb(start = EM4_obj$par, objective = EM4_obj$fn, gradient  = EM4_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM4_obj$gr(mle_spatial_EM4$par))
                               h = optimHess(mle_spatial_EM4$par, fn = EM4_obj$fn, gr = EM4_obj$gr)
                               mle_spatial_EM4$par = mle_spatial_EM4$par - solve(h,g)
                               mle_spatial_EM4$objective = EM4_obj$fn(mle_spatial_EM4$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation EM4 ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    EM4_rep = EM4_obj$report(mle_spatial_EM4$par)
    
    ## EM4
    ## estimate F-init and init age devs
    EM5_data$estimate_F_init = 0
    EM5_data$estimate_init_age_devs = 0
    EM5_data$n_init_age_devs = 1
    na_map_EM5 = fix_pars(par_list = EM5_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_extra_survey_cv", "ln_sigma_init_age_devs", "ln_F_init", "ln_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
    ## build EM 2 
    EM5_obj <- MakeADFun(EM5_data, EM5_pars, map = na_map_EM5, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM5_obj)
    mle_spatial_EM5 = nlminb(start = EM5_obj$par, objective = EM5_obj$fn, gradient  = EM5_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM5_obj$gr(mle_spatial_EM5$par))
                               h = optimHess(mle_spatial_EM5$par, fn = EM5_obj$fn, gr = EM5_obj$gr)
                               mle_spatial_EM5$par = mle_spatial_EM5$par - solve(h,g)
                               mle_spatial_EM5$objective = EM5_obj$fn(mle_spatial_EM5$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation EM5 ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    EM5_rep = EM5_obj$report(mle_spatial_EM5$par)
    
    ## EM6
    ## estimate F-init and init age devs
    EM6_data$estimate_F_init = 0
    EM6_data$estimate_init_age_devs = 0
    EM6_data$n_init_age_devs = 1
    na_map_EM6 = fix_pars(par_list = EM6_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_extra_survey_cv", "ln_sigma_init_age_devs", "ln_F_init", "ln_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
    ## build EM 2 
    EM6_obj <- MakeADFun(EM6_data, EM6_pars, map = na_map_EM6, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM6_obj)
    mle_spatial_EM6 = nlminb(start = EM6_obj$par, objective = EM6_obj$fn, gradient  = EM6_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM6_obj$gr(mle_spatial_EM6$par))
                               h = optimHess(mle_spatial_EM6$par, fn = EM6_obj$fn, gr = EM6_obj$gr)
                               mle_spatial_EM6$par = mle_spatial_EM6$par - solve(h,g)
                               mle_spatial_EM6$objective = EM6_obj$fn(mle_spatial_EM6$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation EM6 ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    EM6_rep = EM6_obj$report(mle_spatial_EM6$par)
    ## EM7
    ## estimate F-init and init age devs
    EM7_data$estimate_F_init = 0
    EM7_data$estimate_init_age_devs = 0
    EM7_data$n_init_age_devs = 1
    na_map_EM7 = fix_pars(par_list = EM7_pars, pars_to_exclude = c("ln_sigma_r", "ln_F", "ln_catch_sd" ,"ln_extra_survey_cv", "ln_sigma_init_age_devs", "ln_F_init", "ln_init_age_devs"))#, "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
    ## build EM 2 
    EM7_obj <- MakeADFun(EM7_data, EM7_pars, map = na_map_EM7, DLL= "AgeStructuredModel", checkParameterOrder = F, silent = T)
    
    ## optimise
    #cbind(start_EM_rep$equilibrium_at_age, start_EM_rep$N[,1])
    #check_gradients(EM7_obj)
    mle_spatial_EM7 = nlminb(start = EM7_obj$par, objective = EM7_obj$fn, gradient  = EM7_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
    try_improve = tryCatch(expr =
                             for(i in 1:2) {
                               g = as.numeric(EM7_obj$gr(mle_spatial_EM7$par))
                               h = optimHess(mle_spatial_EM7$par, fn = EM7_obj$fn, gr = EM7_obj$gr)
                               mle_spatial_EM7$par = mle_spatial_EM7$par - solve(h,g)
                               mle_spatial_EM7$objective = EM7_obj$fn(mle_spatial_EM7$par)
                             }
                           , error = function(e){e})
    
    if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
      cat("Failed simulation EM7 ", sim_iter, "\n")
      convergence[sim_iter, start_year_ndx] = F
      next;
    }
    EM7_rep = EM7_obj$report(mle_spatial_EM7$par)
    ## savve output
    mle_lst_EM1[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM_rep
    mle_lst_EM2[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM2_rep
    mle_lst_EM3[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM3_rep
    mle_lst_EM4[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM4_rep
    mle_lst_EM5[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM5_rep
    mle_lst_EM6[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM6_rep
    mle_lst_EM7[[as.character(start_year[start_year_ndx])]][[as.character(sim_iter)]] = EM7_rep
    
  }
}
cat("Simulation run time = ", Sys.time() - sim_start_time, " minutes\n")
unique(names(EM_obj$par))
if(F) {
  saveRDS(object = mle_lst_EM1, file = file.path(output_data, "mle_lst_EM1.RDS"))
  saveRDS(object = mle_lst_EM2, file = file.path(output_data, "mle_lst_EM2.RDS"))
  saveRDS(object = mle_lst_EM3, file = file.path(output_data, "mle_lst_EM3.RDS"))
  saveRDS(object = mle_lst_EM4, file = file.path(output_data, "mle_lst_EM4.RDS"))
  saveRDS(object = mle_lst_EM5, file = file.path(output_data, "mle_lst_EM5.RDS"))
  saveRDS(object = mle_lst_EM6, file = file.path(output_data, "mle_lst_EM6.RDS"))
  saveRDS(object = mle_lst_EM7, file = file.path(output_data, "mle_lst_EM7.RDS"))
  
  saveRDS(object = OM_pars, file = file.path(output_data, "OM_pars.RDS"))
  saveRDS(object = OM_obj, file = file.path(output_data, "OM_obj.RDS"))
}
if(F) {
  mle_lst_EM1 = readRDS(file = file.path(output_data, "mle_lst_EM1.RDS"))
  mle_lst_EM2 = readRDS(file = file.path(output_data, "mle_lst_EM2.RDS"))
  mle_lst_EM3 = readRDS(file = file.path(output_data, "mle_lst_EM3.RDS"))
  mle_lst_EM4 = readRDS(file = file.path(output_data, "mle_lst_EM4.RDS"))
  
  mle_lst_EM5 = readRDS(file = file.path(output_data, "mle_lst_EM5.RDS"))
  mle_lst_EM6 = readRDS(file = file.path(output_data, "mle_lst_EM6.RDS"))  
  mle_lst_EM7 = readRDS(file = file.path(output_data, "mle_lst_EM7.RDS"))
  OM_pars = readRDS(file = file.path(output_data, "OM_pars.RDS"))
  OM_obj = readRDS(file = file.path(output_data, "OM_obj.RDS"))
  
  OM_rep = OM_obj$report(par = OM_obj$env$last.par.best)
  OM_sim = OM_obj$simulate(par = OM_obj$env$last.par.best, complete = T)
  
}
## save R objects
saveRDS(object = MLE_pars_df, file = file.path(output_data, "MLE_pars_df.RDS"))
saveRDS(object = fishing_mortality_df, file = file.path(output_data, "fishing_mortality_df.RDS"))
saveRDS(object = full_nll_df, file = file.path(output_data, "full_nll_df.RDS"))
saveRDS(object = SSB_df, file = file.path(output_data, "SSB_df.RDS"))
saveRDS(object = recruit_df, file = file.path(output_data, "recruit_df.RDS"))
saveRDS(object = reference_df, file = file.path(output_data, "reference_df.RDS"))
saveRDS(object = depletion_df, file = file.path(output_data, "depletion_df.RDS"))
saveRDS(object = simple_metrics, file = file.path(output_data, "simple_metrics.RDS"))
saveRDS(object = survey_abundance_df, file = file.path(output_data, "survey_abundance_df.RDS"))
saveRDS(object = survey_select_df, file = file.path(output_data, "survey_select_df.RDS"))
saveRDS(object = fishery_select_df, file = file.path(output_data, "fishery_select_df.RDS"))
saveRDS(object = survey_AF_df, file = file.path(output_data, "survey_AF_df.RDS"))
saveRDS(object = survey_mean_age_df, file = file.path(output_data, "survey_mean_age_df.RDS"))
saveRDS(object = fishery_AF_df, file = file.path(output_data, "fishery_AF_df.RDS"))
saveRDS(object = fishery_mean_age_df, file = file.path(output_data, "fishery_mean_age_df.RDS"))
saveRDS(object = one_yr, file = file.path(output_data, "one_yr.RDS"))
saveRDS(object = two_yr, file = file.path(output_data, "two_yr.RDS"))
saveRDS(object = three_yr, file = file.path(output_data, "three_yr.RDS"))
saveRDS(object = four_yr, file = file.path(output_data, "four_yr.RDS"))
saveRDS(object = five_yr, file = file.path(output_data, "five_yr.RDS"))
saveRDS(object = six_yr, file = file.path(output_data, "six_yr.RDS"))
saveRDS(object = seven_yr, file = file.path(output_data, "seven_yr.RDS"))
saveRDS(object = eight_yr, file = file.path(output_data, "eight_yr.RDS"))

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
  theme_bw() +
  facet_wrap(~start_year)
# Plot RE in SSB
SSB_df = SSB_df %>% group_by(years, sim_iter, start_year) %>% mutate(RE = (OM_SSB - EM_SSB) / OM_SSB * 100)
quant_RE_ssb = get_df_quantiles(SSB_df, group_vars = c("years", "start_year"), y_value = "RE", quants = c(0.025, 0.25, 0.4, 0.5, 0.6, 0.75, 0.975))
                             
ggplot(data = SSB_df, aes(x = years, y = RE, group = years)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  #ylim(-50,50) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme_bw();

ggplot(quant_RE_ssb, aes(x = years)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  facet_wrap(~start_year) +
  #ylim(-50,50) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="black", linewidth = 1.1) +
  theme_bw() +
  labs(y = "SSB", x = "Years")

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

##########
# Generate OM plots
##########
# initial n age
plot(TMB_data$ages, OM_sim$equilibrium_at_age, type = "l", lwd = 2, ylab = "Initial numbers", xlab = "Age")
lines(TMB_data$ages, OM_sim$N[,1], type = "l",  col = "red", lty = 2, lwd = 2)
# SSB
png(filename = file.path(OM_fig_dir, "ssb_start_years.png"), width = 7, height = 5, res = 250, units = "in")
plot(c(min(TMB_data$years)- 1,TMB_data$years), OM_sim$ssb/1000, type = "l", lwd = 2, ylab = "SSB (000's t)", xlab = "Years", ylim = c(0, 2000))
abline(v = start_year, lty = "dashed", col = "blue", lwd = 2)
legend('topright', col = c("black", "blue"), legend = c("OM", "EM start years"), lty = c(1,2), lwd = 3)
dev.off()
# recruitment
plot(TMB_data$years, OM_sim$ycs, type = "l", lwd = 2, ylab = "Year class strengths", xlab = "Years", ylim = c(0, 5))
# F
plot(TMB_data$years, OM_sim$annual_Fs, type = "l", lwd = 2, ylab = "Fully selected F", xlab = "Years", ylim = c(0, 0.25))
abline(h = F_msy$F_msy, lty = "dashed", col = "red", lwd =2)
legend('topright', col = c("black", "red"), legend = c("OM", "Fmsy"), lty = c(1,2), lwd = 3)

# Catch scenario
png(filename = file.path(OM_fig_dir, "Catch_example1.png"), width = 7.5, height = 5, res = 250, units = "in")
plot(c(TMB_data$years), OM_sim$catches/1000, type = "l", lwd = 3, ylab = "Catch", xlab = "Years", ylim = c(0, 120))
legend('topright', col = c("black", "gold", "purple", "red", "blue"), legend = c("OM", "EM1a", "EM1b", "EM2-4","Initial year"), lty = c(1,2, 2, 2,2), lwd = 3)
dev.off()

png(filename = file.path(OM_fig_dir, "Catch_example2.png"), width = 7.5, height = 5, res = 250, units = "in")
plot(c(TMB_data$years), OM_sim$catches/1000, type = "l", lwd = 3, ylab = "Catch", xlab = "Years", ylim = c(0, 120))
legend('topright', col = c("black", "gold", "purple", "red", "blue"), legend = c("OM", "EM1a", "EM1b", "EM2-4","Initial year"), lty = c(1,2, 2, 2,2), lwd = 3)
abline(v = start_year[3], lty = "dashed", col = "blue", lwd = 3)
polygon(x = c(1950, 1971, 1971, 1950), y = c(-10,-10, 130,130), fill = "gray60", angle = 45, density = 14)
dev.off()


png(filename = file.path(OM_fig_dir, "Catch_example3.png"), width = 7.5, height = 5, res = 250, units = "in")
plot(c(TMB_data$years), OM_sim$catches/1000, type = "l", lwd = 3, ylab = "Catch", xlab = "Years", ylim = c(0, 120))
lines(c(EM_data$years), EM_data$catches/1000,  lwd = 3, col = "red", lty = 2)
legend('topright', col = c("black", "gold", "purple", "red", "blue"), legend = c("OM", "EM1a", "EM1b", "EM2-4","Initial year"), lty = c(1,2, 2, 2,2), lwd = 3)
abline(v = start_year[3], lty = "dashed", col = "blue", lwd = 3)
polygon(x = c(1950, 1971, 1971, 1950), y = c(-10,-10, 130,130), fill = "gray60", angle = 45, density = 14)
dev.off()

png(filename = file.path(OM_fig_dir, "Catch_example4.png"), width = 7.5, height = 5, res = 250, units = "in")
plot(c(TMB_data$years), OM_sim$catches/1000, type = "l", lwd = 3, ylab = "Catch", xlab = "Years", ylim = c(0, 120))
lines(c(EM_data$years), EM_data$catches/1000,  lwd = 3, col = "red", lty = 2)
lines(c(EM5_data$years), EM5_data$catches/1000,  lwd = 3, col = "gold", lty = 2)
lines(c(EM6_data$years), EM6_data$catches/1000,  lwd = 3, col = "purple", lty = 2)
legend('topright', col = c("black", "gold", "purple", "red", "blue"), legend = c("OM", "EM1a", "EM1b", "EM2-4","Initial year"), lty = c(1,2, 2, 2,2), lwd = 3)
abline(v = start_year[3], lty = "dashed", col = "blue", lwd = 3)
polygon(x = c(1950, 1971, 1971, 1950), y = c(-10,-10, 130,130), fill = "gray60", angle = 45, density = 14)
dev.off()

# Catch scenario
png(filename = file.path(OM_fig_dir, "Catch_example5.png"), width = 7.5, height = 5, res = 250, units = "in")
plot(c(TMB_data$years), OM_sim$catches/1000, type = "l", lwd = 3, ylab = "Catch", xlab = "Years", ylim = c(0, 120))
abline(v = start_year, lty = "dashed", col = "blue", lwd = 2)
legend('topright', col = c("black", "gold", "purple", "red", "blue"), legend = c("OM", "EM1a", "EM1b", "EM2-4","Initial year"), lty = c(1,2, 2, 2,2), lwd = 3)
dev.off()


#########
## Plot across all OM & EMs
#########
get_multiple_Bzeros <- function(mle_lst) {
  mod_labs = names(mle_lst)
  B0_df = NULL
  for(i in 1:length(mod_labs)) {
    tmp_B0 = data.frame(B0 = mle_lst[[i]]$B0, label = mod_labs[i])
    B0_df = rbind(B0_df, tmp_B0)
  }
  return(B0_df)
}


get_multiple_scalar_vals <- function(mle_lst, component = "B0") {
  mod_labs = names(mle_lst)
  if(!component %in% names(mle_lst[[1]]))
    stop(paste0("Could not find ", component, " in mle_lst. Check spelling of 'component' parameter"))
  full_df = NULL
  for(i in 1:length(mod_labs)) {
    val = get(component, mle_lst[[i]])
    tmp_df = data.frame(value = val, label = mod_labs[i])
    full_df = rbind(full_df, tmp_df)
  }
  return(full_df)
}

get_multi_ssb <- function(mle_lst) {
  mod_labs = names(mle_lst)
  full_df = NULL
  for(i in 1:length(mod_labs)) {
    val = get("ssb", mle_lst[[i]])
    tmp_df = data.frame(SSB = val, years = c(min(mle_lst[[i]]$years) - 1, mle_lst[[i]]$years), depletion = val / mle_lst[[i]]$B0 * 100)
    full_df = rbind(full_df, tmp_df)
  }
  return(full_df)
}
get_multiple_vector_vals <- function(mle_lst, component = "B0", element = "last") {
  mod_labs = names(mle_lst)
  if(!component %in% names(mle_lst[[1]]))
    stop(paste0("Could not find ", component, " in mle_lst. Check spelling of 'component' parameter"))
  full_df = NULL
  for(i in 1:length(mod_labs)) {
    vec = get(component, mle_lst[[i]])
    val = NA
    if(element == "first")
      val = vec[1]
    if(element == "last")
      val = vec[length(vec)]
    tmp_df = data.frame(value = val, label = mod_labs[i])
    full_df = rbind(full_df, tmp_df)
  }
  return(full_df)
}
## Get a range of reference points
full_Bzero = NULL
for(start_year_ndx in 1:length(start_year)) {
  
  #EM1_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]], "B0")
  #EM1_B0$start_year = start_year[start_year_ndx]
  #EM1_B0$model = "EM1"
  #EM1_B0$RE = (EM1_B0$value - OM_rep$B0) / OM_rep$B0 * 100

  EM2_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]], "B0")
  EM2_B0$start_year = start_year[start_year_ndx]
  EM2_B0$model = "EM2"
  EM2_B0$RE = (EM2_B0$value - OM_rep$B0) / OM_rep$B0 * 100

  EM3_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]], "B0")
  EM3_B0$start_year = start_year[start_year_ndx]
  EM3_B0$model = "EM3"
  EM3_B0$RE = (EM3_B0$value - OM_rep$B0) / OM_rep$B0 * 100
  
  EM4_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]], "B0")
  EM4_B0$start_year = start_year[start_year_ndx]
  EM4_B0$model = "EM4"
  EM4_B0$RE = (EM4_B0$value - OM_rep$B0) / OM_rep$B0 * 100
  
  EM5_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]], "B0")
  EM5_B0$start_year = start_year[start_year_ndx]
  EM5_B0$model = "EM1a"
  EM5_B0$RE = (EM5_B0$value - OM_rep$B0) / OM_rep$B0 * 100
  
  
  EM6_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]], "B0")
  EM6_B0$start_year = start_year[start_year_ndx]
  EM6_B0$model = "EM1b"
  EM6_B0$RE = (EM6_B0$value - OM_rep$B0) / OM_rep$B0 * 100
  
  
  EM7_B0 = get_multiple_scalar_vals(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]], "B0")
  EM7_B0$start_year = start_year[start_year_ndx]
  EM7_B0$model = "EM1"
  EM7_B0$RE = (EM7_B0$value - OM_rep$B0) / OM_rep$B0 * 100
  
  full_Bzero = rbind(full_Bzero, EM2_B0, EM3_B0, EM4_B0, EM5_B0, EM6_B0, EM7_B0)
}
## Get a range of reference points
full_survey_q = NULL
OM_val = OM_rep$survey_Q
for(start_year_ndx in 1:length(start_year)) {
  #EM1_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]], "F_init")
  #EM1_survey_q$start_year = start_year[start_year_ndx]
  #EM1_survey_q$model = "EM1"
  
  EM2_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]], "survey_Q")
  EM2_survey_q$start_year = start_year[start_year_ndx]
  EM2_survey_q$model = "EM2"
  EM2_survey_q$RE = (EM2_survey_q$value - OM_val) / OM_val * 100
  
  EM3_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]], "survey_Q")
  EM3_survey_q$start_year = start_year[start_year_ndx]
  EM3_survey_q$model = "EM3"
  EM3_survey_q$RE = (EM3_survey_q$value - OM_val) / OM_val * 100
  
  EM4_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]], "survey_Q")
  EM4_survey_q$start_year = start_year[start_year_ndx]
  EM4_survey_q$model = "EM4"
  EM4_survey_q$RE = (EM4_survey_q$value - OM_val) / OM_val * 100
  
  EM5_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]], "survey_Q")
  EM5_survey_q$start_year = start_year[start_year_ndx]
  EM5_survey_q$model = "EM1a"
  EM5_survey_q$RE = (EM5_survey_q$value - OM_val) / OM_val * 100
  
  
  EM6_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]], "survey_Q")
  EM6_survey_q$start_year = start_year[start_year_ndx]
  EM6_survey_q$model = "EM1b"
  EM6_survey_q$RE = (EM6_survey_q$value - OM_val) / OM_val * 100
  
  EM7_survey_q = get_multiple_scalar_vals(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]], "survey_Q")
  EM7_survey_q$start_year = start_year[start_year_ndx]
  EM7_survey_q$model = "EM1"
  EM7_survey_q$RE = (EM7_survey_q$value - OM_val) / OM_val * 100
  
  full_survey_q = rbind(full_survey_q, EM2_survey_q, EM3_survey_q, EM4_survey_q, EM5_survey_q, EM6_survey_q, EM7_survey_q)
}
## Get a range of reference points
full_Finit = NULL
OM_val = OM_rep$F_init
for(start_year_ndx in 1:length(start_year)) {
  #EM1_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]], "F_init")
  #EM1_Finit$start_year = start_year[start_year_ndx]
  #EM1_Finit$model = "EM1"

  EM2_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]], "F_init")
  EM2_Finit$start_year = start_year[start_year_ndx]
  EM2_Finit$model = "EM2"

  EM3_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]], "F_init")
  EM3_Finit$start_year = start_year[start_year_ndx]
  EM3_Finit$model = "EM3"

  EM4_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]], "F_init")
  EM4_Finit$start_year = start_year[start_year_ndx]
  EM4_Finit$model = "EM4"
  
  EM5_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]], "F_init")
  EM5_Finit$start_year = start_year[start_year_ndx]
  EM5_Finit$model = "EM1a"
  
  
  EM6_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]], "F_init")
  EM6_Finit$start_year = start_year[start_year_ndx]
  EM6_Finit$model = "EM1b"
  
  EM7_Finit = get_multiple_scalar_vals(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]], "F_init")
  EM7_Finit$start_year = start_year[start_year_ndx]
  EM7_Finit$model = "EM1"
  
  full_Finit = rbind(full_Finit, EM2_Finit, EM3_Finit, EM4_Finit, EM5_Finit, EM6_Finit, EM7_Finit)
}
## Get a range of reference points
full_terminal_depletion = NULL
OM_val = OM_rep$depletion[length(OM_rep$depletion)]
for(start_year_ndx in 1:length(start_year)) {
  #EM1_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  #EM1_dep$start_year = start_year[start_year_ndx]
  #EM1_dep$model = "EM1"
  #EM1_dep$RE = (EM1_dep$value - OM_val) / OM_val * 100
  
  EM2_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  EM2_dep$start_year = start_year[start_year_ndx]
  EM2_dep$model = "EM2"
  EM2_dep$RE = (EM2_dep$value - OM_val) / OM_val * 100
  
  EM3_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  EM3_dep$start_year = start_year[start_year_ndx]
  EM3_dep$model = "EM3"
  EM3_dep$RE = (EM3_dep$value - OM_val) / OM_val * 100
  
  EM4_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  EM4_dep$start_year = start_year[start_year_ndx]
  EM4_dep$model = "EM4"
  EM4_dep$RE = (EM4_dep$value - OM_val) / OM_val * 100
  
  EM5_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  EM5_dep$start_year = start_year[start_year_ndx]
  EM5_dep$model = "EM1a"
  EM5_dep$RE = (EM5_dep$value - OM_val) / OM_val * 100
  
  EM6_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  EM6_dep$start_year = start_year[start_year_ndx]
  EM6_dep$model = "EM1b"
  EM6_dep$RE = (EM6_dep$value - OM_val) / OM_val * 100
  
  EM7_dep = get_multiple_vector_vals(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]], "depletion", "last")
  EM7_dep$start_year = start_year[start_year_ndx]
  EM7_dep$model = "EM1"
  EM7_dep$RE = (EM7_dep$value - OM_val) / OM_val * 100
  
  full_terminal_depletion = rbind(full_terminal_depletion, EM2_dep, EM3_dep, EM4_dep,EM5_dep, EM6_dep, EM7_dep)
}
## SPR_F40
full_SPR_F40 = NULL
OM_val = F_40$F_ref
for(start_year_ndx in 1:length(start_year)) {
  #EM1_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  #EM1_SPR_F40$start_year = start_year[start_year_ndx]
  #EM1_SPR_F40$model = "EM1"
  #EM1_SPR_F40$RE = (EM1_SPR_F40$value - OM_val) / OM_val * 100
  
  EM2_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  EM2_SPR_F40$start_year = start_year[start_year_ndx]
  EM2_SPR_F40$model = "EM2"
  EM2_SPR_F40$RE = (EM2_SPR_F40$value - OM_val) / OM_val * 100
  
  EM3_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  EM3_SPR_F40$start_year = start_year[start_year_ndx]
  EM3_SPR_F40$model = "EM3"
  EM3_SPR_F40$RE = (EM3_SPR_F40$value - OM_val) / OM_val * 100
  
  EM4_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  EM4_SPR_F40$start_year = start_year[start_year_ndx]
  EM4_SPR_F40$model = "EM4"
  EM4_SPR_F40$RE = (EM4_SPR_F40$value - OM_val) / OM_val * 100
  
  EM5_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  EM5_SPR_F40$start_year = start_year[start_year_ndx]
  EM5_SPR_F40$model = "EM1a"
  EM5_SPR_F40$RE = (EM5_SPR_F40$value - OM_val) / OM_val * 100
  
  
  EM6_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  EM6_SPR_F40$start_year = start_year[start_year_ndx]
  EM6_SPR_F40$model = "EM1b"
  EM6_SPR_F40$RE = (EM6_SPR_F40$value - OM_val) / OM_val * 100
  
  EM7_SPR_F40 = get_multiple_scalar_vals(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]], "SPR_F40")
  EM7_SPR_F40$start_year = start_year[start_year_ndx]
  EM7_SPR_F40$model = "EM1"
  EM7_SPR_F40$RE = (EM7_SPR_F40$value - OM_val) / OM_val * 100
  
  full_SPR_F40 = rbind(full_SPR_F40, EM2_SPR_F40, EM3_SPR_F40, EM4_SPR_F40, EM5_SPR_F40, EM6_SPR_F40, EM7_SPR_F40)
}
## F_0_1
full_F_0_1 = NULL
OM_val = F_0.1$F_0.1
for(start_year_ndx in 1:length(start_year)) {
  #EM1_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]], "F_0_1")
  #EM1_F_0_1$start_year = start_year[start_year_ndx]
  #EM1_F_0_1$model = "EM1"
  #EM1_F_0_1$RE = (EM1_F_0_1$value - OM_val) / OM_val * 100
  
  EM2_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]], "F_0_1")
  EM2_F_0_1$start_year = start_year[start_year_ndx]
  EM2_F_0_1$model = "EM2"
  EM2_F_0_1$RE = (EM2_F_0_1$value - OM_val) / OM_val * 100
  
  EM3_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]], "F_0_1")
  EM3_F_0_1$start_year = start_year[start_year_ndx]
  EM3_F_0_1$model = "EM3"
  EM3_F_0_1$RE = (EM3_F_0_1$value - OM_val) / OM_val * 100
  
  EM4_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]], "F_0_1")
  EM4_F_0_1$start_year = start_year[start_year_ndx]
  EM4_F_0_1$model = "EM4"
  EM4_F_0_1$RE = (EM4_F_0_1$value - OM_val) / OM_val * 100
  
  
  EM5_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]], "F_0_1")
  EM5_F_0_1$start_year = start_year[start_year_ndx]
  EM5_F_0_1$model = "EM1a"
  EM5_F_0_1$RE = (EM5_F_0_1$value - OM_val) / OM_val * 100
  
  EM6_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]], "F_0_1")
  EM6_F_0_1$start_year = start_year[start_year_ndx]
  EM6_F_0_1$model = "EM1b"
  EM6_F_0_1$RE = (EM6_F_0_1$value - OM_val) / OM_val * 100
  
  EM7_F_0_1 = get_multiple_scalar_vals(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]], "F_0_1")
  EM7_F_0_1$start_year = start_year[start_year_ndx]
  EM7_F_0_1$model = "EM1"
  EM7_F_0_1$RE = (EM7_F_0_1$value - OM_val) / OM_val * 100
  full_F_0_1 = rbind(full_F_0_1, EM2_F_0_1, EM3_F_0_1, EM4_F_0_1,EM5_F_0_1, EM6_F_0_1, EM7_F_0_1)
}

## ssbs
full_ssb = NULL
OM_val = OM_rep$ssb
years = c(min(TMB_data$years) - 1, TMB_data$years)
OM_ssb_df = data.frame(SSB = OM_val, years = years, depletion = OM_val / OM_rep$B0 * 100)
for(start_year_ndx in 1:length(start_year)) {
  #EM1_ssb = get_multi_ssb(mle_lst = mle_lst_EM1[[as.character(start_year[start_year_ndx])]])
  #EM1_ssb$start_year = start_year[start_year_ndx]
  #EM1_ssb$model = "EM1"
  #EM1_ssb$RE = (EM1_ssb$value - OM_val) / OM_val * 100
  
  EM2_ssb = get_multi_ssb(mle_lst = mle_lst_EM2[[as.character(start_year[start_year_ndx])]])
  EM2_ssb$start_year = start_year[start_year_ndx]
  EM2_ssb$model = "F-init"
  EM2_ssb = EM2_ssb %>% inner_join(OM_ssb_df, by = "years")
  EM2_ssb$RE = (EM2_ssb$SSB.x - EM2_ssb$SSB.y) / EM2_ssb$SSB.y * 100
  EM2_ssb$dep_RE = (EM2_ssb$depletion.x - EM2_ssb$depletion.y) / EM2_ssb$depletion.y * 100
  
  EM3_ssb = get_multi_ssb(mle_lst = mle_lst_EM3[[as.character(start_year[start_year_ndx])]])
  EM3_ssb$start_year = start_year[start_year_ndx]
  EM3_ssb$model = "EM3"
  EM3_ssb = EM3_ssb %>% inner_join(OM_ssb_df, by = "years")
  EM3_ssb$RE = (EM3_ssb$SSB.x - EM3_ssb$SSB.y) / EM3_ssb$SSB.y * 100
  EM3_ssb$dep_RE = (EM3_ssb$depletion.x - EM3_ssb$depletion.y) / EM3_ssb$depletion.y * 100
  
  EM4_ssb = get_multi_ssb(mle_lst = mle_lst_EM4[[as.character(start_year[start_year_ndx])]])
  EM4_ssb$start_year = start_year[start_year_ndx]
  EM4_ssb$model = "EM4"
  EM4_ssb = EM4_ssb %>% inner_join(OM_ssb_df, by = "years")
  EM4_ssb$RE = (EM4_ssb$SSB.x - EM4_ssb$SSB.y) / EM4_ssb$SSB.y * 100
  EM4_ssb$dep_RE = (EM4_ssb$depletion.x - EM4_ssb$depletion.y) / EM4_ssb$depletion.y * 100
  
  EM5_ssb = get_multi_ssb(mle_lst = mle_lst_EM5[[as.character(start_year[start_year_ndx])]])
  EM5_ssb$start_year = start_year[start_year_ndx]
  EM5_ssb$model = "overreported"
  EM5_ssb = EM5_ssb %>% inner_join(OM_ssb_df, by = "years")
  EM5_ssb$RE = (EM5_ssb$SSB.x - EM5_ssb$SSB.y) / EM5_ssb$SSB.y * 100
  EM5_ssb$dep_RE = (EM5_ssb$depletion.x - EM5_ssb$depletion.y) / EM5_ssb$depletion.y * 100
  
  EM6_ssb = get_multi_ssb(mle_lst = mle_lst_EM6[[as.character(start_year[start_year_ndx])]])
  EM6_ssb$start_year = start_year[start_year_ndx]
  EM6_ssb$model = "underreported"
  EM6_ssb = EM6_ssb %>% inner_join(OM_ssb_df, by = "years")
  EM6_ssb$RE = (EM6_ssb$SSB.x - EM6_ssb$SSB.y) / EM6_ssb$SSB.y * 100
  EM6_ssb$dep_RE = (EM6_ssb$depletion.x - EM6_ssb$depletion.y) / EM6_ssb$depletion.y * 100
  
  EM7_ssb = get_multi_ssb(mle_lst = mle_lst_EM7[[as.character(start_year[start_year_ndx])]])
  EM7_ssb$start_year = start_year[start_year_ndx]
  EM7_ssb$model = "Correct catch"
  EM7_ssb = EM7_ssb %>% inner_join(OM_ssb_df, by = "years")
  EM7_ssb$RE = (EM7_ssb$SSB.x - EM7_ssb$SSB.y) / EM7_ssb$SSB.y * 100
  EM7_ssb$dep_RE = (EM7_ssb$depletion.x - EM7_ssb$depletion.y) / EM7_ssb$depletion.y * 100
  
  full_ssb = rbind(full_ssb, EM2_ssb, EM3_ssb, EM4_ssb, EM5_ssb, EM6_ssb, EM7_ssb)
}


## boxplot of RE for survey Q's
ggplot(full_survey_q, aes(x = model, y = RE)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative error in Survey Q", x = "Model")
ggsave(filename = file.path(OM_fig_dir, "Survey_q_error.png"), width = 8, height = 6)


## boxplot of depletion
ggplot(full_terminal_depletion, aes(x = model, y = RE)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative error (%) for terminal depletion", x = "Model")
ggsave(filename = file.path(OM_fig_dir, "Depletion.png"), width = 8, height = 6)

## boxplot of depletion
ggplot(full_SPR_F40, aes(x = model, y = RE)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative error (%) for F 40%spr", x = "Model")
ggsave(filename = file.path(OM_fig_dir, "RE_F_40_percent.png"), width = 8, height = 6)

ggplot(full_F_0_1, aes(x = model, y = RE)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative error in F 0.1", x = "Model")
ggsave(filename = file.path(OM_fig_dir, "RE_F_0_1_percent.png"), width = 8, height = 6)

ggplot(full_F_0_1, aes(x = model, y = value)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  geom_hline(yintercept =  F_0.1$F_0.1, linetype = "dashed", col ="red") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "F 0.1", x = "Model")
ggsave(filename = file.path(OM_fig_dir, "F_0_1.png"), width = 8, height = 6)

## plot correlation between Q and B0 and Finit
scalar_cor_df = full_Bzero
scalar_cor_df = scalar_cor_df %>% inner_join( full_survey_q, by = c("label", "start_year", "model"))
scalar_cor_df = scalar_cor_df %>% inner_join( full_F_0_1, by = c("label", "start_year", "model"))

## B0 vs Q
ggplot(scalar_cor_df, aes(x = value.y, y = value.x/1000)) +
  geom_point() +
  facet_grid(start_year~model, scales = "free")+
  theme_bw() +
  labs(x = "Q", y = "B0") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, angle = 90))
## B0 vs F-init
ggplot(scalar_cor_df %>% filter(model %in% c("EM2", "EM4")), aes(x = value, y = value.x/1000)) +
  geom_point() +
  facet_grid(start_year~model, scales = "free")+
  theme_bw() +
  labs(x = "F-init", y = "B0") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

ggplot(full_Bzero, aes(x = model, y = RE)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative error (%) for B0", x = "Model") 
ggsave(filename = file.path(OM_fig_dir, "RE_B0.png"), width = 8, height = 6)

ggplot(full_Finit, aes(x = model, y = value)) +
  geom_boxplot() +
  facet_wrap(~start_year) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90)) +
  labs(y = "F-init parameter", x = "Model")

ggsave(filename = file.path(OM_fig_dir, "F_init.png"), width = 8, height = 6)

## Relative error in SSB and depletion
RE_ssb = get_df_quantiles(df = full_ssb, group_vars = c("years", "start_year", "model"), y_value = "RE", quants  = c(0,0.025, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.975, 0.9, 1))
RE_ssb_depletion = get_df_quantiles(df = full_ssb, group_vars = c("years", "start_year", "model"), y_value = "dep_RE", quants  = c(0,0.025, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.975, 0.9, 1))
# alternative plot
ggplot(RE_ssb %>% filter(model %in% c("Correct catch",  "overreported",  "underreported", "F-init")), aes(x = years)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(start_year~model) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  )
ggsave(filename = file.path(OM_fig_dir, "RE_SSB.png"), width =10, height = 8)

# alternative plot
ggplot(RE_ssb_depletion %>% filter(model %in% c("Correct catch",  "overreported",  "underreported", "F-init")), aes(x = years)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in depletion") +
  theme_bw() +
  facet_grid(start_year~model) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90),
    strip.text = element_text(size = 15)
  )
ggsave(filename = file.path(OM_fig_dir, "RE_Depletion.png"), width =10, height = 8)

## Do a profile on F-init
## see how its correlated with q and B0
