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

OM_label = "Medium_OM3"
output_data = file.path(DIR$data, OM_label)
if(!dir.exists(output_data))
  dir.create(output_data)

n_last_ycs_to_estimate = 1;
n_first_ycs_to_estimate = 1; ## you may want to shift this depending on when data starts recruits are observed
n_first_ycs_to_estimate_for_historic_models = 1; ## assume the first 10 year YCS are fixed at YCS

## to help stabilize the EM's

initial_level = 100
rebuild_level = 50
#fig_dir = file.path(DIR$fig, paste0("MediumBiology_", initial_level, "_", rebuild_level))
fig_dir = file.path(DIR$fig, "MediumBiology")
if(!dir.exists(fig_dir))
  dir.create(fig_dir)

this_bio = readRDS(file = file.path(DIR$data, "Medium_biology.RDS"))

## data period 
n_years_historic = 20
n_years_data = 40
n_years = n_years_historic + n_years_data
years = (2020 - n_years + 1):2020
full_years = (min(years) - 1):max(years)
## observation temporal frequency
survey_year_obs = years[(n_years_historic + 1):length(years)]
survey_ages = this_bio$ages
fishery_year_obs = years[(n_years_historic + 1):length(years)]
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
TMB_data$max_age_ndx_for_AFs = 30
TMB_data$survey_AF_obs = matrix(5, nrow = TMB_data$max_age_ndx_for_AFs, ncol = sum(TMB_data$survey_year_indicator))

TMB_data$fishery_year_indicator = array(as.integer(TMB_data$years %in% fishery_year_obs), dim = c(TMB_data$n_years, TMB_data$n_fisheries))
TMB_data$fishery_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, length(fishery_year_obs), TMB_data$n_fisheries))

TMB_data$catches = array(1000, dim = c(TMB_data$n_years, TMB_data$n_fisheries))# this will be overriden in the simulate() call
TMB_data$F_method = 0
TMB_data$F_iterations = 4
TMB_data$F_max = 3

TMB_data$catch_indicator = array(1, dim = c(TMB_data$n_years, TMB_data$n_fisheries))
TMB_data$ycs_estimated = c(rep(1, n_years))
TMB_data$standardise_ycs = 0;
TMB_data$ycs_bias_correction = rep(1, n_years)
## don't apply bias correction to the last 4 years because of lack of data
TMB_data$ycs_bias_correction[(n_years - 3):n_years] = 0

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
TMB_data$n_init_age_devs = max(TMB_data$ages)

TMB_data$rec_devs_sum_to_zero = 0
TMB_data$Q_r_for_sum_to_zero = Q_sum_to_zero_QR(length(TMB_data$years))
## iniital fishery_probs

fishery_probs = c(seq(from = 0.01, to = this_bio$M * 1.5, length = n_years_historic), rep(this_bio$M, n_years_data) )
plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)

## The same parameters as OM, to check for consistency
OM_pars = list(
  ln_R0 = log(this_bio$R0),
  ln_ycs_est =  rnorm(sum(TMB_data$ycs_estimated),  -0.5 * this_bio$sigma_r * this_bio$sigma_r, this_bio$sigma_r),
  ln_sigma_r = log( this_bio$sigma_r),
  #ln_extra_survey_cv = log(0.0001),
  ln_F_init = log(0.01),
  ln_init_age_devs = rep(0, TMB_data$n_init_age_devs),
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
na_map = fix_pars(par_list = OM_pars, pars_to_exclude = c("ln_catch_sd", "ln_sigma_r", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs"))
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

##
SR_code = ifelse(TMB_data$stockRecruitmentModelCode == 0, 0, 1)
F_30 = find_F_depletion(target_depletion = 30, fishery_sel = fishery_sel, M= M, waa= waa, paa = paa, ages = ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
F_35 = find_F_depletion(target_depletion = 35, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
F_40 = find_F_depletion(target_depletion = 40, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
##
F_75 = find_F_depletion(target_depletion = 75, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
F_50 = find_F_depletion(target_depletion = 50, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
F_25 = find_F_depletion(target_depletion = 25, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
F_20 = find_F_depletion(target_depletion = 20, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)
F_15 = find_F_depletion(target_depletion = 15, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness)

init_F_75 = find_F_depletion(target_depletion = 75, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness, n_years = n_years_historic)
init_F_50 = find_F_depletion(target_depletion = 50, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness, n_years = n_years_historic)
init_F_25 = find_F_depletion(target_depletion = 25, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 0.5),SRmodel = SR_code, SR_pars = TMB_data$steepness, n_years = n_years_historic)

###############
## Simulation
###############
## Create the short versio
data_years = years[(n_years_historic + 1):length(years)]
data_year_ndx = which(years %in% data_years)
data_full_years = (min(data_years) - 1):max(data_years)
#sim_survey_year_obs = sim_survey_year_obs[sim_survey_year_obs %in% data_years]
#sim_fishery_year_obs = sim_fishery_year_obs[sim_fishery_year_obs %in% data_years]
EM_short_data = TMB_data
EM_short_data$F_method = 1
EM_short_data$years = data_years
EM_short_data$n_years = length(data_years)
EM_short_data$ycs_bias_correction = rep(1, EM_short_data$n_years)
EM_short_data$ycs_bias_correction[(EM_short_data$n_years - 3):EM_short_data$n_years] = 0
short_survey_year_ndx = data_years %in% survey_year_obs 
short_survey_year_obs = data_years[short_survey_year_ndx]
short_fishery_year_ndx = data_years %in% fishery_year_obs
short_fishery_year_obs = data_years[short_fishery_year_ndx]
EM_short_data$survey_year_indicator = as.integer(short_survey_year_ndx)
EM_short_data$fishery_year_indicator = array(as.integer(short_fishery_year_ndx), dim = c(EM_short_data$n_years, EM_short_data$n_fisheries))
EM_short_data$survey_obs = rep(1, sum(short_survey_year_ndx))
EM_short_data$survey_cv = rep(0.15, sum(short_survey_year_ndx))
EM_short_data$survey_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, sum(short_survey_year_ndx)))
EM_short_data$fishery_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, sum(short_fishery_year_ndx), EM_short_data$n_fisheries))
EM_short_data$catches = array(EM_short_data$catches[EM_short_data$n_years,], dim = c(EM_short_data$n_years, EM_short_data$n_fisheries))
EM_short_data$ycs_estimated = rep(0, EM_short_data$n_years)
EM_short_data$ycs_estimated[n_first_ycs_to_estimate:(length(EM_short_data$ycs_estimated) - n_last_ycs_to_estimate )] = 1

EM_short_data$catch_indicator = array(1, dim = c(EM_short_data$n_years, EM_short_data$n_fisheries))
EM_short_data$stockMeanLength = EM_short_data$catchMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = EM_short_data$n_years, nrow = EM_short_data$n_ages)
EM_short_data$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = EM_short_data$n_years, nrow = EM_short_data$n_ages)
EM_short_data$propZ_ssb = EM_short_data$propZ_survey = rep(0.5, EM_short_data$n_years)
EM_short_pars = OM_pars
EM_short_pars$ln_ycs_est = rep(0, sum(EM_short_data$ycs_estimated))
EM_short_pars$ln_F = array(log(0.1), dim = c(EM_short_data$n_fisheries, EM_short_data$n_years))

## Data starting at the end of historical period and going until 2000
data_years_00 = years[(n_years_historic + 1):which(years == 2000)]
data_year_00_ndx = which(years %in% data_years_00)
data_full_00_years = (min(data_years_00) - 1):max(data_years_00)
EM_short_data_00 = TMB_data
EM_short_data_00$F_method = 1
EM_short_data_00$years = data_years_00
EM_short_data_00$n_years = length(data_years_00)
EM_short_data_00$ycs_bias_correction = rep(1, EM_short_data_00$n_years)
EM_short_data_00$ycs_bias_correction[(EM_short_data_00$n_years - 3):EM_short_data_00$n_years] = 0
survey_year_ndx_00 = data_years_00 %in% survey_year_obs
survey_year_obs_00 = data_years_00[survey_year_ndx_00]
fishery_year_ndx_00 = data_years_00 %in% fishery_year_obs 
fishery_year_obs_00 = data_years[fishery_year_ndx_00]
EM_short_data_00$survey_year_indicator = as.integer(survey_year_ndx_00)
EM_short_data_00$fishery_year_indicator = array(as.integer(fishery_year_ndx_00), dim = c(EM_short_data_00$n_years, EM_short_data_00$n_fisheries))
EM_short_data_00$survey_obs = rep(1, sum(EM_short_data_00$survey_year_indicator))
EM_short_data_00$survey_cv = rep(0.15, sum(EM_short_data_00$survey_year_indicator))
EM_short_data_00$survey_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, sum(EM_short_data_00$survey_year_indicator)))
EM_short_data_00$fishery_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, sum(fishery_year_ndx_00), EM_short_data_00$n_fisheries))
EM_short_data_00$catches = array(EM_short_data_00$catches[EM_short_data_00$n_years,], dim = c(EM_short_data_00$n_years, EM_short_data_00$n_fisheries))
EM_short_data_00$ycs_estimated = rep(0, EM_short_data_00$n_years)
EM_short_data_00$ycs_estimated[n_first_ycs_to_estimate:(length(EM_short_data_00$ycs_estimated) - n_last_ycs_to_estimate )] = 1
EM_short_data_00$catch_indicator = array(1, dim = c(EM_short_data_00$n_years, EM_short_data_00$n_fisheries))
EM_short_data_00$stockMeanLength = EM_short_data_00$catchMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = EM_short_data_00$n_years, nrow = EM_short_data_00$n_ages)
EM_short_data_00$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = EM_short_data_00$n_years, nrow = EM_short_data_00$n_ages)
EM_short_data_00$propZ_ssb = EM_short_data_00$propZ_survey = rep(0.5, EM_short_data_00$n_years)
EM_short_00_pars = OM_pars
EM_short_00_pars$ln_ycs_est = rep(0, sum(EM_short_data_00$ycs_estimated))
EM_short_00_pars$ln_F = array(log(0.1), dim = c(EM_short_data_00$n_fisheries, EM_short_data_00$n_years))
## Data starting at the beginning of the historical period and going until 2000
data_hist_years_00 = years[1:which(years == 2000)]
data_hist_year_00_ndx = which(years %in% data_hist_years_00)
data_full_hist_00_years = (min(data_hist_years_00) - 1):max(data_hist_years_00)
EM_hist_data_00 = TMB_data
EM_hist_data_00$F_method = 1
EM_hist_data_00$years = data_hist_years_00
EM_hist_data_00$n_years = length(data_hist_years_00)
EM_hist_data_00$ycs_bias_correction = rep(1, EM_hist_data_00$n_years)
EM_hist_data_00$ycs_bias_correction[(EM_hist_data_00$n_years - 3):EM_hist_data_00$n_years] = 0
survey_year_ndx_hist_00 = data_hist_years_00 %in% survey_year_obs
survey_year_obs_hist_00 = data_hist_years_00[survey_year_ndx_hist_00]
fishery_year_ndx_hist_00 = data_hist_years_00 %in%  fishery_year_obs 
fishery_year_obs_hist_00 = data_hist_years_00[fishery_year_ndx_hist_00]
EM_hist_data_00$survey_year_indicator = as.integer(survey_year_ndx_hist_00)
EM_hist_data_00$fishery_year_indicator = array(as.integer(fishery_year_ndx_hist_00), dim = c(EM_hist_data_00$n_years, EM_hist_data_00$n_fisheries))
EM_hist_data_00$survey_obs = rep(1, sum(EM_hist_data_00$survey_year_indicator))
EM_hist_data_00$survey_cv = rep(0.15, sum(EM_hist_data_00$survey_year_indicator))
EM_hist_data_00$survey_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, sum(EM_hist_data_00$survey_year_indicator)))
EM_hist_data_00$fishery_AF_obs = array(5, dim = c(TMB_data$max_age_ndx_for_AFs, sum(fishery_year_ndx_hist_00), EM_hist_data_00$n_fisheries))
EM_hist_data_00$catches = array(EM_hist_data_00$catches[EM_hist_data_00$n_years,], dim = c(EM_hist_data_00$n_years, EM_hist_data_00$n_fisheries))
EM_hist_data_00$ycs_estimated = rep(0, EM_hist_data_00$n_years)
EM_hist_data_00$ycs_estimated[n_first_ycs_to_estimate:(length(EM_hist_data_00$ycs_estimated) - n_last_ycs_to_estimate )] = 1
EM_hist_data_00$catch_indicator = array(1, dim = c(EM_hist_data_00$n_years, EM_hist_data_00$n_fisheries))
EM_hist_data_00$stockMeanLength = EM_hist_data_00$catchMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = EM_hist_data_00$n_years, nrow = EM_hist_data_00$n_ages)
EM_hist_data_00$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = EM_hist_data_00$n_years, nrow = EM_hist_data_00$n_ages)
EM_hist_data_00$propZ_ssb = EM_hist_data_00$propZ_survey = rep(0.5, EM_hist_data_00$n_years)
EM_hist_00_pars = OM_pars
EM_hist_00_pars$ln_ycs_est = rep(0, sum(EM_hist_data_00$ycs_estimated))
EM_hist_00_pars$ln_F = array(log(0.1), dim = c(EM_hist_data_00$n_fisheries, EM_hist_data_00$n_years))
## change YCS estimated for full_year models
EM1_data = TMB_data 
EM1_data$F_method = 1
EM1_data$ycs_estimated = rep(0, EM1_data$n_years)
EM1_data$ycs_estimated[n_first_ycs_to_estimate_for_historic_models:(length(EM1_data$ycs_estimated) - n_last_ycs_to_estimate )] = 1
EM_pars = OM_pars
EM_pars$ln_ycs_est = rep(0, sum(EM1_data$ycs_estimated))

n_first_ycs_to_estimate_for_historic_models

## test
OM_short<- MakeADFun(EM_short_data, EM_short_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_short$fn()
## Data only until 00
OM_short<- MakeADFun(EM_short_data_00, EM_short_00_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_short$fn()
## Data only until 00
OM_short<- MakeADFun(EM_hist_data_00, EM_hist_00_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_short$fn()

## OM simulations
## --> Initial
## ----> rebuild 
## ------> Simulations
n_sims = 100
SSB_df = recruit_df = depletion_df = NULL
inital_levels = c(25, 100)# c(25, 50, 75, 100)
rebuild_levels = c(20, 50) # c(20, 35, 50)
EM1a_data = EM1b_data = EM1_data
EM2_data = EM3_data = EM_short_data
EM2_data_00 = EM3_data_00 = EM_short_data_00
EM1_data_00 = EM1a_data_00 = EM1b_data_00 = EM_hist_data_00
EM1_data$estimate_F_init = EM1a_data$estimate_F_init = EM1b_data$estimate_F_init = 0
EM1_data$estimate_init_age_devs = EM1a_data$estimate_init_age_devs = EM1b_data$estimate_init_age_devs = 0
EM1_data_00$estimate_F_init = EM1a_data_00$estimate_F_init = EM1b_data_00$estimate_F_init = 0
EM1_data_00$estimate_init_age_devs = EM1a_data_00$estimate_init_age_devs = EM1b_data_00$estimate_init_age_devs = 0
EM2_data$estimate_F_init = EM2_data_00$estimate_F_init = EM3_data$estimate_F_init = EM3_data_00$estimate_F_init = 1
EM2_data$estimate_init_age_devs = EM2_data_00$estimate_init_age_devs = 0
EM3_data$estimate_init_age_devs = EM3_data_00$estimate_init_age_devs = 1
## sort out parameters
na_EM2_pars = fix_pars(EM_short_pars, pars_to_exclude = c("ln_sigma_r",  "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM3_pars = fix_pars(EM_short_pars, pars_to_exclude = c("ln_sigma_r",  "ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM2_pars_00 = fix_pars(EM_short_00_pars, pars_to_exclude = c("ln_sigma_r",  "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM3_pars_00 = fix_pars(EM_short_00_pars, pars_to_exclude = c("ln_sigma_r",  "ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM1_pars = fix_pars(EM_pars, pars_to_exclude =  c("ln_sigma_r",  "ln_F_init", "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM1a_pars = fix_pars(EM_pars, pars_to_exclude = c("ln_sigma_r",  "ln_F_init", "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM1b_pars = fix_pars(EM_pars, pars_to_exclude = c("ln_sigma_r",  "ln_F_init", "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM1_pars_00 = fix_pars(EM_hist_00_pars, pars_to_exclude =  c("ln_sigma_r",  "ln_F_init", "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM1a_pars_00 = fix_pars(EM_hist_00_pars, pars_to_exclude = c("ln_sigma_r",  "ln_F_init", "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_EM1b_pars_00 = fix_pars(EM_hist_00_pars, pars_to_exclude = c("ln_sigma_r",  "ln_F_init", "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
## test these pars
test <- MakeADFun(EM2_data, EM_short_pars, map = na_EM2_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM3_data, EM_short_pars, map = na_EM3_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM2_data_00, EM_short_00_pars, map = na_EM2_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM3_data_00, EM_short_00_pars, map = na_EM3_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM1_data, EM_pars, map = na_EM1_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM1a_data, EM_pars, map = na_EM1a_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM1b_data, EM_pars, map = na_EM1b_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM1_data_00, EM_hist_00_pars, map = na_EM1_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM1a_data_00, EM_hist_00_pars, map = na_EM1a_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T)
test <- MakeADFun(EM1b_data_00, EM_hist_00_pars, map = na_EM1b_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T)

EM1_convergence = EM1a_convergence = EM1b_convergence = EM1_00_convergence = EM1a_00_convergence = EM1b_00_convergence = 
  EM2_convergence = EM3_convergence = EM2_00_convergence = EM3_00_convergence = array(T, dim = c(length(inital_levels), length(rebuild_levels), n_sims))
mle_lst_EM1 = mle_lst_EM1a = mle_lst_EM1b = mle_lst_EM2 = mle_lst_EM3 = list()
mle_lst_EM1_00 = mle_lst_EM1a_00 = mle_lst_EM1b_00 = mle_lst_EM2_00 = mle_lst_EM3_00 = list()
OM_sim_lst = OM_rep_lst = list()
under_over_reporting_fraction = 0.25 ## 25%
EM_short_pars$ln_F_init = log(0.07)
EM_short_00_pars$ln_F_init = log(0.07)
#####
# Start Simulations
# 
#
#####
for(init_ndx in 1:length(inital_levels)) {
  cat("init ndx ", init_ndx, "\n")
  hist_F = NULL
  if(inital_levels[init_ndx] == 25) {
    hist_F = init_F_25$F_ref
  } else if (inital_levels[init_ndx] == 50) {
    hist_F = init_F_50$F_ref
  } else if (inital_levels[init_ndx] == 75) {
    hist_F = init_F_75$F_ref
  } else if (inital_levels[init_ndx] == 100) {
    hist_F = 0.0000001;
  }
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    cat("rebuild ndx ", rebuild_ndx, "\n")
    data_F = NULL;
    if(rebuild_levels[rebuild_ndx] == 20) {
      data_F = F_20$F_ref
    } else if (rebuild_levels[rebuild_ndx] == 35) {
      data_F = F_35$F_ref
    } else if (rebuild_levels[rebuild_ndx] == 50) {
      data_F = F_50$F_ref
    }
    OM_Fs = c(rep(hist_F, n_years_historic), rep(data_F, n_years_data))
    ## set OM F pars
    OM_pars$ln_F = array(log(OM_Fs), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
    ## Run simulations
    for(sim_iter in 1:n_sims) {
      if(sim_iter %% 10 == 0)
        cat("sim iter ", sim_iter, "\n")
      ## simulate YCS parameters
      OM_pars$ln_ycs_est = rnorm(sum(TMB_data$ycs_estimated), -0.5 * exp(OM_pars$ln_sigma_r) * exp(OM_pars$ln_sigma_r), exp(OM_pars$ln_sigma_r))
      ## Build OM and simulate parameters
      OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
      OM_sim <- OM_obj$simulate(complete = T)
      OM_sim_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = OM_sim
      ## Set catch for EM's
      EM1b_data$catches = EM1a_data$catches = EM1_data$catches = OM_sim$catches     ## self test
      EM1a_data$catches[1:n_years_historic, ] = matrix(OM_sim$catches[1:n_years_historic, ] * (1 + under_over_reporting_fraction), ncol = 1) ## over-reported
      EM1b_data$catches[1:n_years_historic, ] = matrix(OM_sim$catches[1:n_years_historic, ] * (1 - under_over_reporting_fraction), ncol = 1) ## under-reported
      EM2_data$catches = EM3_data$catches = matrix(OM_sim$catches[data_year_ndx, ], ncol = 1)
      EM2_data_00$catches = EM3_data_00$catches = matrix(OM_sim$catches[data_year_00_ndx, ], ncol = 1)
      EM1_data_00$catches = EM1a_data_00$catches =  EM1b_data_00$catches = matrix(OM_sim$catches[data_hist_year_00_ndx, ], ncol = 1)
      EM1a_data_00$catches[1:n_years_historic, ] = matrix(EM1a_data_00$catches[1:n_years_historic, ] * (1 + under_over_reporting_fraction), ncol = 1)
      EM1b_data_00$catches[1:n_years_historic, ] = matrix(EM1b_data_00$catches[1:n_years_historic, ] * (1 - under_over_reporting_fraction), ncol = 1)
      ## observations
      EM1_data$survey_obs = EM1a_data$survey_obs = EM1b_data$survey_obs = OM_sim$survey_obs
      EM1_data$survey_cv = EM1a_data$survey_cv = EM1b_data$survey_cv = OM_sim$survey_cv
      EM1_data$survey_AF_obs = EM1a_data$survey_AF_obs = EM1b_data$survey_AF_obs = OM_sim$survey_AF_obs
      EM1_data$fishery_AF_obs = EM1a_data$fishery_AF_obs = EM1b_data$fishery_AF_obs = OM_sim$fishery_AF_obs
      
      EM2_data$survey_obs = EM3_data$survey_obs = OM_sim$survey_obs[short_survey_year_ndx]
      EM2_data$survey_cv = EM3_data$survey_cv = OM_sim$survey_cv[short_survey_year_ndx]
      EM2_data$survey_AF_obs = EM3_data$survey_AF_obs = OM_sim$survey_AF_obs[,short_survey_year_ndx]
      EM2_data$fishery_AF_obs = EM3_data$fishery_AF_obs = array(OM_sim$fishery_AF_obs[,short_fishery_year_ndx,], dim = c(dim(OM_sim$fishery_AF_obs)[1],sum(short_fishery_year_ndx), OM_sim$n_fisheries))
      
      EM2_data_00$survey_obs = EM3_data_00$survey_obs = OM_sim$survey_obs[survey_year_obs %in% survey_year_obs_00]
      EM2_data_00$survey_cv = EM3_data_00$survey_cv = OM_sim$survey_cv[survey_year_obs %in% survey_year_obs_00]
      EM2_data_00$survey_AF_obs = EM3_data_00$survey_AF_obs = OM_sim$survey_AF_obs[,survey_year_obs %in% survey_year_obs_00]
      EM2_data_00$fishery_AF_obs = EM3_data_00$fishery_AF_obs = array(OM_sim$fishery_AF_obs[,fishery_year_obs %in% fishery_year_obs_00,], dim = c(dim(OM_sim$fishery_AF_obs)[1],sum(fishery_year_ndx_00), OM_sim$n_fisheries))
      
      EM1_data_00$survey_obs = EM1a_data_00$survey_obs = EM1b_data_00$survey_obs = OM_sim$survey_obs[survey_year_obs %in% survey_year_obs_hist_00]
      EM1_data_00$survey_cv = EM1a_data_00$survey_cv = EM1b_data_00$survey_cv = OM_sim$survey_cv[survey_year_obs %in% survey_year_obs_hist_00]
      EM1_data_00$survey_AF_obs = EM1a_data_00$survey_AF_obs = EM1b_data_00$survey_AF_obs = OM_sim$survey_AF_obs[,survey_year_obs %in% survey_year_obs_hist_00]
      EM1_data_00$fishery_AF_obs = EM1a_data_00$fishery_AF_obs = EM1b_data_00$fishery_AF_obs = array(OM_sim$fishery_AF_obs[,fishery_year_obs %in% fishery_year_obs_hist_00,], dim = c(dim(OM_sim$fishery_AF_obs)[1],sum(fishery_year_ndx_hist_00), OM_sim$n_fisheries))
      ## estimation
      EM2_obj <- MakeADFun(EM2_data, EM_short_pars, map = na_EM2_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      #SpatialSablefishAssessment::check_gradients(EM2_obj)
      EM3_obj <- MakeADFun(EM3_data, EM_short_pars, map = na_EM3_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      #SpatialSablefishAssessment::check_gradients(EM3_obj)
      OM_tmp_obj <- MakeADFun(OM_sim, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      OM_tmp_rep = OM_tmp_obj$report()
      OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = OM_tmp_rep
      EM2_00_obj <- MakeADFun(EM2_data_00, EM_short_00_pars, map = na_EM2_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM3_00_obj <- MakeADFun(EM3_data_00, EM_short_00_pars, map = na_EM3_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM1_obj <- MakeADFun(EM1_data, EM_pars, map = na_EM1_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM1a_obj <- MakeADFun(EM1a_data, EM_pars, map = na_EM1a_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM1b_obj <- MakeADFun(EM1b_data, EM_pars, map = na_EM1b_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM1_00_obj <- MakeADFun(EM1_data_00, EM_hist_00_pars, map = na_EM1_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM1a_00_obj <- MakeADFun(EM1a_data_00, EM_hist_00_pars, map = na_EM1a_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      EM1b_00_obj <- MakeADFun(EM1b_data_00, EM_hist_00_pars, map = na_EM1b_pars_00, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)

      ## Estimate
      mle_EM1 = nlminb(start = EM1_obj$par, objective = EM1_obj$fn, gradient  = EM1_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM1_obj$gr(mle_EM1$par))
                                 h = optimHess(mle_EM1$par, fn = EM1_obj$fn, gr = EM1_obj$gr)
                                 mle_EM1$par = mle_EM1$par - solve(h,g)
                                 mle_EM1$objective = EM1_obj$fn(mle_EM1$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM1, sim ", sim_iter, "\n")
        EM1_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM1_rep = EM1_obj$report(mle_EM1$par)
        mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM1_rep
      }

      ## EM1a
      mle_EM1a = nlminb(start = EM1a_obj$par, objective = EM1a_obj$fn, gradient  = EM1a_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM1a_obj$gr(mle_EM1a$par))
                                 h = optimHess(mle_EM1a$par, fn = EM1a_obj$fn, gr = EM1a_obj$gr)
                                 mle_EM1a$par = mle_EM1a$par - solve(h,g)
                                 mle_EM1a$objective = EM1a_obj$fn(mle_EM1a$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM1a, sim ", sim_iter, "\n")
        EM1a_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM1a_rep = EM1a_obj$report(mle_EM1a$par)
        mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM1a_rep
      }
      ## EM1b
      mle_EM1b = nlminb(start = EM1b_obj$par, objective = EM1b_obj$fn, gradient  = EM1b_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM1b_obj$gr(mle_EM1b$par))
                                 h = optimHess(mle_EM1b$par, fn = EM1b_obj$fn, gr = EM1b_obj$gr)
                                 mle_EM1b$par = mle_EM1b$par - solve(h,g)
                                 mle_EM1b$objective = EM1b_obj$fn(mle_EM1b$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM1b, sim ", sim_iter, "\n")
        EM1b_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM1b_rep = EM1b_obj$report(mle_EM1b$par)
        mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM1b_rep
      }
     
      ## EM2
      mle_EM2 = nlminb(start = EM2_obj$par, objective = EM2_obj$fn, gradient  = EM2_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM2_obj$gr(mle_EM2$par))
                                 h = optimHess(mle_EM2$par, fn = EM2_obj$fn, gr = EM2_obj$gr)
                                 mle_EM2$par = mle_EM2$par - solve(h,g)
                                 mle_EM2$objective = EM2_obj$fn(mle_EM2$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM2, sim ", sim_iter, "\n")
        EM2_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM2_rep = EM2_obj$report(mle_EM2$par)
        mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM2_rep
      }
      

      ## EM3
      mle_EM3 = nlminb(start = EM3_obj$par, objective = EM3_obj$fn, gradient  = EM3_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM3_obj$gr(mle_EM3$par))
                                 h = optimHess(mle_EM3$par, fn = EM3_obj$fn, gr = EM3_obj$gr)
                                 mle_EM3$par = mle_EM3$par - solve(h,g)
                                 mle_EM3$objective = EM3_obj$fn(mle_EM3$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM3, sim ", sim_iter, "\n")
        EM3_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM3_rep = EM3_obj$report(mle_EM3$par)
        mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM3_rep
      }
      
      ## Estimate
      mle_EM1_00 = nlminb(start = EM1_00_obj$par, objective = EM1_00_obj$fn, gradient  = EM1_00_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM1_00_obj$gr(mle_EM1_00$par))
                                 h = optimHess(mle_EM1_00$par, fn = EM1_00_obj$fn, gr = EM1_00_obj$gr)
                                 mle_EM1_00$par = mle_EM1_00$par - solve(h,g)
                                 mle_EM1_00$objective = EM1_00_obj$fn(mle_EM1_00$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM1_00, sim ", sim_iter, "\n")
        EM1_00_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM1_00_rep = EM1_00_obj$report(mle_EM1_00$par)
        mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM1_00_rep
      }
      
      ## EM1a_00
      mle_EM1a_00 = nlminb(start = EM1a_00_obj$par, objective = EM1a_00_obj$fn, gradient  = EM1a_00_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM1a_00_obj$gr(mle_EM1a_00$par))
                                 h = optimHess(mle_EM1a_00$par, fn = EM1a_00_obj$fn, gr = EM1a_00_obj$gr)
                                 mle_EM1a_00$par = mle_EM1a_00$par - solve(h,g)
                                 mle_EM1a_00$objective = EM1a_00_obj$fn(mle_EM1a_00$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM1a_00, sim ", sim_iter, "\n")
        EM1a_00_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM1a_00_rep = EM1a_00_obj$report(mle_EM1a_00$par)
        mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM1a_00_rep
      }
      ## EM1b_00
      mle_EM1b_00 = nlminb(start = EM1b_00_obj$par, objective = EM1b_00_obj$fn, gradient  = EM1b_00_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM1b_00_obj$gr(mle_EM1b_00$par))
                                 h = optimHess(mle_EM1b_00$par, fn = EM1b_00_obj$fn, gr = EM1b_00_obj$gr)
                                 mle_EM1b_00$par = mle_EM1b_00$par - solve(h,g)
                                 mle_EM1b_00$objective = EM1b_00_obj$fn(mle_EM1b_00$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM1b_00, sim ", sim_iter, "\n")
        EM1b_00_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM1b_00_rep = EM1b_00_obj$report(mle_EM1b_00$par)
        mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM1b_00_rep
      }
      
      ## EM2_00
      mle_EM2_00 = nlminb(start = EM2_00_obj$par, objective = EM2_00_obj$fn, gradient  = EM2_00_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM2_00_obj$gr(mle_EM2_00$par))
                                 h = optimHess(mle_EM2_00$par, fn = EM2_00_obj$fn, gr = EM2_00_obj$gr)
                                 mle_EM2_00$par = mle_EM2_00$par - solve(h,g)
                                 mle_EM2_00$objective = EM2_00_obj$fn(mle_EM2_00$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM2_00, sim ", sim_iter, "\n")
        EM2_00_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM2_00_rep = EM2_00_obj$report(mle_EM2_00$par)
        mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM2_00_rep
      }
      
      
      ## EM3_00
      mle_EM3_00 = nlminb(start = EM3_00_obj$par, objective = EM3_00_obj$fn, gradient  = EM3_00_obj$gr, control = list(iter.max = 10000, eval.max = 10000))
      try_improve = tryCatch(expr =
                               for(i in 1:2) {
                                 g = as.numeric(EM3_00_obj$gr(mle_EM3_00$par))
                                 h = optimHess(mle_EM3_00$par, fn = EM3_00_obj$fn, gr = EM3_00_obj$gr)
                                 mle_EM3_00$par = mle_EM3_00$par - solve(h,g)
                                 mle_EM3_00$objective = EM3_00_obj$fn(mle_EM3_00$par)
                               }
                             , error = function(e){e})
      
      if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
        cat("Failed simulation EM3_00, sim ", sim_iter, "\n")
        EM3_00_convergence[init_ndx, rebuild_ndx, sim_iter] = F
        next;
      } else {
        ## save the output
        EM3_00_rep = EM3_00_obj$report(mle_EM3_00$par)
        mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = EM3_00_rep
      }
      if(F) {
        ## plot SSBs
        plot(full_years , OM_tmp_rep$ssb, type= "l", ylim = c(0, 1400000))
        lines(c(min(full_years),data_hist_years_00), EM1_00_rep$ssb, col = "red")
        lines(full_years, EM1_rep$ssb, col = "blue", lty = 3)
        
        plot(years , OM_tmp_rep$catches, type= "l", ylim = c(0, 200000), lwd = 2)
        lines(data_hist_years_00, EM1_00_rep$catches, col = "red", lwd = 2, lty = 2)
        lines(years, EM1_rep$catches, col = "blue", lty = 3, lwd = 2)
        
        plot(survey_year_obs , OM_tmp_rep$survey_index_fitted, type= "l", ylim = c(0, 200000), lwd = 2)
        lines(survey_year_obs_00, EM1_00_rep$survey_index_fitted, col = "red", lwd = 2, lty = 2)
        lines(survey_year_obs, EM1_rep$survey_index_fitted, col = "blue", lty = 3, lwd = 2)
        
        plot(survey_year_obs , OM_tmp_rep$survey_obs, type= "l", ylim = c(0, 200000), lwd = 2)
        lines(survey_year_obs_00, EM1_00_rep$survey_obs, col = "red", lwd = 2, lty = 2)
        lines(survey_year_obs, EM1_rep$survey_obs, col = "blue", lty = 3, lwd = 2)
      }
    }
  }
}

saveRDS(object = mle_lst_EM1, file = file.path(output_data, "mle_lst_EM1.RDS"))
saveRDS(object = mle_lst_EM1a, file = file.path(output_data, "mle_lst_EM1a.RDS"))
saveRDS(object = mle_lst_EM1b, file = file.path(output_data, "mle_lst_EM1b.RDS"))
saveRDS(object = mle_lst_EM2, file = file.path(output_data, "mle_lst_EM2.RDS"))
saveRDS(object = mle_lst_EM3, file = file.path(output_data, "mle_lst_EM3.RDS"))
saveRDS(object = mle_lst_EM1_00, file = file.path(output_data, "mle_lst_EM1_00.RDS"))
saveRDS(object = mle_lst_EM1a_00, file = file.path(output_data, "mle_lst_EM1a_00.RDS"))
saveRDS(object = mle_lst_EM1b_00, file = file.path(output_data, "mle_lst_EM1b_00.RDS"))
saveRDS(object = mle_lst_EM2_00, file = file.path(output_data, "mle_lst_EM2_00.RDS"))
saveRDS(object = mle_lst_EM3_00, file = file.path(output_data, "mle_lst_EM3_00.RDS"))
#saveRDS(object = OM_pars, file = file.path(output_data, "OM_pars.RDS"))
#saveRDS(object = OM_obj, file = file.path(output_data, "OM_obj.RDS"))
saveRDS(object = OM_sim_lst, file = file.path(output_data, "OM_sim_lst.RDS"))
saveRDS(object = OM_rep_lst, file = file.path(output_data, "OM_rep_lst.RDS"))

## save convergence rates
saveRDS(object = EM1_convergence, file = file.path(output_data, "convergence_EM1.RDS"))
saveRDS(object = EM1a_convergence, file = file.path(output_data, "convergence_EM1a.RDS"))
saveRDS(object = EM1b_convergence, file = file.path(output_data, "convergence_EM1b.RDS"))
saveRDS(object = EM2_convergence, file = file.path(output_data, "convergence_EM2.RDS"))
saveRDS(object = EM3_convergence, file = file.path(output_data, "convergence_EM3.RDS"))
saveRDS(object = EM1_00_convergence, file = file.path(output_data, "convergence_EM1_00.RDS"))
saveRDS(object = EM1a_00_convergence, file = file.path(output_data, "convergence_EM1a_00.RDS"))
saveRDS(object = EM1b_00_convergence, file = file.path(output_data, "convergence_EM1b_00.RDS"))
saveRDS(object = EM2_00_convergence, file = file.path(output_data, "convergence_EM2_00.RDS"))
saveRDS(object = EM3_00_convergence, file = file.path(output_data, "convergence_EM3_00.RDS"))

