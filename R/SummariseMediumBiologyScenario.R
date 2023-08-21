#'
#' MediumBiologyScenario
#'
#'

source("AuxillaryFunctions.R")
library(dplyr)
library(tidyr)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
library(ggplot2)
library(reshape2)
library(TMB)
library(tidyr)
library(purrr)

## Pass the OM generated data to the TMB model
#sink(file = "compile_output.txt")
compile(file = file.path(DIR$tmb, "AgeStructuredModel_tmp.cpp"), flags = "-Wignored-attributes -O3")
#sink()
#dyn.unload(dynlib(file.path(DIR$tmb, "AgeStructuredModel_tmp")))
dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModel_tmp")))
#setwd(DIR$R)

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

OM_label = "Medium_OM5"
output_data = file.path(DIR$data, OM_label)
if(!dir.exists(output_data))
  dir.create(output_data)
output_fig_dir = file.path(DIR$fig, OM_label)
if(!dir.exists(output_fig_dir))
  dir.create(output_fig_dir)


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
TMB_data$survey_AF_obs = matrix(10, nrow = TMB_data$max_age_ndx_for_AFs, ncol = sum(TMB_data$survey_year_indicator))

TMB_data$fishery_year_indicator = array(as.integer(TMB_data$years %in% fishery_year_obs), dim = c(TMB_data$n_years, TMB_data$n_fisheries))
TMB_data$fishery_AF_obs = array(10, dim = c(TMB_data$max_age_ndx_for_AFs, length(fishery_year_obs), TMB_data$n_fisheries))

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
TMB_data$n_init_age_devs = 20

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
OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)

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
full_short_years = (min(data_years) - 1):max(data_years)
data_year_bool_ndx = (years %in% data_years)
short_year_ndx_for_dep = c(TRUE, data_year_bool_ndx)
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
short_survey_year_ndx = survey_year_obs %in% data_years
short_survey_year_obs = survey_year_obs[short_survey_year_ndx]
short_fishery_year_ndx = fishery_year_obs %in% data_years
short_fishery_year_obs = fishery_year_obs[short_fishery_year_ndx]
EM_short_data$survey_year_indicator = as.integer(short_survey_year_ndx)
EM_short_data$fishery_year_indicator = array(as.integer(short_fishery_year_ndx), dim = c(EM_short_data$n_years, EM_short_data$n_fisheries))
EM_short_data$survey_obs = rep(1, sum(short_survey_year_ndx))
EM_short_data$survey_cv = rep(0.15, sum(short_survey_year_ndx))
EM_short_data$survey_AF_obs = array(5, dim = c(EM_short_data$n_ages, sum(short_survey_year_ndx)))
EM_short_data$fishery_AF_obs = array(5, dim = c(EM_short_data$n_ages, sum(short_fishery_year_ndx), EM_short_data$n_fisheries))
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
data_year_00_bool_ndx = (years %in% data_years_00)
data_year_00_ndx_for_dep = c(TRUE, data_year_00_bool_ndx)
data_year_00_ndx = which(years %in% data_years_00)
data_full_00_years = (min(data_years_00) - 1):max(data_years_00)
EM_short_data_00 = TMB_data
EM_short_data_00$F_method = 1
EM_short_data_00$years = data_years_00
EM_short_data_00$n_years = length(data_years_00)
EM_short_data_00$ycs_bias_correction = rep(1, EM_short_data_00$n_years)
EM_short_data_00$ycs_bias_correction[(EM_short_data_00$n_years - 3):EM_short_data_00$n_years] = 0
survey_year_ndx_00 = survey_year_obs %in% data_years_00
survey_year_obs_00 = survey_year_obs[survey_year_ndx_00]
fishery_year_ndx_00 = fishery_year_obs %in% data_years_00
fishery_year_obs_00 = fishery_year_obs[fishery_year_ndx_00]
EM_short_data_00$survey_year_indicator = as.integer(survey_year_ndx_00)
EM_short_data_00$fishery_year_indicator = array(as.integer(fishery_year_ndx_00), dim = c(EM_short_data_00$n_years, EM_short_data_00$n_fisheries))
EM_short_data_00$survey_obs = rep(1, sum(EM_short_data_00$survey_year_indicator))
EM_short_data_00$survey_cv = rep(0.15, sum(EM_short_data_00$survey_year_indicator))
EM_short_data_00$survey_AF_obs = array(5, dim = c(EM_short_data_00$n_ages, sum(EM_short_data_00$survey_year_indicator)))
EM_short_data_00$fishery_AF_obs = array(5, dim = c(EM_short_data_00$n_ages, sum(fishery_year_ndx_00), EM_short_data_00$n_fisheries))
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
data_hist_year_00_bool_ndx = (years %in% data_hist_years_00)
data_hist_year_00_ndx_for_dep = c(TRUE, data_hist_year_00_bool_ndx)
data_hist_year_00_ndx = which(years %in% data_hist_years_00)
data_full_hist_00_years = (min(data_hist_years_00) - 1):max(data_hist_years_00)
EM_hist_data_00 = TMB_data
EM_hist_data_00$F_method = 1
EM_hist_data_00$years = data_hist_years_00
EM_hist_data_00$n_years = length(data_hist_years_00)
EM_hist_data_00$ycs_bias_correction = rep(1, EM_hist_data_00$n_years)
EM_hist_data_00$ycs_bias_correction[(EM_hist_data_00$n_years - 3):EM_hist_data_00$n_years] = 0
survey_year_ndx_hist_00 = survey_year_obs %in% data_hist_years_00
survey_year_obs_hist_00 = survey_year_obs[survey_year_ndx_hist_00]
fishery_year_ndx_hist_00 = fishery_year_obs %in% data_hist_years_00
fishery_year_obs_hist_00 = fishery_year_obs[fishery_year_ndx_hist_00]
EM_hist_data_00$survey_year_indicator = as.integer(survey_year_ndx_hist_00)
EM_hist_data_00$fishery_year_indicator = array(as.integer(fishery_year_ndx_hist_00), dim = c(EM_hist_data_00$n_years, EM_hist_data_00$n_fisheries))
EM_hist_data_00$survey_obs = rep(1, sum(EM_hist_data_00$survey_year_indicator))
EM_hist_data_00$survey_cv = rep(0.15, sum(EM_hist_data_00$survey_year_indicator))
EM_hist_data_00$survey_AF_obs = array(5, dim = c(EM_hist_data_00$n_ages, sum(EM_hist_data_00$survey_year_indicator)))
EM_hist_data_00$fishery_AF_obs = array(5, dim = c(EM_hist_data_00$n_ages, sum(fishery_year_ndx_hist_00), EM_hist_data_00$n_fisheries))
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
OM_short<- MakeADFun(EM_short_data, EM_short_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
OM_short$fn()
## Data only until 00
OM_short<- MakeADFun(EM_short_data_00, EM_short_00_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
OM_short$fn()
## Data only until 00
OM_short<- MakeADFun(EM_hist_data_00, EM_hist_00_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
OM_short$fn()

## OM simulations
## --> Initial
## ----> rebuild 
## ------> Simulations
n_sims = 100
SSB_df = recruit_df = depletion_df = NULL
inital_levels =c(25, 50, 75, 100)
rebuild_levels =  c(20, 35, 50)
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
test <- MakeADFun(EM2_data, EM_short_pars, map = na_EM2_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM3_data, EM_short_pars, map = na_EM3_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM2_data_00, EM_short_00_pars, map = na_EM2_pars_00, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM3_data_00, EM_short_00_pars, map = na_EM3_pars_00, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM1_data, EM_pars, map = na_EM1_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM1a_data, EM_pars, map = na_EM1a_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM1b_data, EM_pars, map = na_EM1b_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM1_data_00, EM_hist_00_pars, map = na_EM1_pars_00, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM1a_data_00, EM_hist_00_pars, map = na_EM1a_pars_00, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)
test <- MakeADFun(EM1b_data_00, EM_hist_00_pars, map = na_EM1b_pars_00, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T)

under_over_reporting_fraction = 0.25 ## 25%
EM_short_pars$ln_F_init = log(0.07)
EM_short_00_pars$ln_F_init = log(0.07)



############
## Load simulated data
############
mle_lst_EM1 = readRDS(file = file.path(output_data, "mle_lst_EM1.RDS"))
mle_lst_EM1a = readRDS(file = file.path(output_data, "mle_lst_EM1a.RDS"))
mle_lst_EM1b = readRDS(file = file.path(output_data, "mle_lst_EM1b.RDS"))
mle_lst_EM2 = readRDS(file = file.path(output_data, "mle_lst_EM2.RDS"))
mle_lst_EM3 = readRDS(file = file.path(output_data, "mle_lst_EM3.RDS"))
mle_lst_EM1_00 = readRDS(file = file.path(output_data, "mle_lst_EM1_00.RDS"))  
mle_lst_EM1a_00 = readRDS(file = file.path(output_data, "mle_lst_EM1a_00.RDS"))
mle_lst_EM1b_00 = readRDS(file = file.path(output_data, "mle_lst_EM1b_00.RDS"))  
mle_lst_EM2_00 = readRDS(file = file.path(output_data, "mle_lst_EM2_00.RDS"))  
mle_lst_EM3_00 = readRDS(file = file.path(output_data, "mle_lst_EM3_00.RDS"))
OM_sim_lst = readRDS(file = file.path(output_data, "OM_sim_lst.RDS"))
OM_rep_lst = readRDS(file = file.path(output_data, "OM_rep_lst.RDS"))

## read in Standard errors
se_lst_EM1 = readRDS(file = file.path(output_data, "se_lst_EM1.RDS"))
se_lst_EM1a = readRDS(file = file.path(output_data, "se_lst_EM1a.RDS"))
se_lst_EM1b = readRDS(file = file.path(output_data, "se_lst_EM1b.RDS"))
se_lst_EM2 = readRDS(file = file.path(output_data, "se_lst_EM2.RDS"))
se_lst_EM3 = readRDS(file = file.path(output_data, "se_lst_EM3.RDS"))

## summarise convergence
EM1_convergence = readRDS(file = file.path(output_data, "convergence_EM1.RDS"))
EM1a_convergence = readRDS(file = file.path(output_data, "convergence_EM1a.RDS"))
EM1b_convergence = readRDS(file = file.path(output_data, "convergence_EM1b.RDS"))
EM2_convergence = readRDS(file = file.path(output_data, "convergence_EM2.RDS"))
EM3_convergence = readRDS(file = file.path(output_data, "convergence_EM3.RDS"))
EM1_00_convergence = readRDS(file = file.path(output_data, "convergence_EM1_00.RDS"))
EM1a_00_convergence = readRDS(file = file.path(output_data, "convergence_EM1a_00.RDS"))
EM1b_00_convergence = readRDS(file = file.path(output_data, "convergence_EM1b_00.RDS"))
EM2_00_convergence = readRDS(file = file.path(output_data, "convergence_EM2_00.RDS"))
EM3_00_convergence = readRDS(file = file.path(output_data, "convergence_EM3_00.RDS"))

dimnames(EM1_convergence) = dimnames(EM1a_convergence) = dimnames(EM1b_convergence) = dimnames(EM2_convergence) = dimnames(EM3_convergence) = 
  dimnames(EM1_00_convergence) = dimnames(EM1a_00_convergence) = dimnames(EM1b_00_convergence) = dimnames(EM2_00_convergence) = dimnames(EM3_00_convergence) =list(inital_levels, rebuild_levels, 1:dim(EM1_convergence)[3])
n_sims = dim(EM1_convergence)[3]
EM1_con = melt(EM1_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM1_con$model = "EM1"
EM1a_con = melt(EM1a_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM1a_con$model = "EM1a"
EM1b_con = melt(EM1b_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM1b_con$model = "EM1b"
EM2_con = melt(EM2_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM2_con$model = "EM2"
EM3_con = melt(EM3_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM3_con$model = "EM3"
EM1_00_con = melt(EM1_00_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM1_00_con$model = "EM1_00"
EM1a_00_con = melt(EM1a_00_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM1a_00_con$model = "EM1a_00"
EM1b_00_con = melt(EM1b_00_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM1b_00_con$model = "EM1b_00"
EM2_00_con = melt(EM2_00_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM2_00_con$model = "EM2_00"
EM3_00_con = melt(EM3_00_convergence, varnames = c("init", "rebuild", "sim"), value.name = "convergence")
EM3_00_con$model = "EM3_00"
# merge them all together
full_con = rbind(EM1_con, EM1a_con,EM1b_con,EM2_con,EM3_con,
                 EM1_00_con, EM1a_00_con,EM1b_00_con,EM2_00_con,EM3_00_con)
## plot convergence table
convergence_tab = full_con %>% group_by(model, init, rebuild) %>% summarise(proportion_converged = sum(convergence) / n_sims * 100) %>%
  pivot_wider(id_cols = c("init","rebuild"), names_from = model, values_from = proportion_converged)
write.table(x = convergence_tab, file = file.path(output_fig_dir, "Convergence_table.txt"), row.names = F, col.names = T, quote = F)

#########
## Plot across all OM & EMs
#########

## Get a range of reference points
full_catch = NULL
OM_catch = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_catch_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = years)
    OM_obs_catch_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    OM_catch_df$observed = OM_catch_df$value
    OM_catch_df$init = paste0("init ", inital_levels[init_ndx])
    OM_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_catch_df$model = "OM"
    OM_catch = rbind(OM_catch, OM_catch_df)
    
    EM1_catch_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = years)
    EM1_obs_catch_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM1_catch_df$observed = EM1_obs_catch_df$value
    EM1_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_catch_df$model = "EM1"
    #EM1_catch_df$RE = (EM1_catch_df$values - OM_catch_df$values) / OM_catch_df$values * 100
    
    
    EM1a_catch_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = years)
    EM1a_obs_catch_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM1a_catch_df$observed = EM1a_obs_catch_df$value
    EM1a_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_catch_df$model = "EM1a"
    #EM1a_catch_df$RE = (EM1a_catch_df$values - OM_catch_df$values) / OM_catch_df$values * 100   
    
    EM1b_catch_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = years)
    EM1b_obs_catch_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM1b_catch_df$observed = EM1b_obs_catch_df$value
    EM1b_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_catch_df$model = "EM1b"
    #EM1b_catch_df$RE = (EM1b_catch_df$values - OM_catch_df$values) / OM_catch_df$values * 100
    
    EM2_catch_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_years)
    EM2_obs_catch_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM2_catch_df$observed = EM2_obs_catch_df$value
    EM2_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_catch_df$model = "EM2"
    #EM2_catch_df$RE = (EM2_catch_df$values - OM_catch_df$values[rep(data_year_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM3_catch_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_years)
    EM3_obs_catch_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM3_catch_df$observed = EM3_obs_catch_df$value
    EM3_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_catch_df$model = "EM3"
    #EM3_catch_df$RE = (EM3_catch_df$values - OM_catch_df$values[rep(data_year_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM1_00_catch_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_hist_years_00)
    EM1_00_obs_catch_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM1_00_catch_df$observed = EM1_00_obs_catch_df$value
    EM1_00_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_catch_df$model = "EM1_00"
    #EM1_00_catch_df$RE = (EM1_00_catch_df$values - OM_catch_df$values[rep(data_hist_year_00_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_hist_year_00_bool_ndx, n_sims)] * 100
    
    EM1a_00_catch_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_hist_years_00)
    EM1a_00_obs_catch_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM1a_00_catch_df$observed = EM1a_00_obs_catch_df$value
    EM1a_00_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_catch_df$model = "EM1a_00"
    #EM1a_00_catch_df$RE = (EM1a_00_catch_df$values - OM_catch_df$values[rep(data_hist_year_00_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_hist_year_00_bool_ndx, n_sims)] * 100
    
    EM1b_00_catch_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_hist_years_00)
    EM1b_00_obs_catch_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM1b_00_catch_df$observed = EM1b_00_obs_catch_df$value
    EM1b_00_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_catch_df$model = "EM1b_00"
    #EM1b_00_catch_df$RE = (EM1b_00_catch_df$values - OM_catch_df$values[rep(data_hist_year_00_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_hist_year_00_bool_ndx, n_sims)] * 100
    
    EM2_00_catch_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_years_00)
    EM2_00_obs_catch_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM2_00_catch_df$observed = EM2_00_obs_catch_df$value
    EM2_00_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_catch_df$model = "EM2_00"
    #EM2_00_catch_df$RE = (EM2_00_catch_df$values - OM_catch_df$values[rep(data_year_00_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_year_00_bool_ndx, n_sims)] * 100
    
    EM3_00_catch_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "pred_catches", element_labs = data_years_00)
    EM3_00_obs_catch_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "catches")
    EM3_00_catch_df$observed = EM3_00_obs_catch_df$value
    EM3_00_catch_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_catch_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_catch_df$model = "EM3_00"
    #EM3_00_catch_df$RE = (EM3_00_catch_df$values - OM_catch_df$values[rep(data_year_00_bool_ndx, n_sims)]) / OM_catch_df$values[rep(data_year_00_bool_ndx, n_sims)] * 100
    
    full_catch = rbind(full_catch, EM1_catch_df, EM1a_catch_df, EM1b_catch_df, EM2_catch_df, EM3_catch_df, EM1_00_catch_df, EM1a_00_catch_df, EM1b_00_catch_df, EM2_00_catch_df, EM3_00_catch_df)
  }
}

## match OM_catch to full catch
full_catch = full_catch %>% left_join(OM_catch, by = c("sim_iter", "init", "rebuild", "names"))
full_catch$init = factor(full_catch$init, levels = paste0("init ", inital_levels), ordered =T)
## calculate RE
full_catch = full_catch %>% group_by(sim_iter, init, rebuild, names, model.x) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
ggplot(full_catch %>% filter(model.x == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(init~rebuild)

ggplot(full_catch %>% filter(model.x == "EM1_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1_00") +
  facet_grid(init~rebuild)

ggplot(full_catch %>% filter(model.x == "EM1a"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1a") +
  facet_grid(init~rebuild)

ggplot(full_catch %>% filter(model.x == "EM1a_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1a") +
  facet_grid(init~rebuild)

ggplot(full_catch %>% filter(model.x == "EM1b"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1b") +
  facet_grid(init~rebuild)

ggplot(full_catch %>% filter(model.x == "EM2"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM2") +
  facet_grid(init~rebuild)

ggplot(full_catch %>% filter(model.x == "EM3"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM3") +
  facet_grid(init~rebuild)

## Get a range of reference points
full_Bzero = NULL
OM_Bzero = NULL

for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_b0_df = get_multiple_Bzeros(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    OM_b0_df$init = paste0("init ", inital_levels[init_ndx])
    OM_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_b0_df$model = "OM"
    OM_Bzero = rbind(OM_Bzero, OM_b0_df)
                     
    EM1_b0_df = get_multiple_Bzeros(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM1_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_b0_df$model = "EM1"
    #EM1_b0_df$RE = (EM1_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM1a_b0_df = get_multiple_Bzeros(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM1a_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_b0_df$model = "EM1a"
    #EM1a_b0_df$RE = (EM1a_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM1b_b0_df = get_multiple_Bzeros(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM1b_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_b0_df$model = "EM1b"
    #EM1b_b0_df$RE = (EM1b_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM2_b0_df = get_multiple_Bzeros(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM2_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_b0_df$model = "EM2"
    #EM2_b0_df$RE = (EM2_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM3_b0_df = get_multiple_Bzeros(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM3_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_b0_df$model = "EM3"
    #EM3_b0_df$RE = (EM3_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM1_00_b0_df = get_multiple_Bzeros(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM1_00_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_b0_df$model = "EM1_00"
    #EM1_00_b0_df$RE = (EM1_00_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM1a_00_b0_df = get_multiple_Bzeros(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM1a_00_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_b0_df$model = "EM1a_00"
    #EM1a_00_b0_df$RE = (EM1a_00_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM1b_00_b0_df = get_multiple_Bzeros(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM1b_00_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_b0_df$model = "EM1b_00"
    #EM1b_00_b0_df$RE = (EM1b_00_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM2_00_b0_df = get_multiple_Bzeros(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM2_00_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_b0_df$model = "EM2_00"
    #EM2_00_b0_df$RE = (EM2_00_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    
    EM3_00_b0_df = get_multiple_Bzeros(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]])
    EM3_00_b0_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_b0_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_b0_df$model = "EM3_00"
    #EM3_00_b0_df$RE = (EM3_00_b0_df$B0 - OM_b0_df$B0) / OM_b0_df$B0 * 100
    #boxplot(cbind(EM1_00_b0_df$B0,EM1_b0_df$B0))
    full_Bzero = rbind(full_Bzero, EM1_b0_df, EM1a_b0_df, EM1b_b0_df, EM2_b0_df, EM3_b0_df,EM1_00_b0_df, EM1a_00_b0_df, EM1b_00_b0_df, EM2_00_b0_df, EM3_00_b0_df)
  }
}

## match OM_catch to full catch
full_Bzero = full_Bzero %>% left_join(OM_Bzero, by = c("sim_iter", "init", "rebuild"))
full_Bzero$init = factor(full_Bzero$init, levels = paste0("init ", inital_levels), ordered =T)
## calculate RE
full_Bzero = full_Bzero %>% group_by(sim_iter, init, rebuild, model.x) %>% mutate(RE = (B0.x - B0.y)/ B0.y * 100)

ggplot(full_Bzero) +
  geom_boxplot(aes(x = model.x, y = RE, group = model.x)) +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 1.1) +
  theme_bw() +
  ggtitle("B0") +
  labs(x ="", y = "Relative error in B0") +
  ylim(-70,70) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(filename = file.path(output_fig_dir, "B0_RE.png"), width =10, height = 8)
## Nominal confidence coverage
OM_bzero = OM_report$B0
OM_CI_df = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    EM1_b0 = get_Bzero_coverage(mle_lst = se_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_bzero)
    EM1a_b0 = get_Bzero_coverage(mle_lst = se_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_bzero)
    EM1b_b0 = get_Bzero_coverage(mle_lst = se_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_bzero)
    EM2_b0 = get_Bzero_coverage(mle_lst = se_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_bzero)
    EM3_b0 = get_Bzero_coverage(mle_lst = se_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_bzero)
    tmp_df = data.frame(rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx]), init = paste0("init ", inital_levels[init_ndx]), model = c("EM1", "EM1a", "EM1b", "EM2", "EM3"), proportion = c(sum(EM1_b0)/length(EM1_b0), sum(EM1a_b0)/length(EM1a_b0), sum(EM1b_b0)/length(EM1b_b0), sum(EM2_b0)/length(EM2_b0), sum(EM3_b0)/length(EM3_b0)))
    OM_CI_df = rbind(OM_CI_df, tmp_df)
  }
}

Bzero_CI = OM_CI_df %>% pivot_wider(id_cols = c("rebuild", "init"), names_from = model, values_from = proportion)
write.table(x = Bzero_CI, file = file.path(output_fig_dir, "Bzero_coverage_table.txt"), row.names = F, col.names = T, quote = F)
Bzero_CI

## Get a range of reference points
full_survey_q = NULL
full_OM_q = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_q_df = get_multiple_scalar_vals(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    OM_q_df$init = paste0("init ", inital_levels[init_ndx])
    OM_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_q_df$model = "OM"
    full_OM_q = rbind(full_OM_q, OM_q_df)
    EM1_q_df = get_multiple_scalar_vals(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM1_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_q_df$model = "EM1"
    #EM1_q_df$RE = (EM1_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM1a_q_df = get_multiple_scalar_vals(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM1a_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_q_df$model = "EM1a"
    #EM1a_q_df$RE = (EM1a_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM1b_q_df = get_multiple_scalar_vals(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM1b_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_q_df$model = "EM1b"
    #EM1b_q_df$RE = (EM1b_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM2_q_df = get_multiple_scalar_vals(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM2_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_q_df$model = "EM2"
    #EM2_q_df$RE = (EM2_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM3_q_df = get_multiple_scalar_vals(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM3_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_q_df$model = "EM3"
    #EM3_q_df$RE = (EM3_q_df$value  - OM_q_df$value) / OM_q_df$value * 100
    
    EM1_00_q_df = get_multiple_scalar_vals(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM1_00_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_q_df$model = "EM1_00"
    #EM1_00_q_df$RE = (EM1_00_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM1a_00_q_df = get_multiple_scalar_vals(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM1a_00_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_q_df$model = "EM1a_00"
    #EM1a_00_q_df$RE = (EM1a_00_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM1b_00_q_df = get_multiple_scalar_vals(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM1b_00_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_q_df$model = "EM1b_00"
    #EM1b_00_q_df$RE = (EM1b_00_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM2_00_q_df = get_multiple_scalar_vals(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM2_00_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_q_df$model = "EM2_00"
    #EM2_00_q_df$RE = (EM2_00_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    EM3_00_q_df = get_multiple_scalar_vals(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_Q")
    EM3_00_q_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_q_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_q_df$model = "EM3_00"
    #EM3_00_q_df$RE = (EM3_00_q_df$value - OM_q_df$value) / OM_q_df$value * 100
    
    full_survey_q = rbind(full_survey_q, EM1_q_df, EM1a_q_df, EM1b_q_df, EM2_q_df, EM3_q_df, EM1_00_q_df, EM1a_00_q_df, EM1b_00_q_df, EM2_00_q_df, EM3_00_q_df)
  }
}
full_survey_q = full_survey_q %>% left_join(full_OM_q, by = c("sim_iter", "init", "rebuild"))
full_survey_q$init = factor(full_survey_q$init, levels = paste0("init ", inital_levels), ordered =T)
## calculate RE
full_survey_q = full_survey_q %>% group_by(sim_iter, init, rebuild, model.x) %>% mutate(RE = (value.x - value.y)/ value.y * 100)

ggplot(full_survey_q) +
  geom_boxplot(aes(x = model.x, y = RE, group = model.x)) +
  facet_grid(init~rebuild) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 1.1) +
  theme_bw() +
  ggtitle("") +
  labs(x ="", y = "Relative error in Q") +
  ylim(-30,30) +
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = file.path(output_fig_dir, "Q_RE.png"), width =10, height = 8)



## Get a range of reference points
full_terminal_depletion = NULL
OM_full_depletion = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_depletion_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = full_years)
    OM_depletion_df$values = OM_depletion_df$values * 100
    OM_depletion_df$init =  paste0("init ", inital_levels[init_ndx])
    OM_depletion_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_depletion_df$model = "OM"
    OM_full_depletion = rbind(OM_full_depletion, OM_depletion_df)
    ## EM1
    EM1_dep_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = full_years)
    EM1_dep_df$values = EM1_dep_df$values * 100
    EM1_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_dep_df$model = "EM1"
    #EM1_dep_df$RE = (EM1_dep_df$values - OM_depletion_df$values) / OM_depletion_df$values * 100

    ## EM1a
    EM1a_dep_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = full_years)
    EM1a_dep_df$values = EM1a_dep_df$values * 100
    EM1a_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_dep_df$model = "EM1a"
    #EM1a_dep_df$RE = (EM1a_dep_df$values - OM_depletion_df$values) / OM_depletion_df$values * 100
    
    ## EM1b
    EM1b_dep_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = full_years)
    EM1b_dep_df$values = EM1b_dep_df$values * 100
    EM1b_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_dep_df$model = "EM1b"
    #EM1b_dep_df$RE = (EM1b_dep_df$values - OM_depletion_df$values) / OM_depletion_df$values * 100
    
    ## EM2
    EM2_dep_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_years))
    EM2_dep_df$values = EM2_dep_df$values * 100
    EM2_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    #EM2_dep_df$RE = (EM2_dep_df$values - OM_depletion_df$values[rep(short_year_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(short_year_ndx_for_dep, n_sims)] * 100
    EM2_dep_df$model = "EM2"
    ## EM3
    EM3_dep_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_years))
    EM3_dep_df$values = EM3_dep_df$values * 100
    EM3_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_dep_df$model = "EM3"
    #EM3_dep_df$RE = (EM3_dep_df$values - OM_depletion_df$values[rep(short_year_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(short_year_ndx_for_dep, n_sims)] * 100
    
    ## EM1
    EM1_00_dep_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_hist_years_00))
    EM1_00_dep_df$values = EM1_00_dep_df$values * 100
    EM1_00_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_dep_df$model = "EM1_00"
    #EM1_00_dep_df$RE = (EM1_00_dep_df$values - OM_depletion_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)] * 100
    
    ## EM1a
    EM1a_00_dep_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_hist_years_00))
    EM1a_00_dep_df$values = EM1a_00_dep_df$values * 100
    EM1a_00_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_dep_df$model = "EM1a_00"
    #EM1a_00_dep_df$RE = (EM1a_00_dep_df$values - OM_depletion_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)] * 100
    
    ## EM1b
    EM1b_00_dep_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_hist_years_00))
    EM1b_00_dep_df$values = EM1b_00_dep_df$values * 100
    EM1b_00_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_dep_df$model = "EM1b_00"
    #EM1b_00_dep_df$RE = (EM1b_00_dep_df$values - OM_depletion_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)] * 100
    
    ## EM2
    EM2_00_dep_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_years_00))
    EM2_00_dep_df$values = EM2_00_dep_df$values * 100
    EM2_00_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    #EM2_00_dep_df$RE = (EM2_00_dep_df$values - OM_depletion_df$values[rep(data_year_00_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(data_year_00_ndx_for_dep, n_sims)] * 100
    EM2_00_dep_df$model = "EM2_00"
    ## EM3
    EM3_00_dep_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = c(min(full_years), data_years_00))
    EM3_00_dep_df$values = EM3_00_dep_df$values * 100
    EM3_00_dep_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_dep_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_dep_df$model = "EM3_00"
    #EM3_00_dep_df$RE = (EM3_00_dep_df$values - OM_depletion_df$values[rep(data_year_00_ndx_for_dep, n_sims)]) / OM_depletion_df$values[rep(data_year_00_ndx_for_dep, n_sims)] * 100
    
    full_terminal_depletion = rbind(full_terminal_depletion, EM1_dep_df, EM1a_dep_df, EM1b_dep_df, EM2_dep_df, EM3_dep_df,EM1_00_dep_df, EM1a_00_dep_df, EM1b_00_dep_df, EM2_00_dep_df, EM3_00_dep_df)
  }
}
## match OM_catch to full catch
full_terminal_depletion = full_terminal_depletion %>% left_join(OM_full_depletion, by = c("sim_iter", "init", "rebuild", "names"))
## calculate RE
full_terminal_depletion = full_terminal_depletion %>% group_by(sim_iter, init, rebuild, names, model.x) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
full_terminal_depletion$init = factor(full_terminal_depletion$init, levels = paste0("init ", inital_levels), ordered =T)
## plot the relative error
ggplot(full_terminal_depletion %>% filter(model.x == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(rebuild~init)

ggplot(full_terminal_depletion %>% filter(model.x == "EM1_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(rebuild~init)

ggplot(full_terminal_depletion %>% filter(model.x == "EM1a"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1a") +
  facet_grid(init~rebuild)

ggplot(full_terminal_depletion %>% filter(model.x == "EM1b"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1b") +
  facet_grid(init~rebuild)

ggplot(full_terminal_depletion %>% filter(model.x == "EM2"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM2") +
  facet_grid(init~rebuild)

ggplot(full_terminal_depletion %>% filter(model.x == "EM3"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM3") +
  facet_grid(init~rebuild)


ggplot(full_terminal_depletion %>% filter(names == 2020), aes(x = model.x, group = model.x)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("") +
  ylab("Relative error in depletion for terminal year") +
  facet_grid(rebuild~init)
ggsave(filename = file.path(output_fig_dir, "RE_terminal_depletion.png"), width =10, height = 8)

## calculate the MARE over simulations
depletion_MARE = full_terminal_depletion %>% group_by(names, init, rebuild, model.x) %>% summarise(MARE = median(abs(RE)))
depletion_MARE %>% group_by(init, rebuild, model.x) %>% summarise(min = min(MARE), max = max(MARE))
## plot MARE over entire time-series
ggplot(depletion_MARE, aes(x = model.x, group = model.x)) +
  geom_boxplot(aes(y = MARE)) + 
  theme_bw() +
  facet_grid(init~rebuild) +
  ylab("Median absolute relative error (depletion all years)") +
  theme(axis.text.x = element_text(angle = 90))
ggsave(filename = file.path(output_fig_dir, "Depletion_MARE.png"), width =10, height = 8)

## Relative error in depletion in terminal year 
ggplot(full_terminal_depletion %>% filter(names == max(full_years)), aes(x = model.x, group = model.x)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  labs(x = "EM", y = "Relative error in depletion for terminal year") +
  facet_grid(rebuild~init) +
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = file.path(output_fig_dir, "RE_terminal_year.png"), width =10, height = 8)

## Nominal confidence coverage
OM_CI_df = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_depletion = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "depletion", element_labs = full_years)
    OM_depletion = OM_depletion %>% filter(names == 2020)
    EM1_b0 = get_terminal_depletion_coverage(mle_lst = se_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_depletion)
    EM1a_b0 = get_terminal_depletion_coverage(mle_lst = se_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_depletion)
    EM1b_b0 = get_terminal_depletion_coverage(mle_lst = se_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_depletion)
    EM2_b0 = get_terminal_depletion_coverage(mle_lst = se_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_depletion)
    EM3_b0 = get_terminal_depletion_coverage(mle_lst = se_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], OM_val = OM_depletion)
    tmp_df = data.frame(rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx]), init = paste0("init ", inital_levels[init_ndx]), model = c("EM1", "EM1a", "EM1b", "EM2", "EM3"), proportion = c(sum(EM1_b0)/length(EM1_b0), sum(EM1a_b0)/length(EM1a_b0), sum(EM1b_b0)/length(EM1b_b0), sum(EM2_b0)/length(EM2_b0), sum(EM3_b0)/length(EM3_b0)), n = c(length(EM1_b0), length(EM1a_b0),length(EM1b_b0), length(EM2_b0), length(EM3_b0)))
    OM_CI_df = rbind(OM_CI_df, tmp_df)
  }
}

Terminal_Depletion_CI = OM_CI_df %>% pivot_wider(id_cols = c("rebuild", "init"), names_from = model, values_from = proportion)
write.table(x = Bzero_CI, file = file.path(output_fig_dir, "Bzero_coverage_table.txt"), row.names = F, col.names = T, quote = F)

## Get a range of reference points
full_ssb = NULL
OM_ssb = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_ssb_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = full_years)
    OM_ssb_df$values = OM_ssb_df$values * 100
    OM_ssb_df$init =  paste0("init ", inital_levels[init_ndx])
    OM_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_ssb_df$model = "OM"
    OM_ssb = rbind(OM_ssb, OM_ssb_df)
    ## EM1
    EM1_ssb_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = full_years)
    EM1_ssb_df$values = EM1_ssb_df$values * 100
    EM1_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_ssb_df$model = "EM1"
    #EM1_ssb_df$RE = (EM1_ssb_df$values - OM_ssb_df$values) / OM_ssb_df$values * 100
    
    ## EM1a
    EM1a_ssb_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = full_years)
    EM1a_ssb_df$values = EM1a_ssb_df$values * 100
    EM1a_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_ssb_df$model = "EM1a"
    #EM1a_ssb_df$RE = (EM1a_ssb_df$values - OM_ssb_df$values) / OM_ssb_df$values * 100
    
    ## EM1b
    EM1b_ssb_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = full_years)
    EM1b_ssb_df$values = EM1b_ssb_df$values * 100
    EM1b_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_ssb_df$model = "EM1b"
    #EM1b_ssb_df$RE = (EM1b_ssb_df$values - OM_ssb_df$values) / OM_ssb_df$values * 100
    
    ## EM2
    EM2_ssb_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_years))
    EM2_ssb_df$values = EM2_ssb_df$values * 100
    EM2_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    #EM2_ssb_df$RE = (EM2_ssb_df$values - OM_ssb_df$values[rep(short_year_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(short_year_ndx_for_dep, n_sims)] * 100
    EM2_ssb_df$model = "EM2"
    ## EM3
    EM3_ssb_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_years))
    EM3_ssb_df$values = EM3_ssb_df$values * 100
    EM3_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_ssb_df$model = "EM3"
    #EM3_ssb_df$RE = (EM3_ssb_df$values - OM_ssb_df$values[rep(short_year_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(short_year_ndx_for_dep, n_sims)] * 100
    
    ## EM1
    EM1_00_ssb_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_hist_years_00))
    EM1_00_ssb_df$values = EM1_00_ssb_df$values * 100
    EM1_00_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_ssb_df$model = "EM1_00"
    #EM1_00_ssb_df$RE = (EM1_00_ssb_df$values - OM_ssb_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)] * 100
    
    ## EM1a
    EM1a_00_ssb_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_hist_years_00))
    EM1a_00_ssb_df$values = EM1a_00_ssb_df$values * 100
    EM1a_00_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_ssb_df$model = "EM1a_00"
    #EM1a_00_ssb_df$RE = (EM1a_00_ssb_df$values - OM_ssb_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)] * 100
    
    ## EM1b
    EM1b_00_ssb_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_hist_years_00))
    EM1b_00_ssb_df$values = EM1b_00_ssb_df$values * 100
    EM1b_00_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_ssb_df$model = "EM1b_00"
    #EM1b_00_ssb_df$RE = (EM1b_00_ssb_df$values - OM_ssb_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(data_hist_year_00_ndx_for_dep, n_sims)] * 100
    
    ## EM2
    EM2_00_ssb_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_years_00))
    EM2_00_ssb_df$values = EM2_00_ssb_df$values * 100
    EM2_00_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    #EM2_00_ssb_df$RE = (EM2_00_ssb_df$values - OM_ssb_df$values[rep(data_year_00_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(data_year_00_ndx_for_dep, n_sims)] * 100
    EM2_00_ssb_df$model = "EM2_00"
    ## EM3
    EM3_00_ssb_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ssb", element_labs = c(min(full_years), data_years_00))
    EM3_00_ssb_df$values = EM3_00_ssb_df$values * 100
    EM3_00_ssb_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_ssb_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_ssb_df$model = "EM3_00"
    #EM3_00_ssb_df$RE = (EM3_00_ssb_df$values - OM_ssb_df$values[rep(data_year_00_ndx_for_dep, n_sims)]) / OM_ssb_df$values[rep(data_year_00_ndx_for_dep, n_sims)] * 100
    
    full_ssb = rbind(full_ssb, EM1_ssb_df, EM1a_ssb_df, EM1b_ssb_df, EM2_ssb_df, EM3_ssb_df,EM1_00_ssb_df, EM1a_00_ssb_df, EM1b_00_ssb_df, EM2_00_ssb_df, EM3_00_ssb_df)
  }
}
## match OM_catch to full catch
full_ssb = full_ssb %>% left_join(OM_ssb, by = c("sim_iter", "init", "rebuild", "names"))
## calculate RE
full_ssb = full_ssb %>% group_by(sim_iter, init, rebuild, names, model.x) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
full_ssb$init = factor(full_ssb$init, levels = paste0("init ", inital_levels), ordered =T)
## plot the relative error
## Relative error in depletion in terminal year 
ggplot(full_ssb %>% filter(names == max(full_years)), aes(x = model.x, group = model.x)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  labs(x = "EM", y = "Relative error in SSB for terminal year") +
  facet_grid(rebuild~init) +
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = file.path(output_fig_dir, "RE_terminal_year_SSB.png"), width =10, height = 8)

ggplot(full_ssb %>% filter(model.x == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(rebuild~init)

ggplot(full_ssb %>% filter(model.x == "EM1_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(rebuild~init)

ggplot(full_ssb %>% filter(model.x == "EM1a"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1a") +
  facet_grid(rebuild~init)

ggplot(full_ssb %>% filter(model.x == "EM1b"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1b") +
  facet_grid(rebuild~init)

ggplot(full_ssb %>% filter(model.x == "EM2"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM2") +
  facet_grid(rebuild~init)

ggplot(full_ssb %>% filter(model.x == "EM3"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM3") +
  facet_grid(rebuild~init)


## Relative error in SSB and depletion
RE_ssb = get_df_quantiles(df = full_ssb, group_vars = c("names", "rebuild", "init", "model.x"), y_value = "RE", quants  = c(0,0.025, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.975, 0.9, 1))
RE_ssb_depletion = get_df_quantiles(df = full_terminal_depletion, group_vars = c("names", "rebuild", "init", "model.x"), y_value = "RE", quants  = c(0,0.025, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.975, 0.9, 1))
# alternative plot
ggplot(RE_ssb %>% filter(model.x %in% c("EM1")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM1")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM1.png"), width =10, height = 8)

ggplot(RE_ssb %>% filter(model.x %in% c("EM1_00")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM1_00")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM1_00.png"), width =10, height = 8)

# alternative plot
ggplot(RE_ssb %>% filter(model.x %in% c("EM1a")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM1a")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM1a.png"), width =10, height = 8)

ggplot(RE_ssb %>% filter(model.x %in% c("EM1a_00")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM1a_00")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM1a_00.png"), width =10, height = 8)

ggplot(RE_ssb %>% filter(model.x %in% c("EM1b")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM1b")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM1b.png"), width =10, height = 8)

ggplot(RE_ssb %>% filter(model.x %in% c("EM1b_00")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM1b_00")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM1b_00.png"), width =10, height = 8)



ggplot(RE_ssb %>% filter(model.x %in% c("EM2")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM2")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM2.png"), width =10, height = 8)

ggplot(RE_ssb %>% filter(model.x %in% c("EM2_00")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM2_00")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM2_00.png"), width =10, height = 8)


ggplot(RE_ssb %>% filter(model.x %in% c("EM3")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM3")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM3.png"), width =10, height = 8)

ggplot(RE_ssb %>% filter(model.x %in% c("EM3_00")), aes(x = names)) +
  geom_line(aes(y=`50%`), linewidth = 1.1) +
  geom_line(aes(y=`2.5%`), colour='red', alpha=0.2) + geom_line(aes(y=`97.5%`), colour='red', alpha=0.2) +
  geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill='red', alpha=0.1) +
  geom_line(aes(y=`25%`), colour='gold', alpha=0.8) + geom_line(aes(y=`75%`), colour='gold', alpha=0.8) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill='gold', alpha=0.5) +
  geom_line(aes(y=`40%`), colour='blue', alpha=0.5) + geom_line(aes(y=`60%`), colour='blue', alpha=0.5) +
  geom_ribbon(aes(ymin=`40%`, ymax=`60%`), fill='blue', alpha=0.4) +
  labs(x = "Years", y = "Relative error (%) in SSB") +
  theme_bw() +
  facet_grid(rebuild~init) +
  geom_hline(yintercept = 0, linetype = "dashed", col ="red") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("EM3_00")
ggsave(filename = file.path(output_fig_dir, "RE_SSB_EM3_00.png"), width =10, height = 8)


## Recruitment
full_ycs = NULL
OM_ycs = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_ycs_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = years)
    OM_ycs_df$init =  paste0("init ", inital_levels[init_ndx])
    OM_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_ycs_df$model = "OM"
    OM_ycs = rbind(OM_ycs, OM_ycs_df)
    ## EM1
    EM1_ycs_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = years)
    EM1_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_ycs_df$model = "EM1"
    #EM1_ycs_df$RE = (EM1_ycs_df$values - OM_ycs_df$values) / OM_ycs_df$values * 100
    
    ## EM1a
    EM1a_ycs_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = years)
    EM1a_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_ycs_df$model = "EM1a"
    #EM1a_ycs_df$RE = (EM1a_ycs_df$values - OM_ycs_df$values) / OM_ycs_df$values * 100
    
    ## EM1b
    EM1b_ycs_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = years)
    EM1b_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_ycs_df$model = "EM1b"
    #EM1b_ycs_df$RE = (EM1b_ycs_df$values - OM_ycs_df$values) / OM_ycs_df$values * 100
    
    ## EM2
    EM2_ycs_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c(data_years))
    EM2_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    #EM2_ycs_df$RE = (EM2_ycs_df$values - OM_ycs_df$values[rep(data_year_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_year_bool_ndx, n_sims)] * 100
    EM2_ycs_df$model = "EM2"
    ## EM3
    EM3_ycs_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c(data_years))
    EM3_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_ycs_df$model = "EM3"
    #EM3_ycs_df$RE = (EM3_ycs_df$values - OM_ycs_df$values[rep(data_year_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_year_bool_ndx, n_sims)] * 100
    
    ## EM1
    EM1_00_ycs_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c(data_hist_years_00))
    EM1_00_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_ycs_df$model = "EM1_00"
    #EM1_00_ycs_df$RE = (EM1_00_ycs_df$values - OM_ycs_df$values[rep(data_hist_year_00_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_hist_year_00_bool_ndx, n_sims)] * 100
    
    ## EM1a
    EM1a_00_ycs_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c(data_hist_years_00))
    EM1a_00_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_ycs_df$model = "EM1a_00"
    #EM1a_00_ycs_df$RE = (EM1a_00_ycs_df$values - OM_ycs_df$values[rep(data_hist_year_00_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_hist_year_00_bool_ndx, n_sims)] * 100
    
    ## EM1b
    EM1b_00_ycs_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c(data_hist_years_00))
    EM1b_00_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_ycs_df$model = "EM1b_00"
    #EM1b_00_ycs_df$RE = (EM1b_00_ycs_df$values - OM_ycs_df$values[rep(data_hist_year_00_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_hist_year_00_bool_ndx, n_sims)] * 100
    
    ## EM2
    EM2_00_ycs_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c(data_years_00))
    EM2_00_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    #EM2_00_ycs_df$RE = (EM2_00_ycs_df$values - OM_ycs_df$values[rep(data_year_00_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_year_00_bool_ndx, n_sims)] * 100
    EM2_00_ycs_df$model = "EM2_00"
    ## EM3
    EM3_00_ycs_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "ycs", element_labs = c( data_years_00))
    EM3_00_ycs_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_ycs_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_ycs_df$model = "EM3_00"
    #EM3_00_ycs_df$RE = (EM3_00_ycs_df$values - OM_ycs_df$values[rep(data_year_00_bool_ndx, n_sims)]) / OM_ycs_df$values[rep(data_year_00_bool_ndx, n_sims)] * 100
    
    full_ycs = rbind(full_ycs, EM1_ycs_df, EM1a_ycs_df, EM1b_ycs_df, EM2_ycs_df, EM3_ycs_df,EM1_00_ycs_df, EM1a_00_ycs_df, EM1b_00_ycs_df, EM2_00_ycs_df, EM3_00_ycs_df)
  }
}
## match OM_catch to full catch
full_ycs = full_ycs %>% left_join(OM_ycs, by = c("sim_iter", "init", "rebuild", "names"))
## calculate RE
full_ycs = full_ycs %>% group_by(sim_iter, init, rebuild, names, model.x) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
full_ycs$init = factor(full_ycs$init, levels = paste0("init ", inital_levels), ordered =T)
## plot the relative error

ggplot(OM_ycs, aes(x = names, group = names,y = values)) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("OM") +
  stat_summary(
    fun = "mean",
    geom = "point",
    size = 2,
    alpha = 0.6,
    shape = 16,
    color = "steelblue") + 
  geom_hline(yintercept = 1, col = "red", linewidth = 1.1, linetype = "dashed") +
  facet_grid(init~rebuild) +
  ylim(0,6)

OM_ycs %>% group_by(names) %>% summarise(mean(values))

ggplot(full_ycs %>% filter(model.x == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(init~rebuild) +
  ylim(-100,100)

ggplot(full_ycs %>% filter(model.x == "EM1_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(init~rebuild) +
  ylim(-100,100)

## look at initial numbers at age
full_nage = NULL
OM_nage = NULL
init_short_year = data_years_00[1]
short_ndx = which(full_years %in% init_short_year)
ndx_00 = which(full_years %in% data_hist_years_00[length(data_hist_years_00)])
year_00 = data_hist_years_00[length(data_hist_years_00)]
short_00_ndx = which(data_years_00 %in% data_hist_years_00[length(data_hist_years_00)])
short_00_year = data_hist_years_00[length(data_hist_years_00)]
hist_shrt_ndx = which(data_hist_years_00 %in% data_years_00[1])
hist_shrt_yr = data_years_00[1]

for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_nage_df = get_numbers_at_age(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,short_ndx,length(full_years), ndx_00), year_label = c(full_years[1], init_short_year, full_years[length(full_years)],year_00))
    OM_nage_df$init =  paste0("init ", inital_levels[init_ndx])
    OM_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_nage_df$model = "OM"
    OM_nage = rbind(OM_nage, OM_nage_df)
    ## EM1
    EM1_nage_df = get_numbers_at_age(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,short_ndx,length(full_years), ndx_00), year_label = c(full_years[1], init_short_year, full_years[length(full_years)],year_00))
    EM1_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_nage_df$model = "EM1"
    ## EM1a
    EM1a_nage_df = get_numbers_at_age(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,short_ndx,length(full_years), ndx_00), year_label = c(full_years[1], init_short_year, full_years[length(full_years)],year_00))
    EM1a_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_nage_df$model = "EM1a"

    ## EM1b
    EM1b_nage_df = get_numbers_at_age(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,short_ndx,length(full_years), ndx_00), year_label = c(full_years[1], init_short_year, full_years[length(full_years)],year_00))
    EM1b_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_nage_df$model = "EM1b"

    ## EM2
    EM2_nage_df = get_numbers_at_age(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,short_00_ndx, length(data_years_00)), year_label = c(data_years_00[1], short_00_year, year_00))
    EM2_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_nage_df$model = "EM2"
    ## EM3
    EM3_nage_df = get_numbers_at_age(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,short_00_ndx, length(data_years_00)), year_label = c(data_years_00[1], short_00_year, year_00))
    EM3_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_nage_df$model = "EM3"

    ## EM1
    EM1_00_nage_df = get_numbers_at_age(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,length(data_hist_years_00), hist_shrt_ndx), year_label = c(data_hist_years_00[1], data_hist_years_00[length(data_hist_years_00)],hist_shrt_yr))
    EM1_00_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_nage_df$model = "EM1_00"

    ## EM1a
    EM1a_00_nage_df = get_numbers_at_age(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,length(data_hist_years_00), hist_shrt_ndx), year_label = c(data_hist_years_00[1], data_hist_years_00[length(data_hist_years_00)],hist_shrt_yr))
    EM1a_00_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_nage_df$model = "EM1a_00"

    ## EM1b
    EM1b_00_nage_df = get_numbers_at_age(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,length(data_hist_years_00), hist_shrt_ndx), year_label = c(data_hist_years_00[1], data_hist_years_00[length(data_hist_years_00)],hist_shrt_yr))
    EM1b_00_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_nage_df$model = "EM1b_00"

    ## EM2
    EM2_00_nage_df = get_numbers_at_age(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,length(data_years_00)), year_label = c(data_years_00[1], data_years_00[length(data_years_00)]))
    EM2_00_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_nage_df$model = "EM2_00"
    ## EM3
    EM3_00_nage_df = get_numbers_at_age(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], year_ndx = c(1,length(data_years_00)), year_label = c(data_years_00[1], data_years_00[length(data_years_00)]))
    EM3_00_nage_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_nage_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_nage_df$model = "EM3_00"

    full_nage = rbind(full_nage, EM1_nage_df, EM1a_nage_df, EM1b_nage_df, EM2_nage_df, EM3_nage_df,EM1_00_nage_df, EM1a_00_nage_df, EM1b_00_nage_df, EM2_00_nage_df, EM3_00_nage_df)
  }
}
## match OM_catch to full catch
full_nage = full_nage %>% left_join(OM_nage, by = c("sim_iter", "init", "rebuild", "year","Age"))
## calculate RE
full_nage = full_nage %>% group_by(sim_iter, init, rebuild, year,Age, model.x) %>% mutate(RE = (nage.x - nage.y)/ nage.y * 100)
full_nage$init = factor(full_nage$init, levels = paste0("init ", inital_levels), ordered =T)
## plot the relative error

## plot the relative error
ggplot(full_nage %>% filter(model.x == "EM3", year == year_00, sim_iter == 1), aes(x = Age)) +
  geom_line(aes(y = nage.x, col = "EM3"), linewidth = 1.1) +
  geom_line(aes(y = nage.y, col = "OM"), linewidth = 1.1) +
  theme_bw() +
  ggtitle("EM3") +
  facet_grid(rebuild ~ init) 

ggplot(full_nage %>% filter(model.x == "EM2", year == data_years_00[1]), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM2") +
  ylim(-50, 50) +
  facet_grid(init~rebuild) 

ggplot(full_nage %>% filter(model.x == "EM3", year == data_years_00[1]), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM3") +
  ylim(-50, 50) +
  facet_grid(init~rebuild) 

ggplot(full_nage %>% filter(model.x %in% c("EM3", "OM"), year == data_years[1], sim_iter == 1), aes(x = Age, group = Age, y = RE)) +
  geom_line(linewidth = 1.1, aes(col = model.x)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM3") +
  #ylim(-50, 50) +
  facet_grid(init~rebuild) 

## Get a range of reference points
full_finit = NULL
full_OM_finit = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
  
    EM2_finit_df = get_multiple_scalar_vals(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_finit")
    EM2_finit_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_finit_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_finit_df$model = "EM2"
    #EM2_finit_df$RE = (EM2_finit_df$value - OM_finit_df$value) / OM_finit_df$value * 100
    
    EM3_finit_df = get_multiple_scalar_vals(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_finit")
    EM3_finit_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_finit_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_finit_df$model = "EM3"
    #EM3_finit_df$RE = (EM3_finit_df$value  - OM_finit_df$value) / OM_finit_df$value * 100
    
    EM2_00_finit_df = get_multiple_scalar_vals(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_finit")
    EM2_00_finit_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_finit_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_finit_df$model = "EM2_00"
    #EM2_00_finit_df$RE = (EM2_00_finit_df$value - OM_finit_df$value) / OM_finit_df$value * 100
    
    EM3_00_finit_df = get_multiple_scalar_vals(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_finit")
    EM3_00_finit_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_finit_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_finit_df$model = "EM3_00"
    #EM3_00_finit_df$RE = (EM3_00_finit_df$value - OM_finit_df$value) / OM_finit_df$value * 100
    
    full_finit = rbind(full_finit, EM2_finit_df, EM3_finit_df, EM2_00_finit_df, EM3_00_finit_df)
  }
}
full_survey_q = full_survey_q %>% left_join(full_OM_q, by = c("sim_iter", "init", "rebuild"))
full_survey_q$init = factor(full_survey_q$init, levels = paste0("init ", inital_levels), ordered =T)
## calculate RE
full_survey_q = full_survey_q %>% group_by(sim_iter, init, rebuild, model.x) %>% mutate(RE = (value.x - value.y)/ value.y * 100)

## look at fits to observations
## Get a range of reference points
full_index = NULL
OM_index_fit = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_index_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs)
    OM_obs_index_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    OM_index_df$observed = OM_index_df$value
    OM_index_df$init = paste0("init ", inital_levels[init_ndx])
    OM_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_index_df$model = "OM"
    OM_index_fit = rbind(OM_index_fit, OM_index_df)
    
    EM1_index_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs)
    EM1_obs_index_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM1_index_df$observed = EM1_obs_index_df$value
    EM1_sd_index_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM1_index_df$normalised_resid = (log(EM1_index_df$observed / EM1_index_df$values) + EM1_sd_index_df$values^2) /  EM1_sd_index_df$values
    EM1_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_index_df$model = "EM1"
    #EM1_index_df$RE = (EM1_index_df$values - OM_index_df$values) / OM_index_df$values * 100
    
    
    EM1a_index_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs)
    EM1a_obs_index_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM1a_index_df$observed = EM1a_obs_index_df$value
    EM1a_sd_index_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM1a_index_df$normalised_resid = (log(EM1a_index_df$observed / EM1a_index_df$values) + EM1a_sd_index_df$values^2) /  EM1a_sd_index_df$values
    EM1a_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_index_df$model = "EM1a"
    #EM1a_index_df$RE = (EM1a_index_df$values - OM_index_df$values) / OM_index_df$values * 100   
    
    EM1b_index_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs)
    EM1b_obs_index_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM1b_index_df$observed = EM1b_obs_index_df$value
    EM1b_sd_index_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM1b_index_df$normalised_resid = (log(EM1b_index_df$observed / EM1b_index_df$values) + EM1b_sd_index_df$values^2) /  EM1b_sd_index_df$values
    EM1b_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_index_df$model = "EM1b"
    #EM1b_index_df$RE = (EM1b_index_df$values - OM_index_df$values) / OM_index_df$values * 100
    
    EM2_index_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs)
    EM2_obs_index_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM2_index_df$observed = EM2_obs_index_df$value
    EM2_sd_index_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM2_index_df$normalised_resid = (log(EM2_index_df$observed / EM2_index_df$values) + EM2_sd_index_df$values^2) /  EM2_sd_index_df$values
    EM2_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_index_df$model = "EM2"
    #EM2_index_df$RE = (EM2_index_df$values - OM_index_df$values[rep(data_year_bool_ndx, n_sims)]) / OM_index_df$values[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM3_index_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs)
    EM3_obs_index_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM3_index_df$observed = EM3_obs_index_df$value
    EM3_sd_index_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM3_index_df$normalised_resid = (log(EM3_index_df$observed / EM3_index_df$values) + EM3_sd_index_df$values^2) /  EM3_sd_index_df$values
    EM3_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_index_df$model = "EM3"
    #EM3_index_df$RE = (EM3_index_df$values - OM_index_df$values[rep(data_year_bool_ndx, n_sims)]) / OM_index_df$values[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM1_00_index_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs_hist_00)
    EM1_00_obs_index_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM1_00_index_df$observed = EM1_00_obs_index_df$value
    EM1_00_sd_index_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM1_00_index_df$normalised_resid = (log(EM1_00_index_df$observed / EM1_00_index_df$values) + EM1_00_sd_index_df$values^2) /  EM1_00_sd_index_df$values
    EM1_00_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_index_df$model = "EM1_00"
    #EM1_00_index_df$RE = (EM1_00_index_df$values - OM_index_df$values[rep(survey_year_ndx_hist_00, n_sims)]) / OM_index_df$values[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1a_00_index_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs_hist_00)
    EM1a_00_obs_index_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM1a_00_index_df$observed = EM1a_00_obs_index_df$value
    EM1a_00_sd_index_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM1a_00_index_df$normalised_resid = (log(EM1a_00_index_df$observed / EM1a_00_index_df$values) + EM1a_00_sd_index_df$values^2) /  EM1a_00_sd_index_df$values
    EM1a_00_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_index_df$model = "EM1a_00"
    #EM1a_00_index_df$RE = (EM1a_00_index_df$values - OM_index_df$values[rep(survey_year_ndx_hist_00, n_sims)]) / OM_index_df$values[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1b_00_index_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs_hist_00)
    EM1b_00_obs_index_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM1b_00_index_df$observed = EM1b_00_obs_index_df$value
    EM1b_00_sd_index_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM1b_00_index_df$normalised_resid = (log(EM1b_00_index_df$observed / EM1b_00_index_df$values) + EM1b_00_sd_index_df$values^2) /  EM1b_00_sd_index_df$values
    EM1b_00_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_index_df$model = "EM1b_00"
    #EM1b_00_index_df$RE = (EM1b_00_index_df$values - OM_index_df$values[rep(survey_year_ndx_hist_00, n_sims)]) / OM_index_df$values[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM2_00_index_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs_00)
    EM2_00_obs_index_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM2_00_index_df$observed = EM2_00_obs_index_df$value
    EM2_00_sd_index_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM2_00_index_df$normalised_resid = (log(EM2_00_index_df$observed / EM2_00_index_df$values) + EM2_00_sd_index_df$values^2) /  EM2_00_sd_index_df$values
    EM2_00_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_index_df$model = "EM2_00"
    #EM2_00_index_df$RE = (EM2_00_index_df$values - OM_index_df$values[rep(survey_year_ndx_00, n_sims)]) / OM_index_df$values[rep(survey_year_ndx_00, n_sims)] * 100
    
    EM3_00_index_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_index_fitted", element_labs = survey_year_obs_00)
    EM3_00_obs_index_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_obs")
    EM3_00_index_df$observed = EM3_00_obs_index_df$value
    EM3_00_sd_index_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_sd")
    EM3_00_index_df$normalised_resid = (log(EM3_00_index_df$observed / EM3_00_index_df$values) + EM3_00_sd_index_df$values^2) /  EM3_00_sd_index_df$values
    EM3_00_index_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_index_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_index_df$model = "EM3_00"
    #EM3_00_index_df$RE = (EM3_00_index_df$values - OM_index_df$values[rep(survey_year_ndx_00, n_sims)]) / OM_index_df$values[rep(survey_year_ndx_00, n_sims)] * 100
    
    full_index = rbind(full_index, EM1_index_df, EM1a_index_df, EM1b_index_df, EM2_index_df, EM3_index_df, EM1_00_index_df, EM1a_00_index_df, EM1b_00_index_df, EM2_00_index_df, EM3_00_index_df)
  }
}
## match OM_catch to full catch
full_index = full_index %>% left_join(OM_index_fit, by = c("sim_iter", "init", "rebuild", "names"))
## calculate RE
full_index = full_index %>% group_by(sim_iter, init, rebuild, names, model.x) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
full_index$init = factor(full_index$init, levels = paste0("init ", inital_levels), ordered =T)

full_index %>% filter(sim_iter == 1, init == "init 100", rebuild == "rebuild 50", names == 1970)

ggplot(full_index %>% filter(sim_iter == 1), aes(x = names,  linetype = model.x, col = model.x)) +
  geom_point(aes(y = observed.x)) + 
  geom_line(aes(y = values.x), linewidth = 1.1) + 
  theme_bw() +
  facet_grid(rebuild~init)

ggplot(full_index %>% filter(model.x == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(rebuild~init)
ylim(-20,20)

ggplot(full_index %>% filter(model.x == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = normalised_resid)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(init~rebuild) +
  ylim(-3,3)

## Look at the age-frequency
full_survey_AF = NULL
OM_survey_AF = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_AF_df = get_age_fit_by_year(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T, element_labs = survey_year_obs, years = c(min(survey_year_obs), max(survey_year_obs)))
    OM_AF_df$init = paste0("init ", inital_levels[init_ndx])
    OM_AF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_AF_df$model = "OM"
    OM_survey_AF = rbind(OM_survey_AF, OM_AF_df)
    
    EM1_surveyAF_df = get_age_fit_by_year(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs, years = c(min(survey_year_obs), max(survey_year_obs)))
    EM1_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_surveyAF_df$model = "EM1"
    #EM1_surveyAF_df$RE = (EM1_surveyAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100
    
    EM1a_surveyAF_df = get_age_fit_by_year(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T, element_labs = survey_year_obs, years = c(min(survey_year_obs), max(survey_year_obs)))
    EM1a_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_surveyAF_df$model = "EM1a"
    #EM1a_surveyAF_df$RE = (EM1a_surveyAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100   
    
    EM1b_surveyAF_df = get_age_fit_by_year(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs, years = c(min(survey_year_obs), max(survey_year_obs)))
    EM1b_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_surveyAF_df$model = "EM1b"
    #EM1b_surveyAF_df$RE = (EM1b_surveyAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100
    
    EM2_surveyAF_df = get_age_fit_by_year(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs, years = c(min(survey_year_obs), max(survey_year_obs)))
    EM2_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_surveyAF_df$model = "EM2"
    #EM2_surveyAF_df$RE = (EM2_surveyAF_df$Ey - OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)]) / OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM3_surveyAF_df = get_age_fit_by_year(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs, years = c(min(survey_year_obs), max(survey_year_obs)))
    EM3_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_surveyAF_df$model = "EM3"
    #EM3_surveyAF_df$RE = (EM3_surveyAF_df$Ey - OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)]) / OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM1_00_surveyAF_df = get_age_fit_by_year(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T, element_labs = survey_year_obs_hist_00, years = c(min(survey_year_obs_hist_00), max(survey_year_obs_hist_00)))
    EM1_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_surveyAF_df$model = "EM1_00"
    #EM1_00_surveyAF_df$RE = (EM1_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1a_00_surveyAF_df = get_age_fit_by_year(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_hist_00, years = c(min(survey_year_obs_hist_00), max(survey_year_obs_hist_00)))
    EM1a_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_surveyAF_df$model = "EM1a_00"
    #EM1a_00_surveyAF_df$RE = (EM1a_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1b_00_surveyAF_df = get_age_fit_by_year(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_hist_00, years = c(min(survey_year_obs_hist_00), max(survey_year_obs_hist_00)))
    EM1b_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_surveyAF_df$model = "EM1b_00"
    #EM1b_00_surveyAF_df$RE = (EM1b_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM2_00_surveyAF_df = get_age_fit_by_year(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_00, years = c(min(survey_year_obs_00), max(survey_year_obs_00)))
    EM2_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_surveyAF_df$model = "EM2_00"
    #EM2_00_surveyAF_df$RE = (EM2_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)] * 100
    
    EM3_00_surveyAF_df = get_age_fit_by_year(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_00, years = c(min(survey_year_obs_00), max(survey_year_obs_00)))
    EM3_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_surveyAF_df$model = "EM3_00"
    #EM3_00_surveyAF_df$RE = (EM3_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)] * 100
    full_survey_AF = rbind(full_survey_AF, EM1_surveyAF_df, EM1a_surveyAF_df, EM1b_surveyAF_df, EM2_surveyAF_df, EM3_surveyAF_df, EM1_00_surveyAF_df, EM1a_00_surveyAF_df, EM1b_00_surveyAF_df, EM2_00_surveyAF_df, EM3_00_surveyAF_df)
  }
}
## match OM_catch to full catch
full_survey_AF = full_survey_AF %>% left_join(OM_survey_AF, by = c("Age","sim_iter", "init", "rebuild", "year"))
## calculate RE
full_survey_AF = full_survey_AF %>% group_by(Age, sim_iter, init, rebuild, year) %>% mutate(RE = (fit.x - fit.y)/ fit.y * 100)
full_survey_AF$init = factor(full_survey_AF$init, levels = paste0("init ", inital_levels), ordered =T)
## plot the relative error
ggplot(full_survey_AF %>% filter(model.x == "EM1", year == 1981), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM1") +
  ylim(-50, 50) +
  facet_grid(init~rebuild) 

ggplot(full_survey_AF %>% filter(model.x == "EM1", year == 2020), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM1") +
  ylim(-50, 50) +
  facet_grid(init~rebuild) 

ggplot(full_survey_AF %>% filter(model.x == "EM1_00", year == 1961), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM1_00") +
  facet_grid(init~rebuild) 

ggplot(full_survey_AF %>% filter(model.x == "EM2", year == 2020), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM2") +
  facet_grid(init~rebuild) 

## look at fits to observations
## Survey mean age
full_survey_mean_age = NULL
OM_mean_age= NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_AF_df = get_multiple_mean_age(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T, element_labs = survey_year_obs)
    OM_AF_df$init = paste0("init ", inital_levels[init_ndx])
    OM_AF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_AF_df$model = "OM"
    OM_mean_age = rbind(OM_mean_age, OM_AF_df)
    
    EM1_surveyAF_df = get_multiple_mean_age(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs)
    EM1_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_surveyAF_df$model = "EM1"
    #EM1_surveyAF_df$RE = (EM1_surveyAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100

    EM1a_surveyAF_df = get_multiple_mean_age(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T, element_labs = survey_year_obs)
    EM1a_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_surveyAF_df$model = "EM1a"
    #EM1a_surveyAF_df$RE = (EM1a_surveyAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100   
    
    EM1b_surveyAF_df = get_multiple_mean_age(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs)
    EM1b_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_surveyAF_df$model = "EM1b"
    #EM1b_surveyAF_df$RE = (EM1b_surveyAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100
    
    EM2_surveyAF_df = get_multiple_mean_age(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs)
    EM2_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_surveyAF_df$model = "EM2"
    #EM2_surveyAF_df$RE = (EM2_surveyAF_df$Ey - OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)]) / OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM3_surveyAF_df = get_multiple_mean_age(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs)
    EM3_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_surveyAF_df$model = "EM3"
    #EM3_surveyAF_df$RE = (EM3_surveyAF_df$Ey - OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)]) / OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM1_00_surveyAF_df = get_multiple_mean_age(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T, element_labs = survey_year_obs_hist_00)
    EM1_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_surveyAF_df$model = "EM1_00"
    #EM1_00_surveyAF_df$RE = (EM1_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1a_00_surveyAF_df = get_multiple_mean_age(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_hist_00)
    EM1a_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_surveyAF_df$model = "EM1a_00"
    #EM1a_00_surveyAF_df$RE = (EM1a_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1b_00_surveyAF_df = get_multiple_mean_age(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_hist_00)
    EM1b_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_surveyAF_df$model = "EM1b_00"
    #EM1b_00_surveyAF_df$RE = (EM1b_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM2_00_surveyAF_df = get_multiple_mean_age(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_00)
    EM2_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_surveyAF_df$model = "EM2_00"
    #EM2_00_surveyAF_df$RE = (EM2_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)] * 100
    
    EM3_00_surveyAF_df = get_multiple_mean_age(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = T,  element_labs = survey_year_obs_00)
    EM3_00_surveyAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_surveyAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_surveyAF_df$model = "EM3_00"
    #EM3_00_surveyAF_df$RE = (EM3_00_surveyAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)] * 100
    full_survey_mean_age = rbind(full_survey_mean_age, EM1_surveyAF_df, EM1a_surveyAF_df, EM1b_surveyAF_df, EM2_surveyAF_df, EM3_surveyAF_df, EM1_00_surveyAF_df, EM1a_00_surveyAF_df, EM1b_00_surveyAF_df, EM2_00_surveyAF_df, EM3_00_surveyAF_df)
  }
}
## match OM_catch to full catch
full_survey_mean_age = full_survey_mean_age %>% left_join(OM_ssb, by = c("sim_iter", "init", "rebuild", "names"))
## calculate RE
full_survey_mean_age = full_survey_mean_age %>% group_by(sim_iter, init, rebuild, names, model.x) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
full_survey_mean_age$init = factor(full_survey_mean_age$init, levels = paste0("init ", inital_levels), ordered =T)
## plot the relative error
ggplot(full_survey_mean_age %>% filter(sim_iter == 1, model == "EM1_00"), aes(x = year,  linetype = model, col = model)) +
  geom_point(aes(y = Oy)) + 
  geom_line(aes(y = Ey), linewidth = 1.1) + 
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(init~rebuild) 

## look at fits to observations
## fishery mean age
full_fishery_mean_age = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_AF_df = get_multiple_mean_age(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F, element_labs = survey_year_obs)
    OM_AF_df$init = paste0("init ", inital_levels[init_ndx])
    OM_AF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_AF_df$model = "OM"
    
    EM1_fisheryAF_df = get_multiple_mean_age(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs)
    EM1_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_fisheryAF_df$model = "EM1"
    EM1_fisheryAF_df$RE = (EM1_fisheryAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100
    
    EM1a_fisheryAF_df = get_multiple_mean_age(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F, element_labs = survey_year_obs)
    EM1a_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_fisheryAF_df$model = "EM1a"
    EM1a_fisheryAF_df$RE = (EM1a_fisheryAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100   
    
    EM1b_fisheryAF_df = get_multiple_mean_age(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs)
    EM1b_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_fisheryAF_df$model = "EM1b"
    EM1b_fisheryAF_df$RE = (EM1b_fisheryAF_df$Ey - OM_AF_df$Ey) / OM_AF_df$Ey * 100
    
    EM2_fisheryAF_df = get_multiple_mean_age(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs)
    EM2_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_fisheryAF_df$model = "EM2"
    EM2_fisheryAF_df$RE = (EM2_fisheryAF_df$Ey - OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)]) / OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM3_fisheryAF_df = get_multiple_mean_age(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs)
    EM3_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_fisheryAF_df$model = "EM3"
    EM3_fisheryAF_df$RE = (EM3_fisheryAF_df$Ey - OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)]) / OM_AF_df$Ey[rep(data_year_bool_ndx, n_sims)] * 100
    
    EM1_00_fisheryAF_df = get_multiple_mean_age(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F, element_labs = survey_year_obs_hist_00)
    EM1_00_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_fisheryAF_df$model = "EM1_00"
    EM1_00_fisheryAF_df$RE = (EM1_00_fisheryAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1a_00_fisheryAF_df = get_multiple_mean_age(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs_hist_00)
    EM1a_00_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_fisheryAF_df$model = "EM1a_00"
    EM1a_00_fisheryAF_df$RE = (EM1a_00_fisheryAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM1b_00_fisheryAF_df = get_multiple_mean_age(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs_hist_00)
    EM1b_00_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_fisheryAF_df$model = "EM1b_00"
    EM1b_00_fisheryAF_df$RE = (EM1b_00_fisheryAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_hist_00, n_sims)] * 100
    
    EM2_00_fisheryAF_df = get_multiple_mean_age(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs_00)
    EM2_00_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_fisheryAF_df$model = "EM2_00"
    EM2_00_fisheryAF_df$RE = (EM2_00_fisheryAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)] * 100
    
    EM3_00_fisheryAF_df = get_multiple_mean_age(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], survey = F,  element_labs = survey_year_obs_00)
    EM3_00_fisheryAF_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_fisheryAF_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_fisheryAF_df$model = "EM3_00"
    EM3_00_fisheryAF_df$RE = (EM3_00_fisheryAF_df$Ey - OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)]) / OM_AF_df$Ey[rep(survey_year_ndx_00, n_sims)] * 100
    
    full_fishery_mean_age = rbind(full_fishery_mean_age, EM1_fisheryAF_df, EM1a_fisheryAF_df, EM1b_fisheryAF_df, EM2_fisheryAF_df, EM3_fisheryAF_df, EM1_00_fisheryAF_df, EM1a_00_fisheryAF_df, EM1b_00_fisheryAF_df, EM2_00_fisheryAF_df, EM3_00_fisheryAF_df)
  }
}
full_fishery_mean_age$init = factor(full_fishery_mean_age$init, levels = paste0("init ", inital_levels), ordered =T)

ggplot(full_fishery_mean_age %>% filter(sim_iter == 1, model == "EM1_00"), aes(x = year,  linetype = model, col = model)) +
  geom_point(aes(y = Oy)) + 
  geom_line(aes(y = Ey), linewidth = 1.1) + 
  theme_bw() +
  ggtitle("EM1") +
  facet_grid(init~rebuild) 

## look at fit to age for first year and last year.



## Check out survey selectivity
full_survey_selectivity = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_survey_selectivity_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = ages)
    OM_survey_selectivity_df$values = OM_survey_selectivity_df$values * 100
    OM_survey_selectivity_df$init =  paste0("init ", inital_levels[init_ndx])
    OM_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_survey_selectivity_df$model = "OM"
    ## EM1
    EM1_survey_selectivity_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = ages)
    EM1_survey_selectivity_df$values = EM1_survey_selectivity_df$values * 100
    EM1_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_survey_selectivity_df$model = "EM1"
    EM1_survey_selectivity_df$RE = (EM1_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM1a
    EM1a_survey_selectivity_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = ages)
    EM1a_survey_selectivity_df$values = EM1a_survey_selectivity_df$values * 100
    EM1a_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_survey_selectivity_df$model = "EM1a"
    EM1a_survey_selectivity_df$RE = (EM1a_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM1b
    EM1b_survey_selectivity_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = ages)
    EM1b_survey_selectivity_df$values = EM1b_survey_selectivity_df$values * 100
    EM1b_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_survey_selectivity_df$model = "EM1b"
    EM1b_survey_selectivity_df$RE = (EM1b_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM2
    EM2_survey_selectivity_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c( ages))
    EM2_survey_selectivity_df$values = EM2_survey_selectivity_df$values * 100
    EM2_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_survey_selectivity_df$RE = (EM2_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    EM2_survey_selectivity_df$model = "EM2"
    ## EM3
    EM3_survey_selectivity_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c(ages))
    EM3_survey_selectivity_df$values = EM3_survey_selectivity_df$values * 100
    EM3_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_survey_selectivity_df$model = "EM3"
    EM3_survey_selectivity_df$RE = (EM3_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM1
    EM1_00_survey_selectivity_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c( ages))
    EM1_00_survey_selectivity_df$values = EM1_00_survey_selectivity_df$values * 100
    EM1_00_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_survey_selectivity_df$model = "EM1_00"
    EM1_00_survey_selectivity_df$RE = (EM1_00_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM1a
    EM1a_00_survey_selectivity_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c(ages))
    EM1a_00_survey_selectivity_df$values = EM1a_00_survey_selectivity_df$values * 100
    EM1a_00_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_survey_selectivity_df$model = "EM1a_00"
    EM1a_00_survey_selectivity_df$RE = (EM1a_00_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM1b
    EM1b_00_survey_selectivity_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c(ages))
    EM1b_00_survey_selectivity_df$values = EM1b_00_survey_selectivity_df$values * 100
    EM1b_00_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_survey_selectivity_df$model = "EM1b_00"
    EM1b_00_survey_selectivity_df$RE = (EM1b_00_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    ## EM2
    EM2_00_survey_selectivity_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c(ages))
    EM2_00_survey_selectivity_df$values = EM2_00_survey_selectivity_df$values * 100
    EM2_00_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_survey_selectivity_df$RE = (EM2_00_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    EM2_00_survey_selectivity_df$model = "EM2_00"
    ## EM3
    EM3_00_survey_selectivity_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "survey_selectivity", element_labs = c(ages))
    EM3_00_survey_selectivity_df$values = EM3_00_survey_selectivity_df$values * 100
    EM3_00_survey_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_survey_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_survey_selectivity_df$model = "EM3_00"
    EM3_00_survey_selectivity_df$RE = (EM3_00_survey_selectivity_df$values - OM_survey_selectivity_df$values) / OM_survey_selectivity_df$values * 100
    
    full_survey_selectivity = rbind(full_survey_selectivity, EM1_survey_selectivity_df, EM1a_survey_selectivity_df, EM1b_survey_selectivity_df, EM2_survey_selectivity_df, EM3_survey_selectivity_df,EM1_00_survey_selectivity_df, EM1a_00_survey_selectivity_df, EM1b_00_survey_selectivity_df, EM2_00_survey_selectivity_df, EM3_00_survey_selectivity_df)
  }
}
init_survey_selectivity_lvls = paste0("init ", inital_levels)
full_survey_selectivity$init = factor(full_survey_selectivity$init, levels = init_survey_selectivity_lvls, ordered =T)

ggplot(full_survey_selectivity %>% filter(model == "EM1_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM3") +
  labs(x = "age", y = "RE") +
  facet_grid(init~rebuild)

## Get a range of reference points
full_fishery_selectivity = NULL
for(init_ndx in 1:length(inital_levels)) {
  for(rebuild_ndx in 1:length(rebuild_levels)) {
    OM_fishery_selectivity_df = get_multiple_vectors(OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = ages)
    OM_fishery_selectivity_df$values = OM_fishery_selectivity_df$values * 100
    OM_fishery_selectivity_df$init =  paste0("init ", inital_levels[init_ndx])
    OM_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    OM_fishery_selectivity_df$model = "OM"
    ## EM1
    EM1_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM1[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = ages)
    EM1_fishery_selectivity_df$values = EM1_fishery_selectivity_df$values * 100
    EM1_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_fishery_selectivity_df$model = "EM1"
    EM1_fishery_selectivity_df$RE = (EM1_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM1a
    EM1a_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM1a[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = ages)
    EM1a_fishery_selectivity_df$values = EM1a_fishery_selectivity_df$values * 100
    EM1a_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_fishery_selectivity_df$model = "EM1a"
    EM1a_fishery_selectivity_df$RE = (EM1a_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM1b
    EM1b_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM1b[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = ages)
    EM1b_fishery_selectivity_df$values = EM1b_fishery_selectivity_df$values * 100
    EM1b_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_fishery_selectivity_df$model = "EM1b"
    EM1b_fishery_selectivity_df$RE = (EM1b_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM2
    EM2_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c( ages))
    EM2_fishery_selectivity_df$values = EM2_fishery_selectivity_df$values * 100
    EM2_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_fishery_selectivity_df$RE = (EM2_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    EM2_fishery_selectivity_df$model = "EM2"
    ## EM3
    EM3_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM3[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c(ages))
    EM3_fishery_selectivity_df$values = EM3_fishery_selectivity_df$values * 100
    EM3_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_fishery_selectivity_df$model = "EM3"
    EM3_fishery_selectivity_df$RE = (EM3_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM1
    EM1_00_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM1_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c( ages))
    EM1_00_fishery_selectivity_df$values = EM1_00_fishery_selectivity_df$values * 100
    EM1_00_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1_00_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1_00_fishery_selectivity_df$model = "EM1_00"
    EM1_00_fishery_selectivity_df$RE = (EM1_00_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM1a
    EM1a_00_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM1a_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c(ages))
    EM1a_00_fishery_selectivity_df$values = EM1a_00_fishery_selectivity_df$values * 100
    EM1a_00_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1a_00_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1a_00_fishery_selectivity_df$model = "EM1a_00"
    EM1a_00_fishery_selectivity_df$RE = (EM1a_00_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM1b
    EM1b_00_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM1b_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c(ages))
    EM1b_00_fishery_selectivity_df$values = EM1b_00_fishery_selectivity_df$values * 100
    EM1b_00_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM1b_00_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM1b_00_fishery_selectivity_df$model = "EM1b_00"
    EM1b_00_fishery_selectivity_df$RE = (EM1b_00_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    ## EM2
    EM2_00_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM2_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c(ages))
    EM2_00_fishery_selectivity_df$values = EM2_00_fishery_selectivity_df$values * 100
    EM2_00_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM2_00_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM2_00_fishery_selectivity_df$RE = (EM2_00_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    EM2_00_fishery_selectivity_df$model = "EM2_00"
    ## EM3
    EM3_00_fishery_selectivity_df = get_multiple_vectors(mle_lst_EM3_00[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]], "fishery_selectivity", element_labs = c(ages))
    EM3_00_fishery_selectivity_df$values = EM3_00_fishery_selectivity_df$values * 100
    EM3_00_fishery_selectivity_df$init = paste0("init ", inital_levels[init_ndx])
    EM3_00_fishery_selectivity_df$rebuild = paste0("rebuild ", rebuild_levels[rebuild_ndx])
    EM3_00_fishery_selectivity_df$model = "EM3_00"
    EM3_00_fishery_selectivity_df$RE = (EM3_00_fishery_selectivity_df$values - OM_fishery_selectivity_df$values) / OM_fishery_selectivity_df$values * 100
    
    full_fishery_selectivity = rbind(full_fishery_selectivity, EM1_fishery_selectivity_df, EM1a_fishery_selectivity_df, EM1b_fishery_selectivity_df, EM2_fishery_selectivity_df, EM3_fishery_selectivity_df,EM1_00_fishery_selectivity_df, EM1a_00_fishery_selectivity_df, EM1b_00_fishery_selectivity_df, EM2_00_fishery_selectivity_df, EM3_00_fishery_selectivity_df)
  }
}
init_fishery_selectivity_lvls = paste0("init ", inital_levels)
full_fishery_selectivity$init = factor(full_fishery_selectivity$init, levels = init_fishery_selectivity_lvls, ordered =T)

ggplot(full_fishery_selectivity %>% filter(model == "EM1_00"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM3") +
  labs(x = "age", y = "RE") +
  facet_grid(init~rebuild)

ggplot(full_fishery_selectivity %>% filter(model == "EM1"), aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EM1") +
  labs(x = "age", y = "RE") +
  facet_grid(init~rebuild)
