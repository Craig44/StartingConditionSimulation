#'
#' Run a profile for B0 and F-init for EM2 model
#' for the medium biology scenario to provide
#' insight into what data is informing these parameters
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

## test
OM_short<- MakeADFun(EM_short_data, EM_short_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_short$fn()

## OM simulations
## --> Initial
## ----> rebuild 
## ------> Simulations
SSB_df = recruit_df = depletion_df = NULL
inital_levels = c(25, 100)# c(25, 50, 75, 100)
rebuild_levels = c(20, 50) # c(20, 35, 50)
EM2_data = EM_short_data
EM2_data$estimate_F_init = 1
## sort out parameters
na_EM2_pars = fix_pars(EM_short_pars, pars_to_exclude = c("ln_sigma_r",  "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
r0_profile_map = fix_pars(EM_short_pars, pars_to_exclude = c("ln_R0", "ln_sigma_r",  "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
F_init_profile_map = fix_pars(EM_short_pars, pars_to_exclude = c("ln_F_init", "ln_sigma_r",  "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))

## test these pars
test <- MakeADFun(EM2_data, EM_short_pars, map = na_EM2_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)

EM2_convergence = array(T, dim = c(length(inital_levels), length(rebuild_levels), n_sims))
mle_lst_EM2  = list()
OM_sim_lst = OM_rep_lst = list()
EM_short_pars$ln_F_init = log(0.07)

#####
# Start Simulations
# 
#
R0_profiles = F_init_profiles = list()
ln_R0_vals = seq(from = 11, to = 16, by = 0.2)
ln_finit_vals = log(seq(from = 0.01, to = 0.1, by = 0.01))

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
    ## Run 1 simulation
    sim_iter = 1
    ## simulate YCS parameters
    OM_pars$ln_ycs_est = rnorm(sum(TMB_data$ycs_estimated), -0.5 * exp(OM_pars$ln_sigma_r) * exp(OM_pars$ln_sigma_r), exp(OM_pars$ln_sigma_r))
    ## Build OM and simulate parameters
    OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
    OM_sim <- OM_obj$simulate(complete = T)
    OM_sim_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = OM_sim
    ## Set catch for EM's
    EM2_data$catches = matrix(OM_sim$catches[data_year_ndx, ], ncol = 1)

    EM2_data$survey_obs =  OM_sim$survey_obs[short_survey_year_ndx]
    EM2_data$survey_cv = OM_sim$survey_cv[short_survey_year_ndx]
    EM2_data$survey_AF_obs = OM_sim$survey_AF_obs[,short_survey_year_ndx]
    EM2_data$fishery_AF_obs = array(OM_sim$fishery_AF_obs[,short_fishery_year_ndx,], dim = c(dim(OM_sim$fishery_AF_obs)[1],sum(short_fishery_year_ndx), OM_sim$n_fisheries))
    ## estimation
    EM2_obj <- MakeADFun(EM2_data, EM_short_pars, map = na_EM2_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
    SpatialSablefishAssessment::check_gradients(EM2_obj)
    OM_tmp_obj <- MakeADFun(OM_sim, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
    OM_tmp_rep = OM_tmp_obj$report()
    OM_rep_lst[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(sim_iter)]] = OM_tmp_rep

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

    } else {
      ## save the output
      EM2_rep = EM2_obj$report(mle_EM2$par)
      mle_lst_EM2[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]] = EM2_rep
    }
    ## run profiles
    r0_pars = EM_short_pars;
    for(r0_ndx in 1:length(ln_R0_vals)) {
      r0_pars$ln_R0 = ln_R0_vals[r0_ndx]
      EM2_obj_R0 <- MakeADFun(EM2_data, r0_pars, map = r0_profile_map, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      mle_EM2_R0 = nlminb(start = EM2_obj_R0$par, objective = EM2_obj_R0$fn, gradient  = EM2_obj_R0$gr, control = list(iter.max = 10000, eval.max = 10000))
      ## check positive definite hessian
      g = as.numeric(EM2_obj_R0$gr(mle_EM2_R0$par))
      hessian = tryCatch(expr = optimHess(mle_EM2_R0$par, fn = EM2_obj_R0$fn, gr = EM2_obj_R0$gr), error = function(e){e})
      if(inherits(hessian, "error") | inherits(hessian, "warning")) {
        cat("R0 ndx = ", r0_ndx, " hessian not PD")
        next;
      }
      EM2_rep_r0 = EM2_obj_R0$report(mle_EM2_R0$par)
      R0_profiles[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(r0_ndx)]] = EM2_rep_r0
    }
    Finit_pars = EM_short_pars;
    for(finit_ndx in 1:length(ln_finit_vals)) {
      Finit_pars$ln_F_init  = ln_finit_vals[finit_ndx]
      EM2_obj_Finit <- MakeADFun(EM2_data, Finit_pars, map = F_init_profile_map, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
      mle_EM2_Finit = nlminb(start = EM2_obj_Finit$par, objective = EM2_obj_Finit$fn, gradient  = EM2_obj_Finit$gr, control = list(iter.max = 10000, eval.max = 10000))
      ## check positive definite hessian
      g = as.numeric(EM2_obj_Finit$gr(mle_EM2_Finit$par))
      hessian = tryCatch(expr = optimHess(mle_EM2_Finit$par, fn = EM2_obj_Finit$fn, gr = EM2_obj_Finit$gr), error = function(e){e})
      if(inherits(hessian, "error") | inherits(hessian, "warning")) {
        cat("Finit ndx = ", finit_ndx, " hessian not PD")
        next;
      }
      EM2_rep_finit = EM2_obj_Finit$report(mle_EM2_Finit$par)
      F_init_profiles[[as.character(inital_levels[init_ndx])]][[as.character(rebuild_levels[rebuild_ndx])]][[as.character(finit_ndx)]] = EM2_rep_finit
    }
  }
}

