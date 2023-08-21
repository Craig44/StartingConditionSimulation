#'
#' See why EM2 and EM3 can't recreate OM's that are depleted
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
compile(file = file.path(DIR$tmb, "AgeStructuredModel_tmp.cpp"), flags = "-Wignored-attributes -O3")

#sink()
#dyn.unload(dynlib(file.path(DIR$tmb, "AgeStructuredModel_tmp")))

dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModel_tmp")))
dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))

#setwd(DIR$R)

fig_dir = file.path(DIR$fig, "illustrateEM2_EM3")
if(!dir.exists(fig_dir))
  dir.create(fig_dir)

this_bio = readRDS(file = file.path(DIR$data, "Medium_biology.RDS"))
##

## data period 
n_years_data = 40
historic_years = seq(from = 25, to  = 150, by = 25) 
EM2_lst = EM3_lst = OM_lst = mle_lst_EM2 = mle_lst_EM3 = list();
## generate random recruitment
rec_devs =  rnorm(sum(TMB_data$ycs_estimated),  -0.5 * this_bio$sigma_r * this_bio$sigma_r, this_bio$sigma_r)
for(hist_year_ndx in 1:length(historic_years)) {
  
  n_years_historic = historic_years[hist_year_ndx]
    
  
  n_years = n_years_historic + n_years_data
  years = (2020 - n_years + 1):2020
  full_years = (min(years) - 1):max(years)
  ## observation temporal frequency
  survey_year_obs = 1980:2020
  survey_ages = this_bio$ages
  fishery_year_obs = 1980:2020
  fishery_ages = this_bio$ages
  
  ############
  ## Build OM
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
  TMB_data$ycs_bias_correction = c(rep(1, n_years))
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
  
  fishery_probs = c(seq(from = 0.01, to = this_bio$M * 1.5, length = n_years_historic), rep(this_bio$M, n_years_data))
  plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)
  
  ## The same parameters as OM, to check for consistency
  OM_pars = list(
    ln_R0 = log(this_bio$R0),
    ln_ycs_est =  rec_devs,
    ln_sigma_r = log( this_bio$sigma_r),
    #ln_extra_survey_cv = log(0.0001),
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
  OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
  OM_report = OM_obj$report()
  
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
  
  F_max = find_F_max(fishery_sel, M, waa_catch, ages)
  F_0.1 = find_F_0.1(fishery_sel, M, waa_catch, ages, derivative_method = 1)
  F_0.1_alt = find_F_0.1(fishery_sel, M, waa_catch, ages, derivative_method = 2)
  F_msy = find_F_msy(fishery_sel = fishery_sel, M = M, catch_waa = waa_catch, ssb_waa = waa, paa = paa, ages = ages, Ninit = OM_report$equilibrium_at_age, plus_group = T, prop_Z = 0.5, n_runs = 200)
  
  ## 25 year data scenario
  fishery_probs = c(seq(from = 0.01, to = F_15$F_ref, length = n_years_historic), seq(from = F_25$F_ref, to = F_40$F_ref, length = 10), rep(F_40$F_ref, n_years_data - 10))
  plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)
  
  TMB_data$fixed_catch_indicator = rep(0, n_years)
  #TMB_data$fixed_catch_indicator[40:n_years] = 1
  OM_pars$ln_F = array(log(fishery_probs), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
  OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
  OM_report = OM_obj$report()
  
  ##
  ## 
  ##
  f_100 = 0.0000001
  fish_100_50 = c(rep(f_100,n_years_historic), rep(F_50$F_ref, n_years_data))
  fish_100_35 = c(rep(f_100,n_years_historic), rep(F_35$F_ref, n_years_data))
  fish_100_20 = c(rep(f_100,n_years_historic), rep(F_20$F_ref, n_years_data))
  
  fish_75_50 = c(rep(init_F_75$F_ref,n_years_historic), rep(F_50$F_ref, n_years_data))
  fish_75_35 = c(rep(init_F_75$F_ref,n_years_historic), rep(F_35$F_ref, n_years_data))
  fish_75_20 = c(rep(init_F_75$F_ref,n_years_historic), rep(F_20$F_ref, n_years_data))
  
  fish_50_50 = c(rep(init_F_50$F_ref,n_years_historic), rep(F_50$F_ref, n_years_data))
  fish_50_35 = c(rep(init_F_50$F_ref,n_years_historic), rep(F_35$F_ref, n_years_data))
  fish_50_20 = c(rep(init_F_50$F_ref,n_years_historic), rep(F_20$F_ref, n_years_data))
  
  fish_25_50 = c(rep(init_F_25$F_ref,n_years_historic), rep(F_50$F_ref, n_years_data))
  fish_25_35 = c(rep(init_F_25$F_ref,n_years_historic), rep(F_35$F_ref, n_years_data))
  fish_25_20 = c(rep(init_F_25$F_ref,n_years_historic), rep(F_20$F_ref, n_years_data))
  
  
  ## start with the extreme case to help understand
  OM_pars$ln_F = array(log(fish_25_35), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
  
  OM_25_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
  OM_25_35_rep = OM_25_35$report()
  sim_data = OM_25_35$simulate(complete = T)
  
  n_last_ycs_to_estimate = 1;
  n_first_ycs_to_estimate = 1; ## you may want to shift this depending on when data starts recruits are observed
  n_first_ycs_to_estimate_for_historic_models = 1; ## assume the first 10 year YCS are fixed at YCS
  
  N_eff = 150
  data_years = years[(n_years_historic + 1):length(years)]
  data_year_ndx = which(years %in% data_years)
  data_full_years = (min(data_years) - 1):max(data_years)
  EM_short_data = TMB_data
  EM_short_data$F_method = 0
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
  EM_short_data$survey_AF_obs = array(N_eff/TMB_data$max_age_ndx_for_AFs, dim = c(TMB_data$max_age_ndx_for_AFs, sum(short_survey_year_ndx)))
  EM_short_data$fishery_AF_obs = array(N_eff/TMB_data$max_age_ndx_for_AFs, dim = c(TMB_data$max_age_ndx_for_AFs, sum(short_fishery_year_ndx), EM_short_data$n_fisheries))
  EM_short_data$catches = array(EM_short_data$catches[EM_short_data$n_years,], dim = c(EM_short_data$n_years, EM_short_data$n_fisheries))
  #EM_short_data$ycs_estimated = rep(0, EM_short_data$n_years)
  #EM_short_data$ycs_estimated[n_first_ycs_to_estimate:(length(EM_short_data$ycs_estimated) - n_last_ycs_to_estimate )] = 1
  
  EM_short_data$catch_indicator = array(1, dim = c(EM_short_data$n_years, EM_short_data$n_fisheries))
  EM_short_data$stockMeanLength = EM_short_data$catchMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = EM_short_data$n_years, nrow = EM_short_data$n_ages)
  EM_short_data$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = EM_short_data$n_years, nrow = EM_short_data$n_ages)
  EM_short_data$propZ_ssb = EM_short_data$propZ_survey = rep(0.5, EM_short_data$n_years)
  EM_short_pars = OM_pars
  EM_short_pars$ln_ycs_est = rep(0, sum(EM_short_data$ycs_estimated))
  EM_short_pars$ln_F = array(log(0.1), dim = c(EM_short_data$n_fisheries, EM_short_data$n_years))
  
  EM2_data = EM3_data = EM_short_data
  
  EM2_data$catches = EM3_data$catches = matrix(sim_data$catches[data_year_ndx, ], ncol = 1)
  EM2_data$survey_obs = EM3_data$survey_obs = sim_data$survey_obs[short_survey_year_ndx]
  EM2_data$survey_cv = EM3_data$survey_cv = sim_data$survey_cv[short_survey_year_ndx]
  EM2_data$survey_AF_obs = EM3_data$survey_AF_obs = sim_data$survey_AF_obs[,short_survey_year_ndx]
  EM2_data$fishery_AF_obs = EM3_data$fishery_AF_obs = array(sim_data$fishery_AF_obs[,short_fishery_year_ndx,], dim = c(dim(sim_data$fishery_AF_obs)[1],sum(short_fishery_year_ndx), sim_data$n_fisheries))
  EM2_data$estimate_F_init = EM3_data$estimate_F_init = 1
  EM2_data$estimate_init_age_devs = 0
  EM3_data$estimate_init_age_devs = 1
  EM3_data$n_init_age_devs = n_years_historic - 2
  ## set true values and have a look
  EM_short_pars$ln_F_init = log(init_F_25$F_ref)
  EM_short_pars$ln_F = matrix(log(rep(F_35$F_ref, n_years_data)), ncol = 1)
  EM_short_pars$ln_ycs_est = OM_pars$ln_ycs_est[(n_years_historic + 1):length(OM_pars$ln_ycs_est)]
  EM3_short_pars = EM_short_pars
  
  ## calcualte true age-deviations
  init_devs = log((sim_data$N[,n_years_historic]/sum(sim_data$N[,n_years_historic] ) ) / (sim_data$equilibrium_at_age / sum(sim_data$equilibrium_at_age)))
  EM3_short_pars$ln_init_age_devs = init_devs[2:(EM3_data$n_init_age_devs + 1)]
  
  
  na_EM2_pars = fix_pars(EM_short_pars, pars_to_exclude = c("ln_sigma_r",  "ln_init_age_devs","ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
  na_EM3_pars = fix_pars(EM3_short_pars, pars_to_exclude = c("ln_sigma_r",  "ln_sigma_init_age_devs", "ln_F", "ln_catch_sd", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
  
  EM2_obj <- MakeADFun(EM2_data, EM_short_pars, map = na_EM2_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T, silent = T)
  EM2_rep = EM2_obj$report()
  EM3_obj <- MakeADFun(EM3_data, EM3_short_pars, map = na_EM3_pars, DLL= "AgeStructuredModel_tmp", checkParameterOrder = T, silent = T)
  EM3_rep = EM3_obj$report()
  par(mfrow = c(3,1))
  plot(full_years, sim_data$depletion, type = "l", lwd = 3, xlab = "years", ylab = "Depletion", ylim = c(0,1))
  lines(data_full_years, EM2_rep$depletion, col = "blue", lty = 2, lwd =3)
  lines(data_full_years, EM3_rep$depletion, col = "red", lty = 3, lwd =3)
  
  plot(ages, sim_data$N[,n_years_historic] / 1000, type = "l", lwd = 3, xlab = "ages", ylab = "frequency", main = n_years_historic)
  lines(ages, sim_data$equilibrium_at_age / 1000, col = "gray", lty = 2, lwd =3)
  lines(ages, EM2_rep$N[,1] / 1000, col = "blue", lty = 2, lwd =3)
  lines(ages, EM3_rep$N[,1] / 1000, col = "red", lty = 3, lwd =3)
  
  plot(ages, sim_data$N[,n_years_historic + 4] / 1000, type = "l", lwd = 3, xlab = "ages", ylab = "frequency", main = n_years_historic)
  lines(ages, EM2_rep$N[,1 + 4] / 1000, col = "blue", lty = 2, lwd =3)
  lines(ages, EM3_rep$N[,1 + 4] / 1000, col = "red", lty = 3, lwd =3)
  
  Sys.sleep(2);
  EM2_lst[[hist_year_ndx]] = EM2_rep
  EM3_lst[[hist_year_ndx]] = EM3_rep
  OM_lst[[hist_year_ndx]] = sim_data
  
  ## now estimate
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
    cat("hist_year_ndx ", hist_year_ndx, "\n")
    
    next;
  } else {
    ## save the output
    EM2_rep = EM2_obj$report(mle_EM2$par)
    mle_lst_EM2[[hist_year_ndx]]  = EM2_rep
    EM2_se = sdreport(EM2_obj)
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
    cat("hist_year_ndx ", hist_year_ndx, "\n")
    next;
  } else {
    ## save the output
    EM3_rep = EM3_obj$report(mle_EM3$par)
    mle_lst_EM3[[hist_year_ndx]] = EM3_rep
    EM3_se = sdreport(EM3_obj)

  }
  plot(ages, sim_data$N[,n_years_historic] / 1000, type = "l", lwd = 3, xlab = "ages", ylab = "frequency", main = n_years_historic)
  lines(ages, EM2_rep$N[,1] / 1000, col = "blue", lty = 2, lwd =3)
  lines(ages, EM3_rep$N[,1] / 1000, col = "red", lty = 3, lwd =3)
  
  #EM2_rep$B0
  #EM3_rep$B0
  #sim_data$B0
  
}

## show how ages being exposed to multiple fishing events is 
## different from an instantanesous parameterisation
F_init = 0.2
Z=F_init*biol_df$fishery_selectivity+this_bio$M
N.pr1=rep(1,nages) #Number of spawners per recruit at age
for (a in 1:(nages-1)) {
  N.pr1[a+1]=N.pr1[a]*exp(-Z[a])
}
N.pr1[nages]=N.pr1[nages]/(1-exp(-Z[nages]))  #Plus group

exp(- (1:nages - 1) * Z[a])

plot(N.pr1, type = "l")

## vs equilibrium initialisation and applying that same z
N.pr0=rep(1,nages) #Number of spawners per recruit at age
for (a in 1:(nages-1)) {
  N.pr0[a+1]=N.pr0[a]*exp(-this_bio$M)
}
N.pr0[nages]=N.pr0[nages]/(1-exp(-this_bio$M))  #Plus group

## apply that Z for 20 years
n_years = 20
N.pr2=N.pr0
for (y in 1:n_years) {
  tmp_part = N.pr2
  tmp_part[1] = 1 ## recruitment
  for (a in 1:(nages-1)) {
    tmp_part[a+1]=N.pr2[a]*exp(-Z[a])
  }
  tmp_part[nages] = tmp_part[nages] + N.pr2[nages]*exp(-Z[nages])
  N.pr2 = tmp_part
}
  

