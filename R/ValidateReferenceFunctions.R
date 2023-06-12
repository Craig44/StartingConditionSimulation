#'
#' Validate the reference point functions used in this project
#'
source("AuxillaryFunctions.R")
library(dplyr)
library(ggplot2)
library(reshape2)
library(TMB)

## Pass the OM generated data to the TMB model
#sink(file = "compile_output.txt")
compile(file = file.path(DIR$tmb, "AgeStructuredModel.cpp"), flags = "-Wignored-attributes -O3")
#sink()
#dyn.unload(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))
dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))
#setwd(DIR$R)


this_bio = readRDS(file = file.path(DIR$data, "Medium_biology.RDS"))
#this_bio = readRDS(file = file.path(DIR$data, "Fast_biology.RDS"))


n_years = 160
years = (2020 - n_years + 1):2020
full_years = (min(years) - 1):max(years)
## observation temporal frequency
survey_year_obs = years
survey_ages = this_bio$ages
fishery_year_obs = years
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
TMB_data$survey_obs = rnorm(sum(TMB_data$survey_year_indicator), 100, 4)
TMB_data$survey_cv = rep(0.15,sum(TMB_data$survey_year_indicator))
TMB_data$survey_sample_time = rep(0.5,sum(TMB_data$survey_year_indicator))
TMB_data$survey_AF_obs = matrix(5, nrow = TMB_data$n_ages, ncol = sum(TMB_data$survey_year_indicator))

TMB_data$fishery_year_indicator = array(as.integer(TMB_data$years %in% fishery_year_obs), dim = c(length(fishery_year_obs), TMB_data$n_fisheries))
TMB_data$fishery_AF_obs = array(5, dim = c(TMB_data$n_ages, length(fishery_year_obs), TMB_data$n_fisheries))

TMB_data$catches = array(1000, dim = c(TMB_data$n_years, TMB_data$n_fisheries))# this will be overriden in the simulate() call
TMB_data$F_method = 0
TMB_data$F_iterations = 4
TMB_data$F_max = 3

TMB_data$catch_indicator = array(1, dim = c(TMB_data$n_years, TMB_data$n_fisheries))
TMB_data$ycs_estimated = rep(1, n_years)
TMB_data$standardise_ycs = 0;

TMB_data$catchMeanLength = TMB_data$stockMeanLength = matrix(vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0), byrow = F, ncol = TMB_data$n_years, nrow = TMB_data$n_ages)
TMB_data$propMat = matrix(logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95), byrow = F, ncol = TMB_data$n_years, nrow = TMB_data$n_ages)
TMB_data$natMor = this_bio$M
TMB_data$steepness = this_bio$h
TMB_data$stockRecruitmentModelCode = 2 ## BH
TMB_data$propZ_ssb = rep(0.5, TMB_data$n_years)
TMB_data$propZ_survey = rep(0.5, TMB_data$n_years)
TMB_data$sel_ato95_bounds = c(0.1,20)
TMB_data$sel_a50_bounds = c(0.1,20)
TMB_data$mean_weight_a = this_bio$a
TMB_data$mean_weight_b = this_bio$b
TMB_data$estimate_F_init = 0
TMB_data$estimate_init_age_devs = 0
TMB_data$n_init_age_devs = 1
TMB_data$rec_devs_sum_to_zero = 0
TMB_data$Q_r_for_sum_to_zero = Q_sum_to_zero_QR(length(TMB_data$years))

## fishery_probs
fishery_probs = c(rep(this_bio$M * 1.5, 20), rep(this_bio$M, 40), rep(this_bio$M, 100) )
plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)

## The same parameters as OM, to check for consistency
OM_pars = list(
  ln_R0 = log(this_bio$R0),
  ln_ycs_est =  rep(0, sum(TMB_data$ycs_estimated)),#rnorm(sum(TMB_data$ycs_estimated),  0, this_bio$sigma_r),
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

waa_old = waa
paa_old = paa
ages = ages
fishery_sel_old = fishery_sel

F_30 = find_F_percent(target_spr = 30, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_30_old = F_30
F_35 = find_F_percent(target_spr = 35, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_40 = find_F_percent(target_spr = 40, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_20 = find_F_percent(target_spr = 20, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_20_old = F_20
F_50 = find_F_percent(target_spr = 50, fishery_sel, M, waa, paa, ages, prop_Z = 0.5)
F_50_old = F_50

F_max = find_F_max(fishery_sel, M, waa_catch, ages)
F_0.1 = find_F_0.1(fishery_sel, M, waa_catch, ages, derivative_method = 1)
F_0.1_alt = find_F_0.1(fishery_sel, M, waa_catch, ages, derivative_method = 2)
F_msy = find_F_msy(fishery_sel = fishery_sel, M = M, catch_waa = waa_catch, ssb_waa = waa, paa = paa, ages = ages, Ninit = OM_report$equilibrium_at_age, plus_group = T, prop_Z = 0.5, n_runs = 200)

## Reset the OM parameters
OM_pars$ln_Fmax = log(F_max$F_max)
OM_pars$ln_F40 = log(F_40$F_ref)
OM_pars$ln_F35 = log(F_35$F_ref)
OM_pars$ln_F30 = log(F_30$F_ref)
OM_pars$ln_Fmsy = log(F_msy$F_msy)
OM_pars$ln_F_0_1 = log(F_0.1$F_0.1)

## re run the OM with true values so 
## we can check the log-likelihood calculations
OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_report = OM_obj$report()

## Check Negative Log-Likelihood
c(F_msy$mle_par$objective, OM_report$F_msy_nll)
c(F_max$mle_par$objective, OM_report$F_max_nll)
c(F_30$mle_par$objective, OM_report$F30_nll)
c(F_35$mle_par$objective, OM_report$F35_nll)
c(F_40$mle_par$objective, OM_report$F40_nll)

F_msy$SPR
OM_report$SPR_msy
get_spr(F_msy$F_msy, fishery_sel = fishery_sel, M = M, waa = waa, paa = paa, ages = ages)
F_msy$YPR
OM_report$YPR_msy


## Plot YPR with F and f-max
get_ypr(F_tilde = 0.001, fishery_sel, M, waa, ages)
trial_Fs = seq(from = 0.001, to = 1, by = 0.001)
temp = sapply(trial_Fs, FUN = get_ypr, fishery_sel, M, waa, ages)
plot(trial_Fs, temp, type = "l", lwd = 2, xlab = "F", ylab = "YPR")
abline(v = F_max$F_max, col = "red", lty =2, lwd = 2)
abline(v = F_0.1$F_0.1, col = "blue", lty =2, lwd = 2)
abline(v = F_msy$F_msy, col = "purple", lty =3, lwd = 2)
abline(v = F_30$F_ref, col = "gold", lty =2, lwd = 2)
abline(v = F_40$F_ref, col = "darkgreen", lty =2, lwd = 2)
legend('topright', legend = c("Fmax","Fmsy", "F0.1", "F30%spr", "F40%spr"), col = c("red", "purple", "blue", "gold", "darkgreen"), lty = 2, lwd =2)

## Use functions to get Fmax and F_0.1
F_0.1_alt$F_0.1
F_0.1$F_0.1
## Use functions to get Fmsy
F_msy$F_msy
F_msy$Bmsy / OM_report$B0 
F_msy$Bmsy
OM_report$B_msy

OM_report$FM

TMB_data$F_method = 0
OM_pars$ln_F = array(log(F_30$F_ref), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
TMB_data$estimate_F_init = 0
TMB_data$estimate_init_age_devs = 0
TMB_data$stockRecruitmentModelCode = 0
OM_pars$ln_R0 = log(14)
OM_pars$ln_ycs_est = rep(0, length(OM_pars$ln_ycs_est))
OM_obj_30 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)

OM_30_report = OM_obj_30$report()
OM_30_report$depletion * 100
OM_30_report$pred_catches

OM_30_report$B0
OM_30_report$Binit
OM_30_report$N[1,]
OM_30_report$ycs

OM_pars$ln_F = array(log(F_35$F_ref), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_obj_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_35_report = OM_obj_35$report()

OM_pars$ln_F = array(log(F_40$F_ref), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_obj_40 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_40_report = OM_obj_40$report()

plot(full_years, OM_30_report$depletion * 100, type = "l", xlab = "Year", ylab = "Depletion", lty = 1, lwd = 2, ylim = c(0,100))
lines(full_years, OM_35_report$depletion * 100, lty = 2, lwd = 2, col = "red")
lines(full_years, OM_40_report$depletion * 100, lty = 2, lwd = 2, col = "blue")
abline(h = 40, lty = 2,col = "blue")
abline(h = 35, lty = 2,col = "red")
abline(h = 30, lty = 2,col = "black")


## 80 years of F_50 then 80 years of F_25
OM_pars$ln_F = array(c(rep(log(F_50$F_ref), TMB_data$n_years /2),rep(log(F_20$F_ref), TMB_data$n_years /2)), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_obj_50_20 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_50_20_report = OM_obj_50_20$report()

plot(full_years, OM_50_20_report$depletion * 100, type = "l", xlab = "Year", ylab = "Depletion", lty = 1, lwd = 2, ylim = c(0,100))
abline(h = 50, lty = 2,col = "blue")
abline(h = 20, lty = 2,col = "red")

## now look at depletion based reference points with 
## stock recruit functions
F_30_dep = find_F_depletion(target_depletion = 30, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, SRmodel = 0, SR_pars = NULL)
F_30_dep_sr = find_F_depletion(target_depletion = 30, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, SRmodel = 1, SR_pars =TMB_data$steepness)
F_30_dep$F_ref
F_30$F_ref
## with SR
F_30_dep_sr$F_ref

TMB_data$stockRecruitmentModelCode = 2
TMB_data$steepness
OM_pars$ln_F = array(log(F_30_dep_sr$F_ref), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_obj_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_35_report = OM_obj_35$report()
tail(OM_35_report$depletion)


## compare our methods with fishmethods
YPR = fishmethods::ypr(incrF  = 0.0001, maxF = 1, age = ages, wgt = waa, partial = fishery_sel, M = rep(this_bio$M, length(ages)), plus = T, oldest = 400)

##
c(F_0.1$F_0.1, YPR$Reference_Points[1,1])
c(F_max$F_max, YPR$Reference_Points[2,1])
