#'
#' Show all the F's for the medium biology run
#' For all the OM factors
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
compile(file = file.path(DIR$tmb, "AgeStructuredModelOM.cpp"), flags = "-Wignored-attributes -O3")

#sink()
#dyn.unload(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))
#dyn.unload(dynlib(file.path(DIR$tmb, "AgeStructuredModelOM")))

dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModel")))
dyn.load(dynlib(file.path(DIR$tmb, "AgeStructuredModelOM")))

#setwd(DIR$R)

fig_dir = file.path(DIR$fig, "MediumBiology")
if(!dir.exists(fig_dir))
  dir.create(fig_dir)

this_bio = readRDS(file = file.path(DIR$data, "Medium_biology.RDS"))

## data period 
n_years_historic = 20
n_years_data = 60

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

fishery_probs = c(seq(from = 0.01, to = this_bio$M * 1.5, length = 20), rep(this_bio$M, 40), rep(this_bio$M, 20) )
plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)

## The same parameters as OM, to check for consistency
OM_pars = list(
  ln_R0 = log(this_bio$R0),
  ln_ycs_est =  rep(0, sum(TMB_data$ycs_estimated)),
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
fishery_probs = c(seq(from = 0.01, to = F_15$F_ref, length = 30), seq(from = F_25$F_ref, to = F_40$F_ref, length = 10), rep(F_40$F_ref, 40))
plot(years, fishery_probs, type = "l", lty = 2, lwd = 2)

TMB_data$fixed_catch_indicator = rep(0, n_years)
#TMB_data$fixed_catch_indicator[40:n_years] = 1
OM_pars$ln_F = array(log(fishery_probs), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_obj <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_report = OM_obj$report()

## n-years during historic period for each initial
# 20 for 100
# 20 for 75
# 20 for 50
# 20 for 25

## F-rebuild
## With a 15 year up between historical and data period
ramp = F
if(ramp) {
  f_100 = 0.0000001
  fish_100_50 = c(rep(f_100,n_years_historic), seq(from = f_100, to = F_50$F_ref, length = 25), rep(F_50$F_ref, 35))
  fish_100_35 = c(rep(f_100,n_years_historic), seq(from = f_100, to = F_35$F_ref, length = 25), rep(F_35$F_ref, 35))
  fish_100_20 = c(rep(f_100,n_years_historic),  seq(from = f_100, to = F_20$F_ref, length = 25), rep(F_20$F_ref, 35))
  
  fish_75_50 = c(rep(init_F_75$F_ref,n_years_historic), seq(from = init_F_75$F_ref, to = F_50$F_ref, length = 25), rep(F_50$F_ref, 35))
  fish_75_35 = c(rep(init_F_75$F_ref,n_years_historic), seq(from = init_F_75$F_ref, to = F_35$F_ref, length = 25), rep(F_35$F_ref, 35))
  fish_75_20 = c(rep(init_F_75$F_ref,n_years_historic), seq(from = init_F_75$F_ref, to = F_20$F_ref, length = 25), rep(F_20$F_ref, 35))
  
  fish_50_50 = c(rep(init_F_50$F_ref,n_years_historic), seq(from = init_F_50$F_ref, to = F_50$F_ref, length = 25), rep(F_50$F_ref, 35))
  fish_50_35 = c(rep(init_F_50$F_ref,n_years_historic), seq(from = init_F_50$F_ref, to = F_35$F_ref, length = 25), rep(F_35$F_ref, 35))
  fish_50_20 = c(rep(init_F_50$F_ref,n_years_historic), seq(from = init_F_50$F_ref, to = F_20$F_ref, length = 25), rep(F_20$F_ref, 35))
  
  fish_25_50 = c(rep(init_F_25$F_ref,n_years_historic), seq(from = init_F_25$F_ref, to = F_50$F_ref, length = 25), rep(F_50$F_ref, 35))
  fish_25_35 = c(rep(init_F_25$F_ref,n_years_historic), seq(from = init_F_25$F_ref, to = F_35$F_ref, length = 25), rep(F_35$F_ref, 35))
  fish_25_20 = c(rep(init_F_25$F_ref,n_years_historic), seq(from = init_F_25$F_ref, to = F_20$F_ref, length = 25), rep(F_20$F_ref, 35))
} else {
  f_100 = 0.0000001
  fish_100_50 = c(rep(f_100,n_years_historic), rep(F_50$F_ref, 60))
  fish_100_35 = c(rep(f_100,n_years_historic), rep(F_35$F_ref, 60))
  fish_100_20 = c(rep(f_100,n_years_historic), rep(F_20$F_ref, 60))
  
  fish_75_50 = c(rep(init_F_75$F_ref,n_years_historic), rep(F_50$F_ref, 60))
  fish_75_35 = c(rep(init_F_75$F_ref,n_years_historic), rep(F_35$F_ref, 60))
  fish_75_20 = c(rep(init_F_75$F_ref,n_years_historic), rep(F_20$F_ref, 60))
  
  fish_50_50 = c(rep(init_F_50$F_ref,n_years_historic), rep(F_50$F_ref, 60))
  fish_50_35 = c(rep(init_F_50$F_ref,n_years_historic), rep(F_35$F_ref, 60))
  fish_50_20 = c(rep(init_F_50$F_ref,n_years_historic), rep(F_20$F_ref, 60))
  
  fish_25_50 = c(rep(init_F_25$F_ref,n_years_historic), rep(F_50$F_ref, 60))
  fish_25_35 = c(rep(init_F_25$F_ref,n_years_historic), rep(F_35$F_ref, 60))
  fish_25_20 = c(rep(init_F_25$F_ref,n_years_historic), rep(F_20$F_ref, 60))
}

## plot Fs
png(filename = file.path(fig_dir, "OM_Fstrategys.png"), units = "in", width = 9, height = 8, res = 250)
par(mfrow = c(2,2), mar = c(4,4.5,2,1), cex.lab = 1.2, cex.axis = 1.2, oma = c(1,1,0,0))
plot(years, fish_100_50, type = "l", lty = 2, lwd = 2, ylim = c(0,0.13), ylab = "", xlab = "", main = "100% initial depletion")
lines(years, fish_100_35, col = "red", lty = 2, lwd = 2)
lines(years, fish_100_20, col = "blue", lty = 2, lwd = 2)
legend('topleft', col = c("black", "red","blue"), legend =c("F 50%", "F 35%", "F 20%"), lty = 2, lwd = 3, title = "Data period F")

plot(years, fish_75_50, type = "l", lty = 2, lwd = 2, ylim = c(0,0.13), ylab = "", xlab = "Years", main = "75% initial depletion")
lines(years, fish_75_35, col = "red", lty = 2, lwd = 2)
lines(years, fish_75_20, col = "blue", lty = 2, lwd = 2)

plot(years, fish_50_50, type = "l", lty = 2, lwd = 2, ylim = c(0,0.13), ylab = "", xlab = "", main = "50% initial depletion")
lines(years, fish_50_35, col = "red", lty = 2, lwd = 2)
lines(years, fish_50_20, col = "blue", lty = 2, lwd = 2)

plot(years, fish_25_50, type = "l", lty = 2, lwd = 2, ylim = c(0,0.13), ylab = "", xlab = "", main = "25% initial depletion")
lines(years, fish_25_35, col = "red", lty = 2, lwd = 2)
lines(years, fish_25_20, col = "blue", lty = 2, lwd = 2)
mtext(text = "Fully selected F",side = 2, outer = T, at = 0.5, line = -1, cex= 1.2)
mtext(text = "Years",side = 1, outer = T, at = 0.5, line = -1, cex= 1.2)
dev.off()

#########
## Now what the SSB trajectories look like with 
## deterministic recruitment for these
########
OM_pars$ln_F = array(log(fish_100_50), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_100_50 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_100_50_rep = OM_100_50$report()

OM_pars$ln_F = array(log(fish_100_35), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_100_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_100_35_rep = OM_100_35$report()

OM_pars$ln_F = array(log(fish_100_20), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_100_20 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_100_20_rep = OM_100_20$report()

OM_pars$ln_F = array(log(fish_75_50), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_75_50 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_75_50_rep = OM_75_50$report()

OM_pars$ln_F = array(log(fish_75_35), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_75_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_75_35_rep = OM_75_35$report()

OM_pars$ln_F = array(log(fish_75_20), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_75_20 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_75_20_rep = OM_75_20$report()

OM_pars$ln_F = array(log(fish_50_50), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_50_50 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_50_50_rep = OM_50_50$report()

OM_pars$ln_F = array(log(fish_50_35), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_50_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_50_35_rep = OM_50_35$report()

OM_pars$ln_F = array(log(fish_50_20), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_50_20 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_50_20_rep = OM_50_20$report()

OM_pars$ln_F = array(log(fish_25_50), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_25_50 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_25_50_rep = OM_25_50$report()

OM_pars$ln_F = array(log(fish_25_35), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_25_35 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_25_35_rep = OM_25_35$report()

OM_pars$ln_F = array(log(fish_25_20), dim = c(TMB_data$n_fisheries,TMB_data$n_years))
OM_25_20 <- MakeADFun(TMB_data, OM_pars, DLL= "AgeStructuredModel", checkParameterOrder = T)
OM_25_20_rep = OM_25_20$report()

#########
## plot SSB's
#########
## 
plot_year1 = min(full_years) - 10
plot_year2 = min(full_years) + n_years_historic

png(filename = file.path(fig_dir, "OM_SSBs.png"), units = "in", width = 9, height = 8, res = 250)
par(mfrow = c(2,2), mar = c(4,4.5,2,1), cex.lab = 1.2, cex.axis = 1.2, oma = c(1,1,0,0))
plot(full_years, OM_100_50_rep$depletion * 100, type = "l", lty = 2, lwd = 2, ylim = c(0,100), ylab = "", xlab = "", main = "100% initial depletion")
lines(full_years, OM_100_35_rep$depletion * 100, col = "red", lty = 2, lwd = 2)
lines(full_years, OM_100_20_rep$depletion * 100, col = "blue", lty = 2, lwd = 2)
polygon(x = c(plot_year1,plot_year2,plot_year2,plot_year1), y = c(120, 120, -10, -10), fill = "gray60", density = 20)
legend('topright', col = c("black", "red","blue"), legend =c("F 50%", "F 35%", "F 20%"), lty = 2, lwd = 3, title = "Data period F")

plot(full_years, OM_75_50_rep$depletion * 100, type = "l", lty = 2, lwd = 2, ylim = c(0,100), ylab = "", xlab = "Years", main = "75% initial depletion")
lines(full_years, OM_75_35_rep$depletion * 100, col = "red", lty = 2, lwd = 2)
lines(full_years, OM_75_20_rep$depletion * 100, col = "blue", lty = 2, lwd = 2)
polygon(x = c(plot_year1,plot_year2,plot_year2,plot_year1), y = c(120, 120, -10, -10), fill = "gray60", density = 20)

plot(full_years, OM_50_50_rep$depletion * 100, type = "l", lty = 2, lwd = 2, ylim = c(0,100), ylab = "", xlab = "", main = "50% initial depletion")
lines(full_years, OM_50_35_rep$depletion * 100, col = "red", lty = 2, lwd = 2)
lines(full_years, OM_50_20_rep$depletion * 100, col = "blue", lty = 2, lwd = 2)
polygon(x = c(plot_year1,plot_year2,plot_year2,plot_year1), y = c(120, 120, -10, -10), fill = "gray60", density = 20)

plot(full_years, OM_25_50_rep$depletion * 100, type = "l", lty = 2, lwd = 2, ylim = c(0,100), ylab = "", xlab = "", main = "25% initial depletion")
lines(full_years, OM_25_35_rep$depletion * 100, col = "red", lty = 2, lwd = 2)
lines(full_years, OM_25_20_rep$depletion * 100, col = "blue", lty = 2, lwd = 2)
polygon(x = c(plot_year1,plot_year2,plot_year2,plot_year1), y = c(120, 120, -10, -10), fill = "gray60", density = 20)

mtext(text = "Fully selected F",side = 2, outer = T, at = 0.5, line = -1, cex= 1.2)
mtext(text = "Years",side = 1, outer = T, at = 0.5, line = -1, cex= 1.2)
dev.off()






