#' validate afe-structured model for 
#' Medium biology
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

n_last_ycs_to_estimate = 5;
n_first_ycs_to_estimate = 1; ## you may want to shift this depending on when data starts recruits are observed


fig_dir = file.path(DIR$fig, "ValidateMediumBiology")
if(!dir.exists(fig_dir))
  dir.create(fig_dir)

this_bio = readRDS(file = file.path(DIR$data, "Medium_biology.RDS"))
n_years_historic = 20
n_years_data = 40
n_years = n_years_historic + n_years_data
years = (2020 - n_years + 1):2020
full_years = (min(years) - 1):max(years) # first year is unfished year
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
  ln_extra_survey_cv = log(0.0001),
  ln_F_init = log(0.01),
  ln_init_age_devs = rep(0, TMB_data$n_init_age_devs),
  ln_sigma_init_age_devs = log(0.6),
  logit_f_a50 = logit_general(rep(this_bio$f_a50, TMB_data$n_fisheries), TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  logit_f_ato95 = logit_general(rep(this_bio$f_ato95, TMB_data$n_fisheries), TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
  logit_survey_a50 = logit_general(this_bio$s_a50, TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  logit_survey_ato95 = logit_general(this_bio$s_ato95, TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
  logit_surveyQ = qlogis(0.2),
  ln_F = array(log(fishery_probs), dim = c(TMB_data$n_fisheries,TMB_data$n_years)),
  ln_catch_sd = log(0.002),
  ln_Fmax = log(0.05),
  ln_F40 = log(0.05),
  ln_F35 = log(0.05),
  ln_F30 = log(0.05),
  ln_Fmsy = log(0.06),
  ln_F_0_1 = log(0.04)
  
)

# these parameters we are not estimating.
na_map = fix_pars(par_list = OM_pars, pars_to_exclude = c("ln_catch_sd", "ln_extra_survey_cv","ln_sigma_r", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_map_alt = fix_pars(par_list = OM_pars, pars_to_exclude = c("ln_catch_sd","ln_F", "ln_extra_survey_cv","ln_sigma_r", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))

OM_obj <- MakeADFun(TMB_data, OM_pars, map = na_map, DLL= "AgeStructuredModel", checkParameterOrder = T)
EM_lst = EM_alt_lst = OM_lst = list()
OM_lst[["OM"]] = OM_obj$report()
n_sims = 100
for(j in 1:n_sims) {
  simdata = OM_obj$simulate(complete = T)
  #simdata$F_method = 1
  EM <- MakeADFun(simdata, OM_pars, map = na_map, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
  mle_EM1 = nlminb(start = EM$par, objective = EM$fn, gradient  = EM$gr, control = list(iter.max = 10000, eval.max = 10000))
  try_improve = tryCatch(expr =
                           for(i in 1:2) {
                             g = as.numeric(EM$gr(mle_EM1$par))
                             h = optimHess(mle_EM1$par, fn = EM$fn, gr = EM$gr)
                             mle_EM1$par = mle_EM1$par - solve(h,g)
                             mle_EM1$objective = EM$fn(mle_EM1$par)
                           }
                         , error = function(e){e})
  if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
    cat("Failed simulation EM, sim ", j, "\n")
  }
  EM_rep = EM$report(mle_EM1$par)
  EM_lst[[as.character(j)]] = EM_rep

  ## change F method
  simdata$F_method = 1
  EM_alt <- MakeADFun(simdata, OM_pars, map = na_map_alt, DLL= "AgeStructuredModel", checkParameterOrder = T, silent = T)
  mle_EM1_alt = nlminb(start = EM_alt$par, objective = EM_alt$fn, gradient  = EM_alt$gr, control = list(iter.max = 10000, eval.max = 10000))
  try_improve = tryCatch(expr =
                           for(i in 1:2) {
                             g = as.numeric(EM_alt$gr(mle_EM1_alt$par))
                             h = optimHess(mle_EM1_alt$par, fn = EM_alt$fn, gr = EM_alt$gr)
                             mle_EM1_alt$par = mle_EM1_alt$par - solve(h,g)
                             mle_EM1_alt$objective = EM_alt$fn(mle_EM1_alt$par)
                           }
                         , error = function(e){e})
  if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
    cat("Failed simulation EM_alt, sim ", j, "\n")
  }
  EM_rep_alt = EM_alt$report(mle_EM1_alt$par)
  EM_alt_lst[[as.character(j)]] = EM_rep_alt
  
}

## Get a range of reference points
full_ssb = NULL
OM_ssb = data.frame(names = full_years, values = OM_lst[["OM"]]$ssb)
## EM1
EM1_ssb_df = get_multiple_vectors(EM_lst, "ssb", element_labs = full_years)
EM1_ssb_df$model = "EM"
#EM1_ssb_df$RE = (EM1_ssb_df$values - OM_lst[["OM"]]$ssb) / OM_lst[["OM"]]$ssb * 100
full_ssb = rbind(full_ssb, EM1_ssb_df)

full_ssb = full_ssb %>% left_join(OM_ssb, by = c("names"))
full_ssb = full_ssb %>% group_by(names) %>% mutate(RE = (values.x - values.y)/ values.y * 100)

ggplot(full_ssb, aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ylim(-100,100)

OM_dep = data.frame(names = full_years, values = OM_lst[["OM"]]$depletion)
EM1_dep_df = get_multiple_vectors(EM_lst, "depletion", element_labs = full_years)
EM1_dep_df$model = "EM"
EM1_dep_df = EM1_dep_df %>% left_join(OM_dep, by = c("names"))
EM1_dep_df = EM1_dep_df %>% group_by(names) %>% mutate(RE = (values.x - values.y)/ values.y * 100)

ggplot(EM1_dep_df, aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ylim(-100,100)

