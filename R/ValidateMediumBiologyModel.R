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
TMB_data$survey_obs = rep(1, sum(TMB_data$survey_year_indicator))
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
TMB_data$ycs_bias_correction = rep(1, n_years)
## don't apply bias correction to the last 5 years because of lack of data
TMB_data$ycs_bias_correction[(n_years - 4):n_years] = 0
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
  ln_catch_sd = log(0.002),
  ln_Fmax = log(0.05),
  ln_F40 = log(0.05),
  ln_F35 = log(0.05),
  ln_F30 = log(0.05),
  ln_Fmsy = log(0.06),
  ln_F_0_1 = log(0.04)
  
)

# these parameters we are not estimating.
na_map = fix_pars(par_list = OM_pars, pars_to_exclude = c("ln_catch_sd", "ln_sigma_r", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))
na_map_alt = fix_pars(par_list = OM_pars, pars_to_exclude = c("ln_catch_sd","ln_F", "ln_sigma_r", "ln_F_init", "ln_init_age_devs", "ln_sigma_init_age_devs", "ln_Fmax", "ln_F40", "ln_F35", "ln_F30", "ln_Fmsy", "ln_F_0_1"))

OM_obj <- MakeADFun(TMB_data, OM_pars, map = na_map, DLL= "AgeStructuredModel", checkParameterOrder = T)
EM_lst = EM_alt_lst = OM_lst = list()
OM_lst[["OM"]] = OM_obj$report()
EM_est_time = EM_alt_est_time = vector()
n_sims = 100
for(j in 1:n_sims) {
  if(j %% 50 == 0)
    cat("j = ", j, "\n");
  simdata = OM_obj$simulate(complete = T)
  #simdata$F_method = 1
  start_time = Sys.time()
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
  EM_est_time[j] = difftime(Sys.time(), start_time, units = "secs")
  ## change F method
  simdata$F_method = 1
  start_time = Sys.time()
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
  EM_alt_est_time[j] = difftime(Sys.time(), start_time, units = "secs")
  
}

## Get a range of reference points
full_ssb = NULL
OM_ssb = data.frame(names = full_years, values = OM_lst[["OM"]]$ssb)
## EM1
EM1_ssb_df = get_multiple_vectors(EM_lst, "ssb", element_labs = full_years)
EM1_ssb_df$model = "EM"
EM1_alt_ssb_df = get_multiple_vectors(EM_alt_lst, "ssb", element_labs = full_years)
EM1_alt_ssb_df$model = "EM-alt"
#EM1_ssb_df$RE = (EM1_ssb_df$values - OM_lst[["OM"]]$ssb) / OM_lst[["OM"]]$ssb * 100
full_ssb = rbind(EM1_ssb_df, EM1_alt_ssb_df)

full_ssb = full_ssb %>% left_join(OM_ssb, by = c("names"))
full_ssb = full_ssb %>% group_by(names) %>% mutate(RE = (values.x - values.y)/ values.y * 100)

ggplot(full_ssb, aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ylim(-100,100) +
  facet_wrap(~model)

OM_dep = data.frame(names = full_years, values = OM_lst[["OM"]]$depletion)
EM1_dep_df = get_multiple_vectors(EM_lst, "depletion", element_labs = full_years)
EM1_dep_df$model = "EM"
EM1_alt_dep_df = get_multiple_vectors(EM_alt_lst, "depletion", element_labs = full_years)
EM1_alt_dep_df$model = "EM-alt"
full_dep = rbind(EM1_dep_df, EM1_alt_dep_df)

full_dep = full_dep %>% left_join(OM_dep, by = c("names"))
full_dep = full_dep %>% group_by(names) %>% mutate(RE = (values.x - values.y)/ values.y * 100)

ggplot(full_dep, aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ylim(-100,100)+
  facet_wrap(~model)

## recruitment
OM_ycs = data.frame(names = years, values = OM_lst[["OM"]]$ycs)
EM1_ycs_df = get_multiple_vectors(EM_lst, "ycs", element_labs = years)
EM1_ycs_df$model = "EM"
EM1_alt_ycs_df = get_multiple_vectors(EM_alt_lst, "ycs", element_labs = years)
EM1_alt_ycs_df$model = "EM-alt"
full_ycs = rbind(EM1_ycs_df, EM1_alt_ycs_df)

full_ycs = full_ycs %>% left_join(OM_ycs, by = c("names"))
full_ycs = full_ycs %>% group_by(names) %>% mutate(RE = (values.x - values.y)/ values.y * 100)
ggplot(full_ycs, aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  ylim(-100,100)+
  facet_wrap(~model)


tail(OM_ycs)
full_ycs %>% filter(names %in% c(2016:2020)) %>% group_by(names, model) %>% summarise(mean(values.x))

ggplot(full_ycs, aes(x = names, group = names)) +
  geom_boxplot(aes(y = values.x)) + 
  geom_hline(yintercept = 1, col = "red", linewidth = 1.1, linetype = "dashed") +
  geom_point(aes(y = values.y), col = "steelblue", size= 1) +
  theme_bw() +
  facet_wrap(~model)

## Selectivities



## look at survey index fit
OM_index_df = data.frame(names = years, values = OM_lst[["OM"]]$survey_index_fitted)
EM1_index_df = get_multiple_vectors(EM_lst, "survey_index_fitted", element_labs = survey_year_obs)
EM1_obs_index_df = get_multiple_vectors(EM_lst, "survey_obs", element_labs = survey_year_obs)
EM1_index_df$observed = EM1_obs_index_df$value
EM1_index_df$model = "EM"
EM_alt_index_df = get_multiple_vectors(EM_alt_lst, "survey_index_fitted", element_labs = survey_year_obs)
EM_alt_obs_index_df = get_multiple_vectors(EM_alt_lst, "survey_obs", element_labs = survey_year_obs)
EM_alt_index_df$observed = EM_alt_obs_index_df$value
EM_alt_index_df$model = "EM-alt"
EM_index = rbind(EM1_index_df, EM_alt_index_df)
## match OM_catch to full catch
EM_index = EM_index %>% left_join(OM_index_df, by = c("names"))
## calculate RE
EM_index = EM_index %>% group_by(names) %>% mutate(RE = (values.x - values.y)/ values.y * 100)

ggplot(EM_index, aes(x = names, group = names)) +
  geom_boxplot(aes(y = RE)) + 
  geom_hline(yintercept = 0, col = "red", linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  facet_wrap(~model) +
  ylim(-15,15)

new_dat = EM_index %>% filter(model == "EM") %>% group_by(names) %>% summarise(avg_EM = mean(values.x), avg_OM = mean(values.y), sim_OM = mean(observed))
ggplot(new_dat, aes(x = names)) +
  geom_line(aes(y = avg_EM, col = "avg_EM", linetype = "avg_EM"), linewidth = 1) +
  geom_line(aes(y = avg_OM, col = "avg_OM", linetype = "avg_OM"), linewidth = 1) +
  geom_line(aes(y = sim_OM, col = "sim_OM", linetype = "sim_OM"), linewidth = 1) +
  theme_bw() +
  xlim(2000, NA) +
  ylim(50000, 100000)

## survey AF
obs = get("survey_AF_obs", OM_lst[["OM"]])
fit = get("survey_AF_fitted", OM_lst[["OM"]])
colnames(obs) = colnames(fit) = survey_year_obs
molten_obs = reshape2::melt(obs, varnames = c("Age", "year"), value.name = "obs")
molten_fit = reshape2::melt(fit, varnames = c("Age", "year"), value.name = "fit")
molten_obs$fit = molten_fit$fit

OM_AF = molten_obs
EM_AF_df = get_age_fit_by_year(EM_lst, survey = T,  element_labs = survey_year_obs, years = survey_year_obs)
EM_AF_df$model = "EM"
EM_alt_AF_df = get_age_fit_by_year(EM_alt_lst, survey = T,  element_labs = survey_year_obs, years = survey_year_obs)
EM_alt_AF_df$model = "EM-alt"
full_survey_AF = rbind(EM_AF_df, EM_alt_AF_df)
## match OM_catch to full catch
full_survey_AF = full_survey_AF %>% left_join(OM_AF, by = c("Age", "year"))
## calculate RE
full_survey_AF = full_survey_AF %>% group_by(Age, year) %>% mutate(RE = (fit.x - fit.y)/ fit.y * 100)

ggplot(full_survey_AF %>% filter(year %in% c(2018:2020)), aes(x = Age, group = Age, y = RE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
  theme_bw() +
  ggtitle("EM1") +
  ylim(-50, 50) +
  facet_grid(year~model) 

full_survey_AF %>% filter(sim_iter == 1, Age == 1, year == 2020)


