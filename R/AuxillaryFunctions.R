#'
#'
#' Auxilary functions
#'
#'
library(numDeriv)

DIR = list()
DIR$base = getwd()
DIR$data = file.path("..", "Data")
DIR$fig = file.path("..", "Figures")
DIR$tmb = file.path("..", "TMB")
#'
#' get_spr 
#' get spawner per recruit for a given F (F_tilde)
#' @param F_tilde annual fishing mortality
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param paa vector of proportions at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @return spawner per recruit
get_spr <- function(F_tilde, fishery_sel, M, waa, paa, ages, prop_Z = 0.5) {
  Za = F_tilde * fishery_sel + M
  Na = 1;
  SPR = Na * exp(-Za[1] * prop_Z) * waa[1] * paa[1]
  for(age_iter in 2:(length(ages)*4)) {
    if(age_iter > length(ages)) {
      SPR = SPR + Na * exp(-Za[length(ages)] * prop_Z) * waa[length(ages)] * paa[length(ages)]
      Na = Na * exp(-Za[length(ages)]);
    } else {
      SPR = SPR + Na * exp(-Za[age_iter] * prop_Z) * waa[age_iter] * paa[age_iter]
      Na = Na * exp(-Za[age_iter]);
    }
    
  }
  return(SPR)
}
#'
#' get_SSBe
#' get equilbrium SSB 
#' @param F_tilde annual fishing mortality
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param paa vector of proportions at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @param plus_group boolean
#' @param prop_Z proportion of mortality for Spawning calculation
#' 
#' @return spawner per recruit
get_SSBe <- function(F_tilde, fishery_sel, M, waa, paa, ages, Ninit, plus_group = T, prop_Z = 0.5, n_runs = 100) {
  Za = F_tilde * fishery_sel + M
  Na = 1;
  SPR = Na * exp(-Za[1] * prop_Z) * waa[1] * paa[1]
  N_equilibrium = array(0, dim = c(length(ages), n_runs))
  N_equilibrium[,1] = Ninit;
  for(year_ndx in 2:n_runs) {
    N_equilibrium[1, year_ndx] = Ninit[1]
    N_equilibrium[2:length(ages), year_ndx] = N_equilibrium[1:(length(ages) - 1), year_ndx - 1] * exp(-Za[1:(length(ages) - 1)])
    if(plus_group)
      N_equilibrium[length(ages), year_ndx] = N_equilibrium[length(ages) - 1, year_ndx - 1] * exp(-Za[length(ages) - 1]) + N_equilibrium[length(ages), year_ndx - 1] * exp(-Za[length(ages)]) 
  }
  SSBe = sum(as.numeric(N_equilibrium[,n_runs]) * exp(-Za * prop_Z) * waa * paa)
  SSBe_previous = sum(as.numeric(N_equilibrium[,n_runs - 1]) * exp(-Za * prop_Z) * waa * paa)
  if(abs(SSBe - SSBe_previous) > (SSBe * 0.01)) {
    warning("SSB equilibrium calculation may need to run for longer than input n_runs. We identified a difference between the last and second to last SSB values greater than 1%")
  }
  return(SSBe)
}
#' get_ypr 
#' get Yeild per recruit for a given F (F_tilde)
#' @param F_tilde annual fishing mortality
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @return spawner per recruit
get_ypr <- function(F_tilde, fishery_sel, M, waa, ages) {
  Za = F_tilde * fishery_sel + M
  Na = 1;
  YPR = (F_tilde * fishery_sel[1]) / Za[1] * (1 - exp(-Za[1])) * Na * waa[1]
  YPR_vec = rep(0, length(ages)*4)
  YPR_vec[1] = YPR
  for(age_iter in 2:(length(ages)*4)) {
    if(age_iter > length(ages)) {
      Na = Na * exp(-Za[length(ages)]);
      YPR = YPR + (F_tilde * fishery_sel[length(ages)]) / Za[length(ages)] * (1 - exp(-Za[length(ages)])) * Na * waa[length(ages)]
      YPR_vec[age_iter] = (F_tilde * fishery_sel[length(ages)]) / Za[length(ages)] * (1 - exp(-Za[length(ages)])) * Na * waa[length(ages)]
    } else {
      Na = Na * exp(-Za[age_iter]);
      YPR = YPR + (F_tilde * fishery_sel[age_iter]) / Za[age_iter] * (1 - exp(-Za[age_iter])) * Na * waa[age_iter]
      YPR_vec[age_iter] = (F_tilde * fishery_sel[age_iter]) / Za[age_iter] * (1 - exp(-Za[age_iter])) * Na * waa[age_iter]
    }
  }
  return(YPR)
}
#' dYPR 
#' get the derivative of the YPR for a given F_tilde
#' @param F_tilde annual fishing mortality
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @return derivative for the Yeild per recruit

dYPR <-function(F_tilde, fishery_sel, M, waa, ages) {
  h = 0.001;
  ln_F_tilde = log(F_tilde)
  v = -get_ypr(exp(ln_F_tilde + 2.0 * h), fishery_sel, M, waa, ages) + 8.0 * get_ypr(exp(ln_F_tilde + h), fishery_sel, M, waa, ages) - 8.0 * get_ypr(exp(ln_F_tilde - h), fishery_sel, M, waa, ages) + get_ypr(exp(ln_F_tilde - 2.0 * h), fishery_sel, M, waa, ages)
  g = v / (12.0 * h);
  return(g / F_tilde);
}
#'
#' find_F_percent 
#' Find a F that achieves some percent B0
#' @param target_spr percentage SSB/B0 (depletion) that we want to find an F for
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param paa vector of proportions at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @param F_range vector of two values specifying the range of F's to optimise over.
#' @return F_% depletion

find_F_percent <- function(target_spr, fishery_sel, M, waa, paa, ages, prop_Z = 0.5, F_range = c(0.00001, 1.5)) {
  nll_spr <- function(ln_F_tilde, target_spr, fishery_sel, M, waa, paa, ages, prop_Z) {
    SPR_F0 = get_spr(0, fishery_sel, M, waa, paa, ages, prop_Z);
    SPR_Ftilde = get_spr(exp(ln_F_tilde), fishery_sel, M, waa, paa, ages, prop_Z);
    return(((SPR_F0/100 * target_spr) - SPR_Ftilde)^2)
  }
  
  # optimes
  mle_par = optimise(interval = c(log(F_range[1]), log(F_range[2])), f = nll_spr, target_spr = target_spr, fishery_sel = fishery_sel, M = M, waa = waa, paa = paa, ages = ages, prop_Z = prop_Z)
  F_ref = exp(mle_par$minimum)
  SPR = get_spr(F_ref, fishery_sel, M, waa, paa, ages, prop_Z);
  SPR_F0 = get_spr(0, fishery_sel, M, waa, paa, ages, prop_Z);
  
  return(list(F_ref = F_ref, target_spr = target_spr, SPR = SPR, SPR_F0 = SPR_F0, percent_SPR = SPR / SPR_F0 * 100, mle_par = mle_par))
}

#'
#' find_F_max 
#' Find a F that maximises YPR
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @param F_range vector of two values specifying the range of F's to optimise over.
#' @return F_max 
find_F_max <- function(fishery_sel, M, waa, ages, F_range = c(0.00001, 1.5)) {
  nll_ypr <- function(ln_F_tilde, fishery_sel, M, waa, paa, ages) {
    YPR = get_ypr(exp(ln_F_tilde), fishery_sel, M, waa, ages);
    return(-1.0 * log(YPR))
  }
  
  # optimes
  mle_par = optimise(interval = c(log(F_range[1]), log(F_range[2])), f = nll_ypr, fishery_sel = fishery_sel, M = M, waa = waa, ages = ages)
  F_max = exp(mle_par$minimum)
  YPR = get_ypr(F_max, fishery_sel, M, waa, ages);
  
  return(list(F_max = F_max, YPR = YPR, mle_par = mle_par))
}
#'
#' find_F_msy 
#' Find a F that achieves MSY
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param catch_waa vector of weight at age for catch/Yield calculation
#' @param ssb_waa vector of weight at age for catch/Yield calculation
#' @param paa vector of proportions mature by age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @param Ninit vector of initial numbers at age
#' @param plus_group boolean
#' @param prop_Z proportion of Z mortality to take account of before calculating SSB
#' @param n_runs number of annual cycles to calculate equilibrium SSB for a given F
#' @param F_range vector of two values specifying the range of F's to optimise over.
#' @return F_msy  and other quantities
find_F_msy <- function(fishery_sel, M, catch_waa, ssb_waa, paa, ages, Ninit, plus_group = T, prop_Z = 0.5, n_runs = 100, F_range = c(0.00001, 1.5)) {
  
  nll_fmsy <- function(ln_F_tilde, fishery_sel, M, catch_waa, ssb_waa, paa, ages, Ninit, plus_group, prop_Z, n_runs) {
    this_F = exp(ln_F_tilde)
    YPR = get_ypr(this_F, fishery_sel, M, catch_waa, ages);
    SSBe = get_SSBe(this_F, fishery_sel, M, ssb_waa, paa, ages, Ninit, plus_group, prop_Z, n_runs)
    SPR = get_spr(this_F, fishery_sel, M, ssb_waa, paa, ages, prop_Z)
    return(-1.0 * log((YPR * SSBe)/SPR))
  }
  
  # optimes
  mle_par = optimise(interval = c(log(F_range[1]), log(F_range[2])), f = nll_fmsy, Ninit = Ninit, plus_group = plus_group, prop_Z = prop_Z, n_runs = n_runs, fishery_sel = fishery_sel, M = M, catch_waa = catch_waa, ssb_waa = ssb_waa, paa = paa, ages = ages)
  F_msy = exp(mle_par$minimum)
  YPR = get_ypr(F_msy, fishery_sel, M, catch_waa, ages);
  Bmsy = get_SSBe(F_msy, fishery_sel, M, ssb_waa, paa, ages, Ninit, plus_group, prop_Z, n_runs);
  SPR = get_spr(F_msy, fishery_sel, M, ssb_waa, paa, ages, prop_Z)
  
  return(list(F_msy = F_msy, Bmsy = Bmsy, SPR = SPR, YPR = YPR, mle_par = mle_par))
}

#'
#' find_F_0.1 
#' Find a F_0.1 which is the F that has a derivative of 0.1 on the YPR curve
#' @param fishery_sel vector of fishing selectivities
#' @param M vector of natural mortality
#' @param waa vector of weight at age
#' @param ages vector of ages - the last age is assumed to be a plus group
#' @param derivative_method 1 = use our dYPR function which is a numerical difference evaluation, or use the dependent thirdparty numDeriv::grad (=2)
#' @param F_range vector of two values specifying the range of F's to optimise over.
#' @detaails the MLE can be sensitive to upper values of F_range
#' @importFrom numDeriv grad
#' @return F_0.1 and other information 
find_F_0.1 <- function(fishery_sel, M, waa, ages, derivative_method = 1, F_range = c(0.00001, 1.5)) {
  nll_ypr <-  
  if(derivative_method == 1) {
    function(ln_F_tilde, fishery_sel, M, waa, paa, ages) {
      dYPR_F0 = dYPR(0.000000001, fishery_sel, M, waa, ages); ## almost zero
      dYPR_F_tilde = dYPR(exp(ln_F_tilde), fishery_sel, M, waa, ages);
      return((0.1 * dYPR_F0 - dYPR_F_tilde)^2)
    }
  } else {
    function(ln_F_tilde, fishery_sel, M, waa, paa, ages) {
      dYPR_F0 = numDeriv::grad(func = get_ypr, x = 0.000000001, fishery_sel = fishery_sel, M = M, waa = waa, ages = ages)
      dYPR_F_tilde = numDeriv::grad(func = get_ypr, x = exp(ln_F_tilde), fishery_sel = fishery_sel, M = M, waa = waa, ages = ages)
      return((0.1 * dYPR_F0 - dYPR_F_tilde)^2)
    }
  }
  # optimes
  mle_par = optimise(interval = c(log(F_range[1]), log(F_range[2])), f = nll_ypr, fishery_sel = fishery_sel, M = M, waa = waa, ages = ages)
  F_0.1 = exp(mle_par$minimum)
  YPR = get_ypr(F_0.1, fishery_sel, M, waa, ages);
  return(list(F_0.1 = F_0.1, YPR = YPR, mle_par = mle_par))
  
}


#' vonbert applies the Von Bertalanffy age-length relationship
#' @param age take an age look in globe scope for the rest of parameters
#' @param L_inf asympototic length
#' @param K growth rate parameter
#' @param t0 age for length = 0
#' @return mean length at age for the Von Bertalanffy growth function
vonbert <- function(age,K,L_inf,t0) {
  return(L_inf * (1-exp(-K*(age -t0))))
}

#'  Logistic selectivity
#' @param x age or length to evaluate selectivity at
#' @param x50 value of x where selectivity is equal to 0.5
#' @param a95 difference in x when the selectivity is equal to 0.5 and 0.95
#' @return selectivity along x
#' @export
logis_sel<- function(X,a50,a95)
{
  1/(1+19^((a50-X)/a95)) 
}

#' Logistic selectivity reparameterised for a50 and slope
#' @param x age or length to evaluate selectivity at
#' @param x50 value of x where selectivity is equal to 0.5
#' @param slope yearly slope
#' @return selectivity along x
#' @export
logistic_sel_alternative = function (X, a50, slope) {
  1/(1 + exp(slope * (X - a50)))
}

#' logit_general bounds X which is between [lb,ub] to -inf -> inf based on the logit transformation
#' @param X scalar range [lb,ub]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @export
#' @return Y to be between [-inf, inf]
logit_general = function(X, lb, ub) {
  X1 = (X - lb) / (ub - lb)
  log(X1/(1 - X1))
}
#' invlogit_general bounds X which is between -inf -> inf to [lb,ub] based on the logit transformation
#' @param Y scalar range [-inf, inf]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @export
#' @return X to be between [lb,ub]
invlogit_general = function(Y, lb, ub) {
  Y1 = 1 / (1 + exp(-Y))
  lb + (ub - lb)*Y1
}


#' fix_pars
#' @author C.Marsh
#' @description TMB helper function this function returns a list of factors used in the map argument of the MakeADFun function
#' values with NA will not be estimated.
#' @param par_list a named list that you give to the par argument in the MakeADFun
#' @param pars_to_exclude a vector of strings with names of parameters you want to FIX in the objective object.
#' @param vec_elements_to_exclude a named list (names %in% pars_to_exclude) with number of elements = length(vec_pars_to_adjust). each list element
#' @param array_elements_to_exclude a named list (names %in% pars_to_exclude) with a matrix each row corresponds to an element with the first column being the array row index and second column being the array column index to fix
#' @param existing_map a named list that already contains NAs from previous fix_pars calls
#' contains a vector of elements that we want to exclude from estimation.
#' @return a list of factors used in the MakeADFun function
fix_pars <- function(par_list, pars_to_exclude, vec_elements_to_exclude = NULL, array_elements_to_exclude = NULL, existing_map = NULL) {
  if (!any(pars_to_exclude %in% names(par_list))) {
    stop(paste0("The parameters ", paste(pars_to_exclude[!pars_to_exclude %in% names(par_list)],collapse = " ")," in exclusion parameters could not be found in the 'par_list', please sort this out"))
  }
  pars = names(par_list)
  mapped_pars = list();
  existing_na_map = F
  if(!is.null(existing_map)) {
    mapped_pars = existing_map
    existing_na_map = T
  }
  
  if (!is.null(vec_elements_to_exclude)) {
    if (!all(names(vec_elements_to_exclude) %in% pars_to_exclude))
      stop("parameters names in vec_elements_to_exclude, need to also be in pars_to_exclude")
  }
  if (!is.null(array_elements_to_exclude)) {
    if (!all(names(array_elements_to_exclude) %in% pars_to_exclude))
      stop("parameters names in array_elements_to_exclude, need to also be in pars_to_exclude")
  }
  param_factor = 1;
  for(i in 1:length(pars)) {
    if(existing_na_map) {
      if(all(is.na(existing_map[[pars[i]]]))) {
        next;
      }
    }
    
    if (pars[i] %in% pars_to_exclude) {
      params_in_this_par = par_list[[pars[i]]];
      if (pars[i] %in% names(vec_elements_to_exclude)) {
        include_element_index = c(1:length(params_in_this_par))[-vec_elements_to_exclude[[pars[i]]]]
        params_vals = factor(rep(NA, length(params_in_this_par)), levels = factor(param_factor:(param_factor + length(include_element_index) - 1)))
        params_vals[include_element_index] = factor(param_factor:(param_factor + length(include_element_index) - 1))#, levels = factor(include_element_index))
        param_factor = param_factor + length(include_element_index)
        mapped_pars[[pars[i]]] = params_vals;
      } else if(pars[i] %in% names(array_elements_to_exclude)) {
        elements_to_drop = array_elements_to_exclude[[pars[i]]]
        mapped_vector = rep(NA, length(params_in_this_par))
        first_param_factor = param_factor
        vec_ndx = 1;
        ## TMB converts arrays to vectors down columns (not by rows)
        ## can handle up to 3-dimension arrays
        if(length(dim(params_in_this_par)) == 2) {
          for(col_ndx in 1:ncol(params_in_this_par)) {
            for(row_ndx in 1:nrow(params_in_this_par)) {
              dropping_this_element = F
              for(drop_ndx in 1:nrow(elements_to_drop)) {
                if(all(c(row_ndx, col_ndx) == elements_to_drop[drop_ndx,])) {
                  dropping_this_element = T
                  break;
                }
              }
              if(!dropping_this_element) {
                mapped_vector[vec_ndx] = param_factor
                param_factor = param_factor + 1
              }
              vec_ndx = vec_ndx + 1;
            }
          }
        } else if (length(dim(params_in_this_par)) == 3) {
          counter = 1;
          for(dim3_ndx in 1:dim(params_in_this_par)[3]) {
            for(dim2_ndx in 1:dim(params_in_this_par)[2]) {
              for(dim1_ndx in 1:dim(params_in_this_par)[1]) {
                ## check if we need to drop this value
                dropping_this_element = F
                for(drop_ndx in 1:nrow(elements_to_drop)) {
                  if(all(c(dim1_ndx, dim2_ndx, dim3_ndx) == elements_to_drop[drop_ndx,])) {
                    dropping_this_element = T
                    break;
                  }
                }
                if(!dropping_this_element) {
                  mapped_vector[vec_ndx] = param_factor
                  param_factor = param_factor + 1
                }
                vec_ndx = vec_ndx + 1;
              }
            }
          }
        } else {
          stop("this function can only deal with 2 or 3 dimensional arrays")
        }
        mapped_vector = factor(mapped_vector, levels = first_param_factor:max(mapped_vector, na.rm = T))
        mapped_pars[[pars[i]]] = mapped_vector;
      } else {
        ## exclude entire parameters
        mapped_pars[[pars[i]]] = rep(factor(NA),length(params_in_this_par));
        n_params_to_exclude = nrow(vec_elements_to_exclude[[pars[i]]])
      }
    } else {
      params_in_this_par = par_list[[pars[i]]];
      params_vals = factor(param_factor:(param_factor + length(params_in_this_par) - 1))
      param_factor = param_factor + length(params_in_this_par)
      mapped_pars[[pars[i]]] = params_vals
    }
  }
  return(mapped_pars);
}
