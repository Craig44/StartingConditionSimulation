#' compare Bai Li's non-equilibrium initialisation method
#' vs ours
#' The crux of the difference is that when there is an F-init, they scale the exp(-Z) age-structure by
#' R_eq not R_0 (which we do), where R_eq is based on the an adjusted R0 from the adjusted spawner per-recruit function with the assumed F-init.
#' 

source("AuxillaryFunctions.R")
library(dplyr)
library(ggplot2)
library(reshape2)
library(TMB)
library(tidyr)
fig_dir = file.path(DIR$fig, "AmalgamBiology")
if(!dir.exists(fig_dir))
  dir.create(fig_dir)

this_bio = readRDS(file = file.path(DIR$data, "Amalgam_biology.RDS"))

biol_df = data.frame(age = this_bio$ages, length_at_age = vonbert(this_bio$ages, this_bio$K, L_inf = this_bio$L_inf, t0 = this_bio$t0))
biol_df$maturity = logis_sel(this_bio$ages, this_bio$m_a50, this_bio$m_ato95)
biol_df$weight_at_age  = this_bio$a * biol_df$length_at_age^this_bio$b
biol_df$fishery_selectivity = logistic_sel_alternative(this_bio$ages, this_bio$f_a50, this_bio$f_ato95)
nages = length(biol_df$maturity)

## EQ E2.1
N.pr0=rep(1,nages) #Number of spawners per recruit at age
for (a in 1:(nages-1)) {
  N.pr0[a+1]=N.pr0[a]*exp(-this_bio$M)
}
N.pr0[nages]=N.pr0[nages]/(1-exp(-this_bio$M))  #Plus group
Phi.0=sum(N.pr0*maturity*mean_weight)     #Spawners per recruit based on mature female biomass


#' based on equation E3.2 from Bai Li et. al 2021 paper
#' returns the spawner per recruit under a given F scenario
phi_F <- function(F, selectivity, M, mean_weight, maturity) {
  phi_F = vector()
  phi_F[1] = 1
  nages = length(mean_weight)
  Za = F * selectivity + M
  for(a in 1:(length(mean_weight) - 1))
    phi_F[a + 1] = phi_F[a] * exp(-Za[a]) 
  ## plus group
  phi_F[nages] = phi_F[nages] / (1 - exp(-Za[nages]))
  
  ## get spawner per recruit
  SPR = sum(phi_F * mean_weight * maturity)
  return(SPR)
}

## test 
phi_F(F = 0, selectivity = biol_df$fishery_selectivity, M = this_bio$M, mean_weight = biol_df$weight_at_age, maturity = biol_df$maturity)
Phi_0 = Phi.0
phi_0.2 = phi_F(F = 0.2, selectivity = biol_df$fishery_selectivity, M = this_bio$M, mean_weight = biol_df$weight_at_age, maturity = biol_df$maturity)

R_eq <- function(R0, h, phi_0, phi_F) {
  return (R0*(4*h*phi_F-(1-h)*phi_0)/((5*h-1)*phi_F))
}



#Initial conditions assumes equilibrium age structure given initial F
N.pr1=rep(1,nages) #Number of spawners per recruit at age
Z=0.2*biol_df$fishery_selectivity+this_bio$M
for (a in 1:(nages-1))
{N.pr1[a+1]=N.pr1[a]*exp(-Z[a])}
N.pr1[nages]=N.pr1[nages]/(1-exp(-Z[nages])) #Plus group
Phi.F=sum(N.pr1*maturity*mean_weight) #Spawners per recruit based on mature female biomass
## intial age-structure
R_eq(this_bio$R0, this_bio$h, Phi_0, phi_0.2)*N.pr1
