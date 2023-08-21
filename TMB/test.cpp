#include <TMB.hpp>
#include "AuxillaryFuns.hpp"

template<class Type>
  Type objective_function<Type>::operator() () {
  /*
    * Declare Namespace
  */
    using namespace density;
  // model dimensions
  DATA_VECTOR(ages);    // assumes min(ages) >= 1
  int n_ages = ages.size();
  DATA_ARRAY(propMat);                // Proportion Mature (used in SSB)    dim = n_ages x n_years
  DATA_ARRAY(mean_weight_at_age);        // Stock Mean weight used in SSB + survey calculation  dim = n_ages x n_years
  DATA_VECTOR(propZ_ssb);             // proportion of Z for SSB, length n_years
  DATA_SCALAR(natMor);                // Instantaneous natural mortality        
  DATA_SCALAR(steepness);             // Instantaneous natural mortality     
  DATA_SCALAR(mean_weight_a);         // a parameter mean weight from length 
  DATA_SCALAR(mean_weight_b);         // b parameter mean weight from length    
  DATA_SCALAR(R0);         // b parameter mean weight from length    
  DATA_SCALAR(f_a50 );            // logistic fishery ogive parameters for each fishery
  DATA_SCALAR(f_ato95);          // logistic fishery ogive parameters for each fishery
  DATA_SCALAR(F_init); 
  PARAMETER(dummy_variable);
  
  array<Type> N(n_ages, 2);
  // Convert mean length to mean weight
  int age_ndx  = 0;
  vector<Type> fishery_selectivity = logistic_ogive_alt(ages, f_a50, f_ato95);

  // Initialise Partition
  for(age_ndx = 0; age_ndx < n_ages; ++age_ndx) 
    N(age_ndx, 0) = R0 * exp(-(ages(age_ndx)) * natMor);
  // plus group
  N(n_ages - 1, 0) = R0 * exp(- ages(n_ages - 1) * natMor) / (1.0 - exp(-natMor));
  /*
   * Calcualte SPR for B0
   */
  vector<Type> Na(ages.size());
  Type F = 0;
  Na(0) = 1.0;
  vector<Type> Za = natMor + fishery_selectivity * F;
  vector<Type> Sa = exp(-Za);
  for(unsigned age_iter = 1; age_iter < ages.size(); ++age_iter) {
    Na(age_iter) = Na(age_iter - 1) * Sa(age_iter - 1);
  }
  // plus group
  Na(ages.size() - 1) = Na(n_ages - 1)  / (1.0 - exp(-Za(n_ages - 1)));
  Type SPR_0 = (Na * mean_weight_at_age.col(0).vec() *  propMat.col(0).vec()).sum();
  // redo with F-init
  F = F_init;
  vector<Type> Na_F(ages.size());
  Na_F(0) = 1.0;
  Za = natMor + fishery_selectivity * F;
  Sa = exp(-Za);
  for(unsigned age_iter = 1; age_iter < ages.size(); ++age_iter) {
    Na_F(age_iter) = Na_F(age_iter - 1) * Sa(age_iter - 1);
  }
  // plus group
  Na_F(ages.size() - 1) = Na_F(n_ages - 1)  / (1.0 - exp(-Za(n_ages - 1)));
  Type SPR_Finit = (Na_F * mean_weight_at_age.col(0).vec() *  propMat.col(0).vec()).sum();
  
  // Try using the function
  Type SPR_0_alt = get_SPR(Type(0.0), fishery_selectivity, natMor, mean_weight_at_age.col(0).vec(), propMat.col(0).vec(),Type(0.0),ages);
  Type SPR_Finit_alt = get_SPR(F_init, fishery_selectivity, natMor, mean_weight_at_age.col(0).vec(), propMat.col(0).vec(),Type(0.0),ages);
  Type SPR_0_alt2 = get_SPR_alt(Type(0.0), fishery_selectivity, natMor, mean_weight_at_age.col(0).vec(), propMat.col(0).vec(),propZ_ssb(0),ages);
  Type SPR_finit_alt2 = get_SPR_alt(F_init, fishery_selectivity, natMor, mean_weight_at_age.col(0).vec(), propMat.col(0).vec(),propZ_ssb(0),ages);
  REPORT(SPR_0);
  REPORT(SPR_0_alt);
  REPORT(SPR_0_alt2);
  REPORT(SPR_Finit);
  REPORT(SPR_Finit_alt);
  REPORT(SPR_finit_alt2);
  REPORT(Na);
  REPORT(Na_F);
  REPORT(Za);
  REPORT(propZ_ssb);
  REPORT(fishery_selectivity);
  return 0.0;
}
