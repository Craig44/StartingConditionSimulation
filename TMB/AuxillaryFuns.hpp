/*
 * isNA
 */
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}
/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

// Zerofun functions see - for a discussion on this https://github.com/kaskr/adcomp/issues/7
template<class Type>
Type posfun(Type x, Type eps, Type &pen) {
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  Type xp = -(x/eps-1);
  return CppAD::CondExpGe(x,eps,x,
                          eps*(1/(1+xp+pow(xp,2)+pow(xp,3)+pow(xp,4)+pow(xp,5))));
}

/*
 template<class Type>
 Type posfun(Type x, Type eps, Type &pen){
 pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
 return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
 }
 */

/*
 * 
 */
template <class Type> 
Type square(Type x){return x*x;}

template <class Type> 
vector<Type> square(vector<Type>& x) {
  return x*x;
}
// transform Y -Inf-Inf -> X bound lb - ub
template <class Type> 
Type invlogit_general(Type& Y, Type& lb, Type& ub) {
  return(lb + (ub - lb) * (1 / (1 + exp(-Y))));
}
// transform Y -Inf-Inf -> X bound lb - ub
template <class Type> 
vector<Type> invlogit_general(vector<Type>& Y, Type& lb, Type& ub) {
  vector<Type> X(Y.size());
  for(int i = 0; i < X.size(); ++i) {
    X(i) = lb + (ub - lb) * (1 / (1 + exp(-Y(i))));
  }
  return(X);
}


// logistic ogive function
template <class Type> 
vector<Type> logistic_ogive(vector<Type> ages, Type sel_50, Type sel_95) {
  std::cout << "logistic_ogive\n";
  int n_ages = ages.size();
  vector<Type> logis(n_ages);
  for (int age = 0;  age < n_ages; ++age) {
    logis[age] = Type(1.0) / (Type(1.0) + pow(Type(19.0), (sel_50 - ages[age]) / sel_95));
  }
  return logis;
}
/*
 * Geometric mean
 */
template <class Type> 
Type geo_mean(vector<Type>& x){
  return exp((log(x).sum())/x.size());
}

/*
 * centred log transform
 */
template <class Type> 
vector<Type> crl(vector<Type>& x) { 
  return log(x / geo_mean(x));
}

// https://github.com/timjmiller/wham/blob/master/src/age_comp_sim.hpp
template<class Type>
vector<Type> rmultinom_alt(vector<Type> prob, Type N)
{
  //multinomial
  int dim = prob.size();
  vector<Type> x(dim);
  int Nint = CppAD::Integer(N);
  x.setZero();
  for(int i = 0; i < Nint; i++)
  {
    Type y = runif(0.0,1.0);
    for(int a = 0; a < dim; a++) if(y < prob.head(a+1).sum())
    {
      x(a) += 1.0;
      break;
    }
  }
  return x;
}

/*
 * rmultinomm - for simulate call
 */
template <class Type> 
vector<Type> rmultinom(vector<Type> prob, Type N) { 
  vector<Type> sim_X(prob.size());
  sim_X.setZero();
  // Now simulate using the uniform random variable
  Type rng_uniform;
  Type cumulative_expect;
  while(N > 0) {
    rng_uniform = runif(Type(0),Type(1));
    //std::cout << rng_uniform << " ";
    cumulative_expect = 0.0;
    for (unsigned i = 0; i < prob.size(); ++i) {
      cumulative_expect += prob[i];
      if (cumulative_expect >= rng_uniform) {
        sim_X[i] += 1.0;
        break;
      }
      //sim_X[prob.size() - 1] += 1.0;
    }
    N -= 1;
  }
  //std::cout << "\n";
  return(sim_X);
}

/*
 *  Simulate a single draw from a multinomial-dirichlet distribution
 */
template <class Type> 
vector<Type> rdirichletmulti(vector<Type> fitted_props, Type& n_eff, Type& theta) {
  vector<Type> dirichlet_draw(fitted_props.size()); 
  for(int ndx = 0; ndx < fitted_props.size(); ndx++) 
    dirichlet_draw(ndx) = rgamma(fitted_props(ndx) * theta * n_eff, (Type)1.0);// shape, rate = 1.0
  
  Type dirich_total = dirichlet_draw.sum();
  dirichlet_draw /= dirich_total;
  return(rmultinom(dirichlet_draw, n_eff));
}
/*
 * inverse centred log transform
 * Up to a constant so standardise
 */
template <class Type> 
vector<Type> inv_crl(vector<Type>& y){ 
  return exp(y) / (exp(y)).sum();
}

// Beverton-Holt SR relationship function
template<class Type>
Type BevertonHolt(Type SSB, Type B0, Type h) {
  Type ssb_ratio = SSB / B0;
  Type part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio));
  return (ssb_ratio / part_2);
}

// Beverton-Holt SR relationship function without equilibrium assumptions
template<class Type>
Type BevertonHoltNoEquil(Type a, Type b, Type SSB) {
  Type Rt = (a + SSB) / (SSB * b);
  return Rt;
}


// Calculate spawner per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param paa vector of proportions at age
 * @param propZ_ssb scalar interpolating SSB
 * @param ages vector of ages
 * @return SPR
 */

template<class Type>
Type get_SPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& paa, Type propZ_ssb, vector<Type>& ages) {
  Type Na = 1;
  vector<Type> Za = natural_mortality + sel * F;
  vector<Type> Sa = exp(-Za);
  Type SPR = Na * exp(-Za(0) * propZ_ssb) * waa(0) * paa(0);
  for(unsigned age_iter = 1; age_iter < (ages.size() * 4); ++age_iter) {
    if(age_iter >= ages.size() ) {
      SPR += Na * exp(-Za(ages.size() - 1) * propZ_ssb) * waa(ages.size() - 1) * paa(ages.size() - 1);
      Na *= Sa(ages.size() - 1);
    } else {
      SPR += Na * exp(-Za(age_iter) * propZ_ssb) * waa(age_iter) * paa(age_iter); 
      Na *= Sa(age_iter);
    }
  }
  return SPR;
}

// Calculate Yeild per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return SPR
 */
template<class Type>
Type get_YPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages) {
  Type Na = 1;
  Type YPR = (sel(0) * F) / (natural_mortality + sel(0) * F) * (Type(1) - exp(-(natural_mortality + sel(0) * F))) * Na * waa(0);
  for(unsigned age_iter = 1; age_iter < ages.size() * 4; ++age_iter) {
    if(age_iter >= ages.size() ) {
      Na *= exp(-(natural_mortality + sel(ages.size() - 1) * F));
      YPR += (sel(ages.size() - 1) * F) / (natural_mortality + sel(ages.size() - 1) * F) * (Type(1) - exp(-(natural_mortality + sel(ages.size() - 1) * F))) * Na * waa(ages.size() - 1);
    } else {
      Na *= exp(-(natural_mortality + sel(age_iter) * F));
      YPR += (sel(age_iter) * F) / (natural_mortality + sel(age_iter) * F) * (Type(1) - exp(-(natural_mortality + sel(age_iter) * F))) * Na * waa(age_iter);
      
    }
  }
  return YPR;
}
// Calculate derivative for Yeild per recruit
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return derivative YPR
 */
template<class Type>
Type get_dYPR(Type F, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages) {
  Type ln_F = log(F);
  Type h = 0.001;
  Type v = -get_YPR(exp(ln_F + 2.0 * h), sel, natural_mortality, waa, ages) + 8.0 * get_YPR(exp(ln_F + h), sel, natural_mortality, waa, ages) - 8.0 * get_YPR(exp(ln_F - h), sel, natural_mortality, waa, ages) + get_YPR(exp(ln_F - 2.0 * h), sel, natural_mortality, waa, ages);
  Type g = v / (12.0 * h);
  return(g / F);
}
// Calculate Equilbrium SSB for a given F
/*
 * @param F F to test
 * @param sel vector of ages selectivity
 * @param natural_mortality natural mortality
 * @param waa vector of weights at age
 * @param ages vector of ages
 * @return SPR
 */
template<class Type>
Type get_SSBe(Type& F, vector<Type>& N_init, vector<Type>& sel, Type& natural_mortality, vector<Type>& waa, vector<Type>& ages, vector<Type>& paa, Type propZ_ssb, int n_runs, int maxAgePlusGroup) {
  
  vector<Type> Fa = F * sel;
  vector<Type> Za = Fa + natural_mortality;
  vector<Type> Sa = exp(-Za);
  array<Type> N_equilibrium(ages.size(), n_runs);
  N_equilibrium.col(0) = N_init;
  
  // Run annual cycle 
  for(int year_ndx = 1; year_ndx < n_runs; ++year_ndx) {
    // Initial recruitment
    N_equilibrium(0, year_ndx) = N_init(0);
    // Ageing + Z
    for(int age_ndx = 1; age_ndx < ages.size() ; ++age_ndx) 
      N_equilibrium(age_ndx, year_ndx) = N_equilibrium(age_ndx - 1, year_ndx - 1) * Sa(age_ndx - 1);
    if(maxAgePlusGroup == 1) {
      N_equilibrium(ages.size() - 1, year_ndx) = N_equilibrium(ages.size()  - 2, year_ndx - 1) * Sa(ages.size()  - 2) +
        N_equilibrium(ages.size()  - 1, year_ndx - 1) * Sa(ages.size()  - 1);
    }  
  }
  
  // Calculate SSBs an interpolation bewtween the year, starting with previous years Paritition
  Type ssbe = 0;
  for(int age_ndx = 0; age_ndx < ages.size(); ++age_ndx) 
    ssbe += N_equilibrium(age_ndx, n_runs - 1) * exp(-Za(age_ndx) * propZ_ssb) * paa(age_ndx) * waa(age_ndx);
  
  return(ssbe);
}
