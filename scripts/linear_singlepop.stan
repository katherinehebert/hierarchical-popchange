data {
  int n; //number of observations in the data
  vector[n] biomass; //vector of length n for the population's biomass
  vector[n] year; //vector of length n for year
  // posterior predictions
  int<lower=0> npost;
}

parameters {
  real <lower=0> alpha; //the intercept parameter
  real beta_year; //slope parameter for year
  real<lower=0> sigma; //model variance parameter
}

model {
  //linear predictor mu
  vector[n] mu;
  
  //write the linear equation
  mu = alpha + beta_year * year;
  
  //likelihood function
  biomass ~ normal(mu, sigma);
  
  // very uncertain priors
  alpha ~ normal(0, .5);
  beta_year ~ normal(0, .1);
  sigma ~ exponential(1);
  
}

generated quantities {
  vector[npost] y_rep;
  
  for(i in 1:npost){
    //replications for the posterior predictive distribution
    y_rep[i] = normal_rng(alpha + beta_year * year[i], sigma);
  }
}
