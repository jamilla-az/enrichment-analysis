//mean and var parameter estimation for Ymaze using Bayesian inference
//num turns
//created by: jamilla a
//created on: 2018-10-08
//updated: 2018-1-16

functions {
// user defined fxns
}

data {
// required data for the model

int <lower=0> N;	//size of data
int <lower=0> L;	//num of treatments
int y[N];	//modeled data
matrix[N,L] X;	//matrix of treatments

}

transformed data {
// data transforms and constants - evaluated once

}

parameters {
// parameters for estimation

vector[L] beta;	//regression parameters
vector[L] gamma;

}

transformed parameters {
// parameter transform - evaluated once per leapfrog step

vector<lower=0>[N] m; //linear predictor (mean)
vector[N] v; //transformed variance
vector<lower=0>[N] phi;	//overdispersion parameter

m = exp(X*beta); //log link
phi = exp(X*gamma);

for (i in 1:N) { 
    v[i] = m[i]+((m[i]*m[i])/phi[i]);   
  }
}

model {
// Priors

beta[1] ~ cauchy(0, 10);   //set broad priors on intercept coefficients
gamma[1] ~ cauchy(0,10);

for (i in 2:L){
    beta[i] ~ cauchy(0,2.5); //set priors on regression coefficients
    gamma[i] ~ cauchy(0,2.5);
}


// define log probability fxn

y ~ neg_binomial_2(m,phi);

}

generated quantities {
// derive quantities based on parameters and data once per sample (post-leapfrog steps)

vector[N] y_rep;

for (i in 1:N) { 
   y_rep[i] = neg_binomial_2_rng(m[i], phi[i]);
}

}
