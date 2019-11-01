//mean and var parameter estimation using Bayesian inference
//normal model - turn bias, switchiness, clumpiness, light-choice
//created by: jamilla a
//created on: 2018-08-29
//updated: 2019-01-10

functions {
// user defined fxns
}

data {
// required data for the model

int <lower=0> N;	//size of data
int <lower=0> L;	//num of treatments
vector[N] y;	//modeled data
matrix[N,L] X;	//matrix of treatments

}

transformed data {
// data transforms and constants - evaluated once
}

parameters {
// parameters for estimation

real a;	//mean intercept
vector[L] be;	//mean coeff for treatment
real <lower=0> v0;	//var intercept
vector[L] ve;	//var coeff for treatment

}

transformed parameters {
// parameter transform - evaluated once per leapfrog step

vector<lower=0>[N] sigma;
sigma = sqrt(v0 + X * ve);

}

model {
// Priors
a ~ cauchy(0, 2.5); //broad prior on intercept and coefs
be ~ cauchy(0, 2.5);
v0 ~ cauchy(0, 2.5);
ve ~ cauchy(0, 2.5);

// define log probability fxn

y ~ normal(a + X * be, sigma);

}

generated quantities {
// derive quantities based on parameters and data once per sample (post-leapfrog steps)

vector[N] y_rep;

for (i in 1:N) { 
   y_rep[i] = normal_rng(a + X[i]*be, sigma[i]);
}

}