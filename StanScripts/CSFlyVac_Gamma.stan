//mean and var parameter estimation for FlyVac using Bayesian inference
//inter trial time
//created by: jamilla a
//created on: 2018-10-15
//updated: 2018-10-22

functions {
// user defined fxns
}

data {
// required data for the model

int <lower=0> N;	//size of data
int <lower=0> L;	//num of treatments
real y[N];	//modeled data
matrix[N,L] X;	//matrix of treatments

}

transformed data {
// data transforms and constants - evaluated once

}

parameters {
// parameters for estimation

vector[L] beta;	//regression parameters
vector[L] gamma;//regression parameters
	
}

transformed parameters {
// parameter transform - evaluated once per leapfrog step

vector[N] m; //linear predictor (mean)
vector[N] v; //linear predictor (var)
vector[N] a; //shape parameter for the gamma distribution
vector[N] b; //rate parameter for the gamma distribution


m = exp(X*beta); //log link
v = exp(X*gamma);

for (i in 1:N) { 
    a[i] = (m[i]*m[i])/v[i];   
  }
for (i in 1:N) { 
    b[i] = m[i]/v[i];   
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

y ~ gamma(a,b);

}

generated quantities {
// derive quantities based on parameters and data once per sample (post-leapfrog steps)

vector[N] y_rep;

for (i in 1:N) { 
   y_rep[i] = gamma_rng(a[i], b[i]);
}

}