//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower = 0> N;
  int<lower = 0> edu[N];
  int<lower = 0, upper = 1> sex[N];
  int<lower = 0, upper = 1> y[N];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real b[3];
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for(i in 1:N){
    y[i] ~ bernoulli_logit(b[1] + b[2]*edu[i] + b[3]*sex[i]);
  }
}

