data {
  int N;
  int X[N];
}

parameters {
  real<lower = 0, upper = 1> q;
}

model{
  for(n in 1:N){
    X[n] ~ bernoulli(q);
  }
}

generated quantities{
  real log_lik[N];
  for(n in 1:N){
    log_lik[n] = bernoulli_lpmf(X[n]|q);
  }
}
