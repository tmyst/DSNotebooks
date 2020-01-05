//beer_temperature by design matrix
data{
  int N;
  int K;
  vector[N] Y;
  matrix[N, K] X;
}

parameters{
  vector[K] b;
  real<lower = 0> sigma;
}

model{
  vector[N] mu = X * b;
  Y ~ normal(mu, sigma);
}
