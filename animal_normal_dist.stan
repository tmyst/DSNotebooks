// 2-5-1
data{
  int N;
  vector[N] animal_num;
}

parameters{
  real<lower = 0> mu;
  real<lower = 0> sigma;
}

model{
  animal_num ~ normal(mu, sigma);
}

generated quantities{
  vector[N] pred;
  for( i in 1:N){
    pred[i] = normal_rng(mu, sigma);
  }
}
