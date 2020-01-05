//beer sales ab stan
data{
  int N;
  vector[N] sales_a;
  vector[N] sales_b;
}

parameters{
  real mu_a;
  real<lower=0> sigma_a;
  real mu_b;
  real<lower=0> sigma_b;
}

model{
  sales_a ~ normal(mu_a, sigma_a);
  sales_b ~ normal(mu_b, sigma_b);
}

generated quantities{
  real diff;
  diff = mu_b - mu_a;
}
