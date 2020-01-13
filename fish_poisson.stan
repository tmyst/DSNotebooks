data{
  int N;
  int y[N];
  vector[N] temperature;
  vector[N] weathersunny;
  // 2.to use design matrix
  // int K;
  // matrix[N, K] X;
}

parameters{
  real Intercept;
  real b_temperature;
  real b_weathersunny;
  // 2.to use design matrix
  // vector[K] b;
}

model{
  vector[N] lambda = exp(Intercept + b_temperature * temperature + b_weathersunny * weathersunny);
  // 1.to make calculation efficient
  // vector[N] lambda_log = Intercept + b_temperature * temperature + b_weathersunny * weathersunny;
  // 2.to use design matrix
  // vector[N] lambda = X * b;
  
  y ~ poisson(lambda);
  // 1 or 2
  // y ~ poisson_log(lambda)
}

