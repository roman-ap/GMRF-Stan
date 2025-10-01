functions {
#include ST_lpdf.stan
}

data {
  int<lower=1> T;
  int<lower=1> S;
  array[T*S] int<lower=0> y;
  //vector[T*S] log_offset;
  vector[S] lambdas;
  real pm_mu;
  real<lower=0> psd_mu;
  real<lower=0> pm_sigma;
  real<lower=0> psd_sigma;
  real<lower=-1, upper=1> pm_rho_time;
  real<lower=0> psd_rho_time;
  real<lower=1/min(lambdas), upper=1/max(lambdas)> pm_rho_space;
  real<lower=0> psd_rho_space;
  real log_det_Ds;
  vector[T*S] ItDs_ii;
  int<lower=1> N_ItAs;
  vector[N_ItAs] ItAs_w;
  array[N_ItAs] int ItAs_v;
  array[T*S+1] int ItAs_u;
  int<lower=1> N_AtDs;
  vector[N_AtDs] AtDs_w;
  array[N_AtDs] int AtDs_v;
  array[T*S+1] int AtDs_u;
  int<lower=1> N_AtAs;
  vector[N_AtAs] AtAs_w;
  array[N_AtAs] int AtAs_v;
  array[T*S+1] int AtAs_u;
  vector[T*S] IItDs_ii;
  int<lower=1> N_IItAs;
  vector[N_IItAs] IItAs_w;
  array[N_IItAs] int IItAs_v;
  array[T*S+1] int IItAs_u;
}

transformed data {
  vector[T*S] zeros;
  zeros = rep_vector(0, T*S);
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=-1, upper=1> rho_time;
  real<lower=1/min(lambdas), upper=1/max(lambdas)> rho_space;
  vector[T*S] phi;
}

transformed parameters {
    vector[T*S] eta;
    eta = mu + sigma*phi;
}

model {
  target += poisson_log_lpmf(y | eta);
  target += std_marcar_lpdf(phi | zeros,
                                  rho_time, rho_space,
                                  log_det_Ds, lambdas,
                                  ItDs_ii,
                                  ItAs_w, ItAs_v, ItAs_u,
                                  AtDs_w, AtDs_v, AtDs_u,
                                  AtAs_w, AtAs_v, AtAs_u,
                                  IItDs_ii, 
                                  IItAs_w, IItAs_v, IItAs_u,
                                  T, S);
  //target += normal_lpdf(rho_time | pm_rho_time, psd_rho_time);
  target += uniform_lpdf(rho_time | -1, 1);
  //target += normal_lpdf(rho_space | pm_rho_space, psd_rho_space);
  target += uniform_lpdf(rho_space | 1/min(lambdas), 1/max(lambdas));
  target += normal_lpdf(sigma | pm_sigma, psd_sigma);
  target += normal_lpdf(mu | pm_mu, psd_mu);
}
