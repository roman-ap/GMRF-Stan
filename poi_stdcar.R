##### Poisson regression with CAR effects 
##### using Sweden's counties as irregular lattice 

### Set-up
library(sf)
library(dplyr)
library(Matrix)
library(cmdstanr)
library(extraDistr)
library(withr)
library(ggplot2)
library(latex2exp)

### Sweden's counties
# Spatial geometries
counties <- read_sf("./Sweden/counties.shp")
counties <- counties[-8,] #Removing Gotland
plot(counties)

### Spatial structure
## CAR(1) structure
# Number of counties
S<-nrow(counties)
# Binary spatial relationships  
S_match<-st_relate(counties,counties,pattern="****11212",sparse=TRUE)
# Sparse matrices
Ds_ii<-sapply(S_match,length)
D<-sparseMatrix(i=c( 1:S ),
                j=c( 1:S ),
                x=c( Ds_ii ),
                dims=c(S,S))
A<-sparseMatrix(i=c( rep(1:S, Ds_ii) ),
                j=c( unlist(S_match) ),
                x=c( rep(1, length(unlist(S_match))) ),
                dims=c(S,S))
C<-sparseMatrix(i=c( rep(1:S, Ds_ii) ),
                j=c( unlist(S_match) ),
                x=c( rep(1/Ds_ii, Ds_ii) ),
                dims=c(S,S))
# Eigenvalues of C_space
lambdas<-Schur(C, vectors = FALSE)$EValues

### Realisation of MARCAR effects
## Hyperpriors for the mean, standard deviation and spatial dependence
pm_mu = 1
psd_mu = 0.1
pm_sigma = 0
psd_sigma = 0.1
pm_rho = 0.5
psd_rho = 0.1
## Realisation
with_seed(87L,{
  mu_obs <- rnorm(n = 1, pm_mu, psd_mu)
  sigma_obs <- rtnorm(n = 1, pm_sigma, psd_sigma, a=0)
  rho_obs <- rtnorm(n = 1, pm_rho, psd_rho, a=1/min(lambdas), b=1/max(lambdas))
  Q = (D-rho_obs*A)
  M = 1
  Q = forceSymmetric(Q)
  L = Cholesky(Q,LDL=F)
  z = Matrix(rnorm(S*M),S,M)
  u_per = solve(L,z,system="Lt")
  u = solve(L,u_per,system="Pt")
  eta = mu_obs + sigma_obs*u
})
## Realization of a Poisson regression with rates exp(eta)
## where eta follows a CAR distribution   
y_obs <- numeric(S)
mu_y <- numeric(S)
for(s in 1:S){
  mu_y[s] = exp(eta[s])
  y_obs[s] = rpois(1, mu_y[s])
}

### Bayesian estimation via HMC (Stan)
# Stan model using a sparse representation for the MARCAR distribution
car_model <- cmdstan_model("./poi_stdcar.stan")
# Data inputs for the Stan models
A_w<-A@x
A_v<-A@i +1
A_u<-A@p +1
log_det_D<-sum(log(Ds_ii))
# Stan inputs
data_list <- list(S = S,
                  y = y_obs,
                  pm_mu = pm_mu,
                  psd_mu = psd_mu,
                  pm_sigma = pm_sigma,
                  psd_sigma = psd_sigma,
                  pm_rho = pm_rho,
                  psd_rho = psd_rho,
                  Ds_ii = Ds_ii,
                  log_det_D = log_det_D,
                  lambdas = lambdas,
                  As_sparse = length(A_w),
                  A_w = A_w,
                  A_v = A_v,
                  A_u = A_u)
## MCMC sampling
# Sampler
car_fit <- car_model$sample(data = data_list,
                            chains = 4,
                            parallel_chains = 4,
                            refresh = 100,
                            iter_warmup = 1000,
                            iter_sampling = 1000
)
# Posterior summaries
car_fit$summary(variables = c("mu", "sigma", "rho"))



### Visualizations (Prior and posterior distributions + observed value)
## Data frame with posterior draws 
car_draws <- car_fit$draws(format = "df")
# Mean
ggplot(data.frame(draw=car_draws$`mu`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu, sd=psd_mu), colour="black", linewidth=1) +
  xlim(pm_mu - 6*psd_mu, pm_mu + 6*psd_mu) +
  geom_vline(xintercept = mu_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu$"),
       y = TeX("Density"), x = TeX("$\\mu$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation
ggplot(data.frame(draw=car_draws$`sigma`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dhnorm, args = list(sigma=psd_sigma), colour="black", linewidth=1) +
  xlim(0, pm_sigma + 6*psd_sigma) +
  geom_vline(xintercept = sigma_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma$"),
       y = TeX("Density"), x = TeX("$\\sigma$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Spatial dependence
ggplot(data.frame(draw=car_draws$`rho`)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_rho, sd=psd_rho), colour="black", linewidth=1) +
  xlim(pm_rho - 6*psd_rho, pm_rho + 6*psd_rho) +
  geom_vline(xintercept = rho_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\rho$"),
       y = TeX("Density"), x = TeX("$\\rho$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
