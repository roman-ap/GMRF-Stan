##### Poisson regression with MARCAR effects 
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

### Spatial and temporal structure
## CAR(1) structure
# Number of counties
S<-nrow(counties)
# Binary spatial relationships  
S_match<-st_relate(counties,counties,pattern="****11212",sparse=TRUE)
# Sparse matrices
Ds_ii<-sapply(S_match,length)
D_space<-sparseMatrix(i=c( 1:S ),
                      j=c( 1:S ),
                      x=c( Ds_ii ),
                      dims=c(S,S))
A_space<-sparseMatrix(i=c( rep(1:S,Ds_ii) ),
                      j=c( unlist(S_match) ),
                      x=c( rep(1,length(unlist(S_match))) ),
                      dims=c(S,S))
C_space<-sparseMatrix(i=c( rep(1:S,Ds_ii) ),
                      j=c( unlist(S_match) ),
                      x=c( rep(1/Ds_ii,Ds_ii) ),
                      dims=c(S,S))
# Eigenvalues of C_space
lambdas<-Schur(C_space, vectors = FALSE)$EValues
## AR(1) structure
# Number of years
T<-22L
# Sparse matrices
I_time<-Diagonal(T, 1)
II_time<-sparseMatrix(i=c( 2:(T-1) ),
                      j=c( 2:(T-1) ),
                      x=c( rep(1, T-2) ),
                      dims=c(T,T))
A_time<-bandSparse(T, 
                   k = -1,
                   diagonals = list(rep(1, T-1)),
                   symmetric = TRUE)

### Realisation of MARCAR effects
## Hyperpriors for the mean, standard deviation, spatial and temporal dependence
pm_mu = -1.3
psd_mu = 0.1
pm_sigma = 0
psd_sigma = 0.1
pm_rho_time = 0.5
psd_rho_time = 0.1
pm_rho_space = 0.5
psd_rho_space = 0.1
## Simulations
with_seed(27L,{
  # Draw of the parameters using the hyperpriors
  mu_obs <- rnorm(n = 1, pm_mu, psd_mu)
  sigma_obs <- rtnorm(n = 1, pm_sigma, psd_sigma, a = 0)
  rho_time_obs <- rtnorm(n = 1, pm_rho_time, psd_rho_time, a = -1, b = 1)
  rho_space_obs <- rtnorm(n = 1, pm_rho_space, psd_rho_space, a = 1/min(lambdas), b = 1/max(lambdas))
  # Draw of the MARCAR effects 
  Q_time = (I_time + rho_time_obs^{2}*II_time - rho_time_obs*A_time)
  Q_space = (D_space - rho_space_obs*A_space)
  Q = Matrix::kronecker(Q_time, Q_space)
  Q = forceSymmetric(Q)
  M = 1
  L = Cholesky(Q,LDL=F)
  z = Matrix(rnorm(T*S*M),T*S, M)
  u_per = solve(L,z,system="Lt")
  u = solve(L,u_per,system="Pt")
  eta = mu_obs + sigma_obs*u
})
## Realization of a Poisson regression with rates exp(eta)
## where eta follows a MARCAR distribution   
mu_y = numeric(T*S)
y_obs = numeric(T*S)
for(i in 1:(T*S)){
  mu_y[i] = exp(eta[i])
  y_obs[i] = rpois(1,mu_y[i])
}

### Bayesian estimation via HMC (Stan)
# Stan model using a sparse representation for the MARCAR distribution
marcar_model <- cmdstan_model("./poi_stdmarcar.stan")
# Kronecker products
ItDs<-Matrix::kronecker(I_time, D_space)
ItAs<-Matrix::kronecker(I_time, A_space)
AtDs<-Matrix::kronecker(A_time, D_space)
AtAs<-Matrix::kronecker(A_time, A_space)
IItDs<-Matrix::kronecker(II_time, D_space)
IItAs<-Matrix::kronecker(II_time, A_space)
# Data inputs for the Stan models
log_det_Ds<-sum(log(Ds_ii))
ItDs_ii = ItDs@x
N_ItAs = length(ItAs@x)
ItAs_w = ItAs@x
ItAs_v = ItAs@i + 1
ItAs_u = ItAs@p + 1
N_AtDs = length(AtDs@x)
AtDs_w = AtDs@x
AtDs_v = AtDs@i + 1
AtDs_u = AtDs@p + 1
N_AtAs = length(AtAs@x)
AtAs_w = AtAs@x
AtAs_v = AtAs@i + 1
AtAs_u = AtAs@p + 1
IItDs_ii = diag(IItDs)
N_IItAs = length(IItAs@x)
IItAs_w = IItAs@x
IItAs_v = IItAs@i + 1
IItAs_u = IItAs@p + 1
# Stan inputs
marcar_data <- list(T = T,
                    S = S,
                    y = y_obs,
                    pm_mu = pm_mu,
                    psd_mu = psd_mu,
                    pm_sigma = pm_sigma,
                    psd_sigma = psd_sigma,
                    pm_rho_time = pm_rho_time,
                    psd_rho_time = psd_rho_time,
                    pm_rho_space = pm_rho_space,
                    psd_rho_space = psd_rho_space,
                    log_det_Ds = log_det_Ds,
                    lambdas = lambdas,
                    ItDs_ii = ItDs_ii,
                    N_ItAs = N_ItAs,
                    ItAs_w = ItAs_w,
                    ItAs_v = ItAs_v,
                    ItAs_u = ItAs_u,
                    N_AtDs = N_AtDs,
                    AtDs_w = AtDs_w,
                    AtDs_v = AtDs_v,
                    AtDs_u = AtDs_u,
                    N_AtAs = N_AtAs,
                    AtAs_w = AtAs_w,
                    AtAs_v = AtAs_v,
                    AtAs_u = AtAs_u,
                    IItDs_ii = IItDs_ii,
                    N_IItAs = N_IItAs,
                    IItAs_w = IItAs_w,
                    IItAs_v = IItAs_v,
                    IItAs_u = IItAs_u)
## MCMC sampling
# Sampler
marcar_fit <- marcar_model$sample(data = marcar_data,
                                  chains = 4,
                                  parallel_chains = 4,
                                  refresh = 100,
                                  iter_warmup = 1000,
                                  iter_sampling = 1000
                                  )
# Posterior summaries
marcar_fit$summary(variables = c("mu", "sigma", "rho_time", "rho_space"))

### Visualizations (Prior and posterior distributions + observed value)
## Data frame with posterior draws 
marcar_draws <- marcar_fit$draws(format = "df")
# Mean
ggplot(data.frame(draw=marcar_draws$mu)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_mu, sd=psd_mu), colour="black", linewidth=1) +
  xlim(pm_mu - 6*psd_mu, pm_mu + 6*psd_mu) +
  geom_vline(xintercept = mu_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\mu$"),
       y = TeX("Density"), x = TeX("$\\mu$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Standard deviation
ggplot(data.frame(draw=marcar_draws$sigma)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dhnorm, args = list(sigma=psd_sigma), colour="black", linewidth=1) +
  xlim(max(0, pm_sigma - 6*psd_sigma), pm_sigma + 6*psd_sigma) +
  geom_vline(xintercept = sigma_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\sigma$"),
       y = TeX("Density"), x = TeX("$\\sigma$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Temporal dependence
ggplot(data.frame(draw=marcar_draws$rho_time)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_rho_time, sd=psd_rho_time), colour="black", linewidth=1) +
  xlim(max(-1, pm_rho_time - 6*psd_rho_time), max(1, pm_rho_time + 6*psd_rho_time)) +
  geom_vline(xintercept = rho_time_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\rho_{T}$"),
       y = TeX("Density"), x = TeX("$\\rho_{T}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
# Spatial dependence
ggplot(data.frame(draw=marcar_draws$rho_space)) +
  geom_density(aes(x=draw, y=after_stat(density)), color="pink", fill="pink") +
  stat_function(fun = dnorm, args = list(mean=pm_rho_space, sd=psd_rho_space), colour="black", linewidth=1) +
  xlim(max(1/min(lambdas), pm_rho_space - 6*psd_rho_space), min(1/max(lambdas), pm_rho_space + 6*psd_rho_space)) +
  geom_vline(xintercept = rho_space_obs, color="blue") +
  theme_minimal() +
  labs(title = TeX("Prior and posterior distribution of $\\rho_{S}$"),
       y = TeX("Density"), x = TeX("$\\rho_{S}$")) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
