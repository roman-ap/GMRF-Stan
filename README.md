The primary goal of this repository is to document the implementation of Gaussian Markov Randon Fields for Spatial, Temporal and Spatio-Temporal effects in Stan. `ST_lpdf.stan` contains the necessary log probability density functions, whereas `ST_lpdf.html` explains the derivation of such functions.  

For illustration purposes, we will use Sweden's counties whose geometries are stored in `Sweden`. 

- `poi_stdmarcar.R` and `poi_stdmarcar.stan` illustrates the usage of spatio-temporal effects following a MARCAR distribution.
- `poi_stdcar.R` and `poi_stdcar.stan` illustrates the usage of spatio-temporal effects following a CAR distribution.
