# SCP_vonmises
Matlab code associated with sequential change-point procedures for von Mises data

The file D_stat.m calculates statistic (6) in Potgieter (2019) which is used to construct the 
sequential change-point location CUSUM. The file uses inputs (X,N,mu,kappa) where X = (X1,...,XN) 
is a sequence of angular observations, N is the length of the sequence, and (mu,kappa) are the 
in-control location and concentration parameters of the von Mises distribution.

The file calibrate_LRT_location.m can be used to calculate calibrated h-values (seqential limits) for
the sequential change-point CUSUM for a change in location. The file uses inputs (kappa,alpha,rep,K,J)
where kappa is the known concentration parameter of the von Mises distribution, alpha = 1/ARL is the 
inverse of the desired Average Run Length, rep is a parameter used to set the random seed (and is useful 
if multiple runs for the same kappa and alpha are implemented, ensuring these are based on different random
number sequences). Finally, K is the number of sequential calibrated h-values to compute and J is the
sample size for the Monte Carlo simulation.

Calibration was performed for alpha = {0.01,0.002,0.001} and for katppa = {0.5,1,2}. A total of 5 runs, 
each based on a sample size of J = 200,000 were calculated for K = 500. These are stored in the file
changept_calb_loc.mat.

The file h_rational_interpolation.m shows how to access the h-sequences in changept_calb_loc. Furthermore,
this file fits a rational linear function to the sequence (j,h_j), j = 1,...,J. The rational linear
function fits very well and the fitted values (rather than the raw Monte Carlo h-values) were used for
in Potgieter (2019) for j > 10.
