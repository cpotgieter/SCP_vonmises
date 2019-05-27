function [h100,h500,h1000] = h_rational_interpolation(kappa)

% Calculates the h-sequence from fitting a rational polynomial model
% This was used in Potgieter (2019) rather than raw h-sequence scores
% obtained from a Monte Carlo simulation

% Sequences h100, h500, h1000 are calibrated to ARLs of 100, 500, 1000.

% kappa is the known in-control concentration parameter of the von Mises
% distribution

% Below keeps the first 10 raw scores and then uses the interpolants (the 
% interpolants don't fit well for small n) and the remaining values
% are from interpolation

% kappa values that were pre-calibrated
% i0 allows for numerical error 

kappa_in = [0.5,1,2];
i0 = find(abs(kappa_in-kappa)<=0.001);

load('changept_calb_loc.mat','h_store')
i = i0;
j = 1; %alpha = 0.01
H = h_store{i,j}';
f = fit((1:500)', H, 'rat22','StartPoint',[max(H),0,0,0,0]);
h100 = [H(1:10);f(11:500)];
j = 2; %alpha = 0.002
H = h_store{i,j}';
f = fit((1:500)', H, 'rat22','StartPoint',[max(H),0,0,0,0]);
h500 = [H(1:10);f(11:500)];
j = 3; %alpha = 0.001
H = h_store{i,j}';
f = fit((1:500)', H, 'rat22','StartPoint',[max(H),0,0,0,0]);
h1000 = [H(1:10);f(11:500)];