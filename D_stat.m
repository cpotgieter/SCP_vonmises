function D = D_stat(X,N,mu,kappa)
% X = [X1,X2,...,XN] is a random sequence defined on a circle
% N = length(X)
% mu and kappa are the location and concentration parameters of the von
% Mises distribution under in-control conditions
% Typically this means mu = 0

% Code computes statistic D_{max,N} based on equation (6) in Potgieter
% (2019)

S = sin(X-mu); C = cos(X-mu);
cs_C = cumsum(C);
cs_S = cumsum(S);

delta_kn = atan2((cs_S(N)-[0,cs_S(1:(N-1))]),(cs_C(N)-[0,cs_C(1:(N-1))]));
D = kappa*max(abs((cos(delta_kn)-1).*(cs_C(N)-[0,cs_C(1:(N-1))])+sin(delta_kn).*(cs_S(N)-[0,cs_S(1:(N-1))])));

end