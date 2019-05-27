function h_seq = calibrate_LRT_location(kappa,alpha,rep,K,J)

% Code calculates sequence of control limits to monitor statistics D_{max,n}
% equation (7) in Potgieter (2019)

% kappa is the known concentration parameter of the von Mises distribution

% alpha = 1/ARL is the desired Average Run Length

% rep is a parameter that is helpful if several versions of the code is
% being run and different random number sequences are required in
% parallelization

% K is the number of consecutive h-values to compute
% Calculations slow down exponentially fast for increasing K

% J is the number of samples used to compute the h-values

% For illustration, set K = 100 and J = 10^5

% Below sets random seed
rng(1900*alpha*rep)

% Generate a sequence of J von Mises random variables under in-control
% conditions
mu0 = 0;
X0 = circ_vmrnd(mu0,kappa,J)';

% Compute J replicates of the statistic for n = 1
n = 1;
for j = 1:J
    X = X0(j)';
    M(n,j) = D_stat(X,n,mu0,kappa);
end
M_new = M(1,:);
% Finds the estimated control limit for n = 1
h = quantile(M_new,1-alpha);
X0(:,(M_new>h))=[];
M(:,(M_new>h))=inf;

% Simulate new von Mises numbers satisfying the condition that [M_new<h]
% back to having J replicates
[~,m] = size(X0);
while m<J
    X_new = circ_vmrnd(mu0, kappa, n);
    M_new = D_stat(X,n,mu0,kappa);
    if M_new<h
        X0 = [X0,X_new];
        [~,m] = size(X0);
    end
end

% Repeat above procedure for all sample sizes up to K-1
for n = 2:(K-1)
    X0 = [X0;circ_vmrnd(mu0,kappa,J)'];
    for j = 1:J
        X = X0(:,j)';
        M_out = D_stat(X,n,mu0,kappa);
        M(n,j) = M_out;
    end
    M_new = M(n,:);
    h(n,1) = quantile(M_new,1-alpha);
    f = find(M_new>h(n));
    X0(:,f)=[];
    M(n,f)=inf;
    [~,m] = size(X0);
    while m<J
        X_new = circ_vmrnd(mu0, kappa, n);
        M_new = D_stat_all(X_new',n,mu0,kappa);
        if sum(M_new<h)<n
            X0 = [X0,X_new];
            [~,m] = size(X0);
        end
    end

end

% Repeat simulation for n = K
% (This step does not require simulation of new random numbers)
n = K;
X0 = [X0;circ_vmrnd(mu0,kappa,J)'];
for j = 1:J
    X = X0(:,j);
    M_out = D_stat(X',n,mu0,kappa);
    M(n,j) = M_out;
end 
M_new = M(n,:);
h(n,1) = quantile(M_new,1-alpha);
%f = find(M_new>h(n));
%X0(:,f)=[];
%M(n,f)=inf;

h_seq = h;

end