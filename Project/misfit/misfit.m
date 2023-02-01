function misfit=misfit(n,strainexp,tn, K, s, E, sigma_dot)
% Calculation of the misfit function with no gradient estimation
[t, strain]=ode45(@forward,[0:.01:1],0,[],K, n, s, E, sigma_dot); % creation of synthetic data% simulation of the response using the parameter g
misfit=quad(@discrepancy,tn(1),tn(end),[],0,strainexp, strain,n,tn); % calculation of the integral over time of the discrepancy model-exp