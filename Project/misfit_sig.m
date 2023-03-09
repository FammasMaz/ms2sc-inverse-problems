function misfit_sig=misfit_sig(n, eps_exp, tn, beam, dsigma)
% Calculation of the misfit function with no gradient estimation
beamtemp = beam;
beamtemp.n = n; % update of the parameter
[t, epsilon]=ode45(@forward_sigma, tn, 0,[], beamtemp, dsigma); % simulation of the response using the parameter g
misfit_sig=quad(@discrepancy,tn(1),tn(end),[],0,eps_exp,epsilon,n,tn); % calculation of the integral over time of the discrepancy model-exp