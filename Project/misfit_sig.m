function misfit_sig=misfit_sig(x, eps_exp, tn, beam, dsigma)
% Calculation of the misfit function with no gradient estimation
beamtemp = beam;
beamtemp.n = x(3); % update of the parameter
beamtemp.E = x(1);
beamtemp.K = x(2);
%[t, epsilon]=ode45(@forward_sigma, tn, 0,[], beamtemp, dsigma);
t = tn;
epsilon = forana(beamtemp, t, dsigma);

% squared error integrated over time
%misfit_sig=quad(@discrepancy,tn(1),tn(end),[],0,eps_exp,epsilon,x(3),tn);

% Calculate RMSE/MSE
misfit_sig = sqrt(mean((epsilon - eps_exp).^2));
