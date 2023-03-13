function misfit_sig=misfit_sig(x, eps_exp, tn, beam, dsigma, lambda)
if nargin == 5
    lambda = 0;
end
% Calculation of the misfit function with no gradient estimation
beamtemp.n = x(3); % update of the parameter
beamtemp.E = x(1);
beamtemp.K = x(2);
%[t, epsilon]=ode45(@forward_sigma, tn, 0,[], beamtemp, dsigma);
epsilon = forana(beamtemp, tn, dsigma);

% squared error integrated over time
%misfit_sig=quad(@discrepancy,tn(1),tn(end),[],0,eps_exp,epsilon,x(3),tn);

% Calculate RMSE/MSE
misfit_sig = sqrt(mean((epsilon - eps_exp).^2));

% Calculate regularization term
%reg_term = lambda * sqrt(sum(x.^2));

% Total objective function
% misfit_sig = misfit_sig + reg_term;

% Define regularization parameter
% Define singular value decomposition of the forward problem
% [U,S,V] = svd(forana(beamtemp, tn, dsigma));
% s = diag(S);

% % Define truncated singular value decomposition
% s_trunc = zeros(size(s));
% s_trunc(s > lambda) = s(s > lambda);
% S_trunc = diag(s_trunc);

% % Regularized misfit function calculation
% reg = mean((U*S_trunc*V' * eps_exp - eps_exp).^2);

% misfit_sig = misfit_sig + reg;
