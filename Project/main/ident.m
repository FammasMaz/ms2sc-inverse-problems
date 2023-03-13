function [misfit_sig, E_opt, K_opt, n_opt] = ident(x0, eps_exp, tn, beam, dsigma)
% Calculation of the misfit function with no gradient estimation
options = optimoptions('fminunc', 'Display', 'off');
[x_opt, misfit_sig] = fminunc(@(x) misfit_fun(x, eps_exp, tn, beam, dsigma), x0, options);
E_opt = x_opt(1);
K_opt = x_opt(2);
n_opt = x_opt(3);
end

function misfit = misfit_fun(x, eps_exp, tn, beam, dsigma)
% Calculates the misfit function for given parameter values
beamtemp = beam;
beamtemp.E = x(1);
beamtemp.K = x(2);
beamtemp.n = x(3);
[t, epsilon] = ode45(@forward_sigma, tn, 0, [], beamtemp, dsigma);
misfit = sqrt(mean((epsilon - eps_exp).^2));
end
