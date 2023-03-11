function [misfit_values, E_opt, K_opt, n_opt] = misfit_mincon(eps_exp, tn, beam, dsigma, E_range, K_range, n_range)
options = optimset('Display','off', 'Algorithm', 'interior-point', 'TolFun', 1e-6);
misfit_min = Inf;
misfit_values = zeros(length(E_range), length(K_range), length(n_range));

% Set up bounds for E, K, and n
lb = [min(E_range) min(K_range) min(n_range)];
ub = [max(E_range) max(K_range) max(n_range)];

% Set up multiple initial guesses for x0
x0_list = [E_range(1) K_range(1) n_range(1);           E_range(end) K_range(end) n_range(end);           mean(E_range) mean(K_range) mean(n_range)];

for i = 1:size(x0_list, 1)
    x0 = x0_list(i, :);
    [x, misfit_temp] = fmincon(@(x) misfit_sig(x, eps_exp, tn, beam, dsigma), x0, [], [], [], [], lb, ub, [], options);
    misfit_values(i) = misfit_temp;
    if misfit_temp < misfit_min
        misfit_min = misfit_temp;
        E_opt = x(1);
        K_opt = x(2);
        n_opt = x(3);
    end
end
end
