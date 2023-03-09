function [misfit_min, E_opt, K_opt, n_opt] = misfit_minunc(eps_exp, tn, beam, dsigma, E_range, K_range, n_range)
    % Define anonymous function to pass additional arguments to objective function
    obj_func = @(x) misfit_sig(x, eps_exp, tn, beam, dsigma);
    % Define initial guess for optimization variables
    x0 = [E_range(1); K_range(1); n_range(1)];
    % Define optimization options
    options = optimoptions('fminunc', 'Display', 'off');
    % Perform optimization
    [x_opt, misfit_min] = fminunc(obj_func, x0, options);
    % Assign optimal values to output variables
    E_opt = x_opt(1);
    K_opt = x_opt(2);
    n_opt = x_opt(3);
end
