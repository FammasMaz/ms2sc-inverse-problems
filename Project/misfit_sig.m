function misfit_sig=misfit_sig(x, eps_exp, tn, beam, dsigma, lambda, arg)
    if nargin == 5
        lambda = 0;
        arg = 'squared';
    elseif nargin == 6
        arg = 'squared';
    elseif nargin == 7
        arg = arg;
    else
        error('Wrong number of input arguments')
    end

    % Calculation of the misfit function with no gradient estimation
    beamtemp.n = x(3); % update of the parameter
    beamtemp.E = x(1);
    beamtemp.K = x(2);
    %[t, epsilon]=ode45(@forward_sigma, tn, 0,[], beamtemp, dsigma);
    epsilon = forana(beamtemp, tn, dsigma);
    if isequal(arg, 'squared')
        misfit_sig = quad(@discrepancy,tn(1),tn(end),[],0,eps_exp,epsilon,x(3),tn);
    elseif isequal(arg, 'absolute')
        misfit_sig = mean(abs(epsilon - eps_exp));
    elseif isequal(arg, 'RMSE')
        misfit_sig = sqrt(mean((epsilon - eps_exp).^2));
    elseif isequal(arg, 'MSE')
        misfit_sig = mean((epsilon - eps_exp).^2);
    end


% Calculate RMSE/MSE

% Calculate regularization term

% Total objective function

